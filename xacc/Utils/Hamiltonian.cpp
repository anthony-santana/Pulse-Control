#include "Hamiltonian.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include "xacc_service.hpp"
#include "expression_parsing_util.hpp"
#include <cassert>
#include "xacc.hpp"
#include "json.hpp"
extern "C" {
#include "interface_xacc_ir.h"
}

namespace {
    std::string removeWhiteSpaces(const std::string& in_str)
    {
        auto result = in_str;
        std::remove_if(result.begin(), result.end(), [](unsigned char c){ return std::isspace(c); });

        return result;
    }

    std::string toUpperCase(const std::string& in_str)
    {
        auto result = in_str;
        transform(result.begin(), result.end(), result.begin(), ::toupper); 
        return result;
    }

    bool isNumberString(const std::string& in_str)
    {
        return std::find_if(in_str.begin(), in_str.end(), [](const auto in_char){ return in_char != '0' && in_char != '1' &&
                                                                                         in_char != '2' && in_char != '3' &&
                                                                                         in_char != '4' && in_char != '5' &&
                                                                                         in_char != '6' && in_char != '7' &&
                                                                                         in_char != '8' && in_char != '9'; }) == in_str.end();
    }

    bool GetLastOperator(const std::string& in_str, QuaC::QubitOp& out_Op, std::string& out_subString)
    {
        // Find the last '*' character
        const auto pos = in_str.find_last_of("*");
        if (pos == std::string::npos)
        {
            return false;
        }

        const auto opName = toUpperCase(in_str.substr(pos + 1));     
        const auto qubitIdxIter = std::find_if(opName.begin(), opName.end(), [](const auto in_char){ return !isalpha(in_char); } );
        
        if (qubitIdxIter == opName.end())
        {
            return false;
        }

        const auto opStr = opName.substr(0, std::distance(opName.begin(), qubitIdxIter));


        if (opStr.length() < 1 ||  QuaC::ConvertOperatorFromString(opStr) == QuaC::Operator::NA)
        {
            return false;
        }

        const auto qubitIdxStr = opName.substr(std::distance(opName.begin(), qubitIdxIter));
        if (!isNumberString(qubitIdxStr))
        {
            return false;
        }

        const auto op = QuaC::ConvertOperatorFromString(opStr);
        const size_t qIdx = std::stoi(qubitIdxStr);
        
        out_Op = std::make_pair(op, qIdx);
        out_subString = in_str.substr(0, pos);
        return true;
    }

    bool TryEvaluateExpression(const std::string& in_exprString, const QuaC::VarsMap& in_vars, double& out_result)
    {
        auto parsingUtil = xacc::getService<xacc::ExpressionParsingUtil>("exprtk");
        std::vector<std::string> varNames;
        std::vector<double> varVals;
        for (const auto& kv : in_vars)
        {
            varNames.emplace_back(kv.first);
            varVals.emplace_back(kv.second);
        }

        if (!parsingUtil->validExpression(in_exprString, varNames))
        {            
            return false;
        }
 
        double evaled = 0.0;
    
        if (!parsingUtil->evaluate(in_exprString, varNames, varVals, evaled))
        {
            return false;
        }

        out_result = evaled;
        return true;
    }
    
}

namespace QuaC { 
std::unique_ptr<HamiltonianTimeDependentTerm> HamiltonianTimeDependentTerm::fromString(const std::string& in_string, const VarsMap& in_vars) 
{
    auto exprStr = removeWhiteSpaces(in_string);
    // Find the special '||' channel separator
    auto separatorPos = exprStr.find("||");
    if (separatorPos == std::string::npos)
    {
        return nullptr;
    }

    const auto channelName = toUpperCase(exprStr.substr(separatorPos + 2));
    // Minimum length: 2
    if (channelName.length() < 2 || (channelName.front() != 'D' && channelName.front() != 'U') || (!isNumberString(channelName.substr(1))))
    {
        return nullptr;
    }
    
    const auto operatorExpression = exprStr.substr(0, separatorPos);

    std::string opConstExpr;
    QubitOp resultOp;
    if (!GetLastOperator(operatorExpression, resultOp, opConstExpr))
    {
       return nullptr;
    }
   
    double evaled = 0.0;
    if (!TryEvaluateExpression(opConstExpr, in_vars, evaled))
    {
        return nullptr;
    }
  
    return std::make_unique<HamiltonianTimeDependentTerm>(channelName, evaled, resultOp);
}

std::unique_ptr<HamiltonianTimeIndependentTerm> HamiltonianTimeIndependentTerm::fromString(const std::string& in_string, const VarsMap& in_vars)
{
    auto exprStr = removeWhiteSpaces(in_string);

    // Don't process time-dependent terms
    auto separatorPos = exprStr.find("||");
    if (separatorPos != std::string::npos)
    {
        return nullptr;
    }
    std::vector<QubitOp> operators;

    QubitOp tempOp;
    std::string remainderStr;
    std::string tempStr = exprStr;
    while (GetLastOperator(tempStr, tempOp, remainderStr))
    {
        operators.emplace_back(tempOp);
        tempStr = remainderStr;
    }
 
    double evaled = 0.0;
    if (!TryEvaluateExpression(remainderStr, in_vars, evaled))
    {
        return nullptr;
    }
    
    // Reverse the vector list since we were parsing operators from the back.
    std::reverse(operators.begin(), operators.end());
    
    return std::make_unique<HamiltonianTimeIndependentTerm>(evaled, operators);
}

std::unique_ptr<HamiltonianSumTerm> HamiltonianSumTerm::fromString(const std::string& in_string, const VarsMap& in_vars)
{
    static const std::string SUM_TERM_PREFIX = "_SUM[";
    auto exprStr = removeWhiteSpaces(in_string);
    if (exprStr.compare(0, SUM_TERM_PREFIX.size(), SUM_TERM_PREFIX) != 0 || exprStr.back() != ']')
    {
        // Not a Sum term
        return nullptr;
    }

    exprStr =  exprStr.substr(SUM_TERM_PREFIX.size());
    exprStr.pop_back();
    
    const auto GetNextCommaSeparatedSubString = [](std::string& io_string) -> std::string {
        auto commaPos = io_string.find(",");
        if (commaPos == std::string::npos)
        {
            return "";
        }
        
        const auto subStr = io_string.substr(0, commaPos);
        io_string = io_string.substr(commaPos + 1);
        return subStr;
    };
    
    
    const auto loopVarName = GetNextCommaSeparatedSubString(exprStr);
    const auto startValStr = GetNextCommaSeparatedSubString(exprStr);
    const auto endValStr = GetNextCommaSeparatedSubString(exprStr);
    const auto loopExpression = exprStr;
    const std::string varFmt = "{" + loopVarName + "}";
    
    if (loopVarName.empty() || startValStr.empty() || endValStr.empty() || loopExpression.empty() || 
            !isNumberString(startValStr) || !isNumberString(endValStr) ||  
            // The expression doesn't contain the loop index var!!!
            loopExpression.find(varFmt) == std::string::npos) 
    {
        return nullptr;
    }

    const int startLoopVal = std::stoi(startValStr);
    const int endLoopVal = std::stoi(endValStr);

    const auto resolveLoopTemplate = [](const std::string& in_string, const std::string& in_template, int in_val) -> std::string {
	    std::string result = in_string;
        // Get the first occurrence
	    size_t pos = in_string.find(in_template);
        const std::string replaceStr = std::to_string(in_val);
	    // Repeat till end is reached
	    while(pos != std::string::npos)
	    {
		    // Replace this occurrence 
		    result.replace(pos, in_template.size(), replaceStr);
		    // Get the next occurrence from the current position
		    pos = result.find(in_template, pos + replaceStr.size());
	    }

        return result;
    };


    if (startLoopVal > endLoopVal)
    {
        return nullptr;
    }

    const auto resolvedExpression = resolveLoopTemplate(loopExpression, varFmt, startLoopVal);
    const auto tryTimeIndependent = HamiltonianTimeIndependentTerm::fromString(resolvedExpression, in_vars);
    const auto tryTimeDependent = HamiltonianTimeDependentTerm::fromString(resolvedExpression, in_vars);

    if (tryTimeDependent == nullptr && tryTimeDependent == nullptr)
    {
        return nullptr;
    }

    const bool isTimeDependent = tryTimeDependent != nullptr;
    const auto parseLoopExpression = [&](const std::string& in_exprStr) ->  std::unique_ptr<HamiltonianTerm> {
        if (isTimeDependent)
        {
            auto result = HamiltonianTimeDependentTerm::fromString(resolvedExpression, in_vars);
            return std::unique_ptr<HamiltonianTerm>(result.release()); 
        }
        else
        {
            auto result = HamiltonianTimeIndependentTerm::fromString(resolvedExpression, in_vars);
            return std::unique_ptr<HamiltonianTerm>(result.release()); 
        }        
    };

    // Note: IBM uses an inclusive loop index (i.e. the end value is included)
    std::vector<std::unique_ptr<HamiltonianTerm>> loopOps;

    for (int i = startLoopVal; i <= endLoopVal; ++i)
    {
        const auto resolvedExpression = resolveLoopTemplate(loopExpression, varFmt, i);
        std::unique_ptr<HamiltonianTerm> result = parseLoopExpression(resolvedExpression);
        assert(result != nullptr);
        loopOps.emplace_back(std::move(result));
    }

    return std::make_unique<HamiltonianSumTerm>(std::move(loopOps), isTimeDependent);
}

void HamiltonianTimeIndependentTerm::apply(IChannelNameResolver* in_channelResolver)
{
    // This constraint can be lifted if necessary, just add API's to the backend. 
    if (m_operators.size() > 2)
    {
        xacc::error("We only support Hamiltonian terms which are products of maximum two operators.");
    }
    
    if (m_operators.size() == 1)
    {
        const auto op = m_operators.front();        
        XACC_QuaC_AddConstHamiltonianTerm1(OperatorToString(op.first).c_str(), op.second, { m_coefficient.real(), m_coefficient.imag()});    
    }
    else if (m_operators.size() == 2)
    {
        const auto op1 = m_operators[0];
        const auto op2 = m_operators[1];
        XACC_QuaC_AddConstHamiltonianTerm2(OperatorToString(op1.first).c_str(), op1.second, OperatorToString(op2.first).c_str(), op2.second, { m_coefficient.real(), m_coefficient.imag() });
    }
}

void HamiltonianTimeDependentTerm::apply(IChannelNameResolver* in_channelResolver)
{
    XACC_QuaC_AddTimeDependentHamiltonianTerm1(OperatorToString(m_operator.first).c_str(), m_operator.second, in_channelResolver->GetChannelId(m_channelName), m_coefficient);
}

void HamiltonianSumTerm::apply(IChannelNameResolver* in_channelResolver)
{
    for (auto& term : m_terms)
    {
        term->apply(in_channelResolver);
    }
}

std::unique_ptr<HamiltonianTerm> HamiltonianParsingUtil::tryParse(const std::string& in_expr, const VarsMap& in_vars)
{
    {
        auto trySum = HamiltonianSumTerm::fromString(in_expr, in_vars);
        if (trySum)
        {
            return std::unique_ptr<HamiltonianTerm>(trySum.release()); 
        }
    }
    
    {
        auto tryTimeDep = HamiltonianTimeDependentTerm::fromString(in_expr, in_vars);
    
        if (tryTimeDep)
        {
            return std::unique_ptr<HamiltonianTerm>(tryTimeDep.release()); 

        }
    }

    {
        auto tryTimeInd = HamiltonianTimeIndependentTerm::fromString(in_expr, in_vars);
    
        if (tryTimeInd)
        {
            return std::unique_ptr<HamiltonianTerm>(tryTimeInd.release()); 

        }
    }
    
    xacc::warning("Cannot parse Hamiltonian string " + in_expr);
    return nullptr;
}

bool HamiltonianParsingUtil::tryParse(const std::string& in_jsonString, std::function<void(HamiltonianTerm&)> in_forEachTermFn)
{
    auto j = nlohmann::json::parse(in_jsonString);
    if (!j.is_object())
    {
        xacc::warning("Hamiltonian JSON must be an object.");
        return false;
    }

    // Hamiltonian strings and vars map
    auto hamStrArray = j["h_str"];
    auto varsArray = j["vars"];
    
    VarsMap vars;
    for (auto varsArrayIter = varsArray.begin(); varsArrayIter != varsArray.end(); ++varsArrayIter)
    {
        vars.emplace(varsArrayIter.key(), varsArrayIter.value().get<double>());
    }

    for (auto hamStrIter = hamStrArray.begin(); hamStrIter != hamStrArray.end(); ++hamStrIter) 
    {
        auto hamStr = (*hamStrIter).get<std::string>();
        // std::cout << "Hamiltonian term: " << hamStr << "\n";
        auto parseResult = tryParse(hamStr, vars);
        if (parseResult == nullptr)
        {
            return false;
        }

        in_forEachTermFn(*parseResult);
    }

    return true;
}
}