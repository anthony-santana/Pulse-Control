#pragma once
#include <sstream>
#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/functional.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/utility.hpp>

using SerializationType = std::stringstream;
using SerializationOutputDataType =  cereal::BinaryOutputArchive;
using SerializationInputDataType = cereal::BinaryInputArchive;

#define DECLARE_CORE_TYPE(parentType, derivedType)\
    CEREAL_REGISTER_TYPE(derivedType)\
    CEREAL_REGISTER_POLYMORPHIC_RELATION(parentType, derivedType)

template<typename T>    
void deserializeBuffer(char* in_buffer, int in_length, T& out_typedVar)
{
    SerializationType ss(std::string(in_buffer, in_buffer + in_length));
    SerializationInputDataType inputSerialization(ss); 
    inputSerialization(out_typedVar);
}

template<typename T>    
std::string serializeObject(const T& in_obj)
{
    SerializationType ss;
    SerializationOutputDataType outputSerializer(ss); 
    outputSerializer(in_obj); 
    return ss.str();
}