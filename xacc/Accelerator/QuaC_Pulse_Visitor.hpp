#pragma one

#include "Identifiable.hpp"
#include "OptionsProvider.hpp"
#include "AcceleratorBuffer.hpp"
#include "Cloneable.hpp"

using namespace xacc;

namespace QuaC {    
    class PulseVisitor : public OptionsProvider, public Cloneable<PulseVisitor>
    {
    public:
        void initialize(std::shared_ptr<AcceleratorBuffer> buffer, const HeterogeneousMap& in_params = {});
        void solve();
        void finalize();
        virtual std::shared_ptr<PulseVisitor> clone() { return std::make_shared<PulseVisitor>(); }
    };    
}