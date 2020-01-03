#include "cppmicroservices/BundleActivator.h"
#include "cppmicroservices/BundleContext.h"
#include "cppmicroservices/ServiceProperties.h"
#include "QuaC_Accelerator.hpp"
#include <iostream>

using namespace cppmicroservices;

class US_ABI_LOCAL QuaC_Activator : public BundleActivator {
public:
  QuaC_Activator() {}

  void Start(BundleContext context) {
    // TODO
    std::cout << ">> DEBUG: Starting QuaC bundle ...\n";
    auto acc = std::make_shared<QuaC::QuaC_Accelerator>();
    // Register the QuaC accelerator with the service registry
    context.RegisterService<xacc::Accelerator>(acc);
  }

  void Stop(BundleContext context) {}
};

CPPMICROSERVICES_EXPORT_BUNDLE_ACTIVATOR(QuaC_Activator)
