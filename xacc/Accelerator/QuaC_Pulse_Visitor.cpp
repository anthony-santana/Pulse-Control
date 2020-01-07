#include "QuaC_Pulse_Visitor.hpp"
#include "xacc.hpp"

extern "C" {
#include "interface_xacc_ir.h"
}

namespace {
   double gaussian_scaling;
   double gaussian_mu;
   double gaussian_sigma;
   double omega_lo;
   class GaussianFunc {
      public:  
      static double value(double in_time) 
      {
         TODO(This implementation is not scalable. Need to handle multiple LOs and frame changes. This will involve QuaC changes.)
         // f(t) = scaling * exp(-(t-mu)^2/(2 * sigma^2)) * cos (omega_LO * t)
         // (drive pulse mixed with LO)
         const double result = gaussian_scaling * exp(-pow(in_time - gaussian_mu, 2)/(2.0 * pow(gaussian_sigma, 2))) * cos(omega_lo * in_time);
         return result;
      }  
   };
}

// Macro to declare *global* params for the Gaussian pulse function
#define DECLARE_GAUSSIAN_PULSE(scaling, mu, sigma, omega) {\
  gaussian_scaling = scaling;\
  gaussian_mu = mu;\
  gaussian_sigma = sigma;\
  omega_lo = omega;\
}

// Get the Gaussian Pulse function (depends on time only)
// its params must be declared by DECLARE_GAUSSIAN_PULSE before use.
// otherwise, it will not be the Gaussian function that you may expect!!!
#define GAUSSIAN_PULSE_FUNC GaussianFunc::value
namespace QuaC {
   void PulseVisitor::initialize(std::shared_ptr<AcceleratorBuffer> buffer, const HeterogeneousMap& in_params) 
   {
      // Debug
      std::cout << "Initialize Pulse simulator \n";
      
      // Gaussian pulse params
      double driveSigma = 75;
      double driveMu = driveSigma*4;
      
      // Initialize some params for testing
      TODO(Get rid of these default params and enforce that they must be set upstream.)
      double dt = 0.1;
      double stopTime = driveSigma*8;
      int stepMax = 10000;
      double nu = 5.0;
      double omega = 2 * M_PI * nu;
      double pulseAmplitude = 1.0;
      // Qubit decay: just use a very small value
      double kappa = 0.000001;
      // Get the params from in_params map
      // TODO: this is a skeleton only atm,
      // i.e. I only retrieve the amplitude from the input HeterogeneousMap (to test Rabi osc)
      // we should pass/retrieve every param from this map (e.g. pulse shapes, pulse params, simulator configs, etc.)
      if (in_params.keyExists<double>("drive_amp")) 
      {
         const double driveAmp = in_params.get<double>("drive_amp");
         if (driveAmp < 0.0)
         {
            xacc::error("Invalid drive signal amplitude!\n");
            return;
         }
         
         // Set the pulse amplitude
         pulseAmplitude = driveAmp;
         // Debug
         std::cout << "Pulse simulator: Set pulse amplitude to " << pulseAmplitude << ". \n";
      }

      XACC_QuaC_InitializePulseSim(buffer->size(), dt, stopTime, stepMax);
      
      // Debug:
      // XACC_QuaC_SetLogVerbosity(DEBUG_DIAG);
      
      // TODO: Parse the QObj to get the Hamiltonian params
      // For now, just hard-coded
      // H = 2*pi*nu(1-sigma_z)/2 + D(t)*sigma_x
      // D(t) is a Gaussian pulse 
      
      // Time-independent terms:
      XACC_QuaC_AddConstHamiltonianTerm1("I", 0, { omega/2.0, 0.0});
      XACC_QuaC_AddConstHamiltonianTerm1("Z", 0, { -omega/2.0, 0.0});

      // Time-dependent term:
      DECLARE_GAUSSIAN_PULSE(pulseAmplitude, driveMu, driveSigma, omega);
      XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, GAUSSIAN_PULSE_FUNC);

      XACC_QuaC_AddQubitDecay(0, kappa);
   }

   void PulseVisitor::solve() 
   {
      std::cout << "Pulse simulator: solving the Hamiltonian. \n";
      double* results = nullptr;
      const auto resultSize = XACC_QuaC_RunPulseSim(&results);
      std::cout << "Final result: ";
      for (int i = 0; i < resultSize; ++i)
      {
         std::cout << results[i] << ", ";
      }
      std::cout << "\n";

      free(results);
   }

   void PulseVisitor::finalize() 
   {     
      std::cout << "Pulse simulator: Finalized. \n";
      XACC_QuaC_Finalize();
   }
}