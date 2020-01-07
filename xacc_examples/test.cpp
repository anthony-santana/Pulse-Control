// !!! TEMP CODE !!! Thien Nguyen: For testing purposes only

#include "xacc.hpp"

int main (int argc, char** argv) {

	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    double pulseAmp = 0.0;
    if (argc <= 1) 
    {
        // No drive pulse amplitude provided in the command line, it will be default as 0.0 (mo drive)
        std::cerr << "Usage: " << argv[0] << " --drive-amp VALUE\n"
            << "\t--drive-amp\tSet the amplitude of the Gaussian Pulse drive.\n"
            << "\tIf not provided, it will be set to 0.0."
            << std::endl;
    }
    else
    {
        for (int i = 1; i < argc; ++i) 
        {
            const std::string arg = argv[i];
            if (arg == "--drive-amp") 
            { 
                if (i + 1 < argc) 
                { 
                    const std::string ampAsString = argv[++i]; 
                    size_t doubleStringSize; 
                    try 
                    {
                        pulseAmp = std::stod(ampAsString, &doubleStringSize);
                    }
                    catch (...) 
                    {
                       // Invalid number
                       std::cerr << "'" << ampAsString << "' cannot be converted to a number." << std::endl;
                       return 1;
                    }
                } 
                else 
                {   // No value is provided
                    std::cerr << "--drive-amp option requires one argument." << std::endl;
                    return 1;
                }  
            } 
        }
    }
    
    // Get the Pulse simulator
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("sim-mode", "Pulse"), std::make_pair("drive_amp", pulseAmp) });    
    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, nullptr);
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}