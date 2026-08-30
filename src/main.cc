#include <iostream>
#include <string>

#include "Util.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <input.root> <output.root> [<Option>] [<Energy Loss Scale>] [<Tracker Energy Loss Scale>] [<Scale Mode: all|s1s2|s2>] [<Reference Point: center|before_tof>] [<Scale Config: none|path>]" << std::endl;
        std::cout << "Option: " << std::endl;
        std::cout << "  -2: Save energy loss information to ROOT file" << std::endl;
        std::cout << "  -1: Save magnetic field information to ROOT file" << std::endl;
        std::cout << "   0: Save beta information to ROOT file" << std::endl;
        std::cout << "   1: Fit per-event TOF and tracker energy-loss scales" << std::endl;
        std::cout << "   2: Benchmark BetaNL::Beta() function execution time" << std::endl;
        std::cout << "   3: Save beta difference information to ROOT file" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int Option = argc > 3 ? atoi(argv[3]) : 0;
    const double energyLossScale = argc > 4 ? atof(argv[4]) : 2;
    const double trackerEnergyLossScale = argc > 5 ? atof(argv[5]) : 1;
    const std::string energyLossScaleModeName = argc > 6 ? argv[6] : "all";
    EnergyLossScaleMode energyLossScaleMode;
    if (energyLossScaleModeName == "all")
        energyLossScaleMode = EnergyLossScaleMode::All;
    else if (energyLossScaleModeName == "s1s2")
        energyLossScaleMode = EnergyLossScaleMode::S1S2;
    else if (energyLossScaleModeName == "s2")
        energyLossScaleMode = EnergyLossScaleMode::S2Only;
    else
    {
        std::cerr << "Error: Scale mode must be 'all', 's1s2', or 's2'" << std::endl;
        return 1;
    }

    const std::string referencePointName = argc > 7 ? argv[7] : "center";
    BetaReferencePoint referencePoint;
    if (referencePointName == "center")
        referencePoint = BetaReferencePoint::AMSCenter;
    else if (referencePointName == "before_tof")
        referencePoint = BetaReferencePoint::BeforeTOF;
    else
    {
        std::cerr << "Error: Reference point must be 'center' or 'before_tof'" << std::endl;
        return 1;
    }

    const std::string scaleConfig = argc > 8 ? argv[8] : "none";

    std::string info = "\nInput file: " + inputFile + "\nOutput file: " + outputFile + "\nOption: " + std::to_string(Option) +
                       "\nEnergy Loss Scale Mode: " + energyLossScaleModeName +
                       "\nBeta Reference Point: " + referencePointName +
                       "\nScale Config: " + scaleConfig + "\n";

    if (Option == 0)
        info += "Energy Loss Scale: " + std::to_string(energyLossScale) +
                "\nTracker Energy Loss Scale: " + std::to_string(trackerEnergyLossScale) + "\n";

    std::cout << info << std::endl;

    switch (Option)
    {
    case -2:
        return Util::saveEnergyLoss(inputFile, outputFile) ? 0 : 1;
    case -1:
        return Util::saveMagneticField(outputFile) ? 0 : 1;
    case 0:
        return Util::saveBeta(inputFile, outputFile, energyLossScale,
                              trackerEnergyLossScale, energyLossScaleMode,
                              referencePoint, scaleConfig)
                   ? 0
                   : 1;
    case 1:
        return Util::saveEnergyLossScale(inputFile, outputFile, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 2:
        return Util::benchmarkBetaNL(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 3:
        return Util::saveBetaDiff(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    default:
        Util::test();
        return 1;
    }
}
