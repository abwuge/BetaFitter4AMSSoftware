#include <iostream>
#include <string>

#include "Util.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <input.root|input.list> <output.root> [<Option>] [<Energy Loss Scale>] [<Scale Mode: all|s1s2|s2>] [<Reference Point: center|before_tof>]" << std::endl;
        std::cout << "Option: " << std::endl;
        std::cout << "  -2: Save energy loss information to ROOT file" << std::endl;
        std::cout << "  -1: Save magnetic field information to ROOT file" << std::endl;
        std::cout << "   0: Save beta information to ROOT file" << std::endl;
        std::cout << "   1: Save per-event energy loss scale information to ROOT file" << std::endl;
        std::cout << "   2: Benchmark BetaNL::Beta() function execution time" << std::endl;
        std::cout << "   3: Save beta difference information to ROOT file" << std::endl;
        std::cout << "   4: Fit one global energy loss scale shared by all MC events" << std::endl;
        std::cout << "      Option 4 arguments: [betaMax=0.9] [scale mode] [reference point] [zetaMin=0] [zetaMax=6]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int Option = argc > 3 ? atoi(argv[3]) : 0;
    double energyLossScale = argc > 4 ? atof(argv[4]) : 2;
    const double globalBetaMax = argc > 4 ? atof(argv[4]) : 0.9;
    const double globalZetaMin = argc > 7 ? atof(argv[7]) : 0.0;
    const double globalZetaMax = argc > 8 ? atof(argv[8]) : 6.0;
    const std::string energyLossScaleModeName = argc > 5 ? argv[5] : "all";
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

    const std::string referencePointName = argc > 6 ? argv[6] : "center";
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

    std::string info = "\nInput file: " + inputFile + "\nOutput file: " + outputFile + "\nOption: " + std::to_string(Option) +
                       "\nEnergy Loss Scale Mode: " + energyLossScaleModeName +
                       "\nBeta Reference Point: " + referencePointName + "\n";

    if (Option == 0)
        info += "Energy Loss Scale: " + std::to_string(energyLossScale) + "\n";
    else if (Option == 4)
        info += "Global Beta Maximum: " + std::to_string(globalBetaMax) +
                "\nGlobal Zeta Range: [" + std::to_string(globalZetaMin) +
                ", " + std::to_string(globalZetaMax) + "]\n";

    std::cout << info << std::endl;

    switch (Option)
    {
    case -2:
        return Util::saveEnergyLoss(inputFile, outputFile) ? 0 : 1;
    case -1:
        return Util::saveMagneticField(outputFile) ? 0 : 1;
    case 0:
        return Util::saveBeta(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 1:
        return Util::saveEnergyLossScale(inputFile, outputFile, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 2:
        return Util::benchmarkBetaNL(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 3:
        return Util::saveBetaDiff(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 4:
        return Util::saveGlobalEnergyLossScale(inputFile, outputFile, globalBetaMax,
                                               globalZetaMin, globalZetaMax,
                                               energyLossScaleMode, referencePoint)
                   ? 0
                   : 1;
    default:
        Util::test();
        return 1;
    }
}
