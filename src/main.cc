#include <iostream>
#include <string>

#include "Util.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <input.root|input.list> <output.root> [<Option>] [<Energy Loss Scale>] [<Tracker Energy Loss Scale>] [<Scale Mode: all|s1s2|s2>] [<Reference Point: center|before_tof>]" << std::endl;
        std::cout << "Option: " << std::endl;
        std::cout << "  -2: Save energy loss information to ROOT file" << std::endl;
        std::cout << "  -1: Save magnetic field information to ROOT file" << std::endl;
        std::cout << "   0: Save beta information to ROOT file" << std::endl;
        std::cout << "   1: Save per-event energy loss scale information to ROOT file" << std::endl;
        std::cout << "   2: Benchmark BetaNL::Beta() function execution time" << std::endl;
        std::cout << "   3: Save beta difference information to ROOT file" << std::endl;
        std::cout << "   4: Jointly fit global TOF and tracker energy-loss scales" << std::endl;
        std::cout << "      Option 4 arguments: [betaMax=0.9] [scale mode] [reference point] [zetaMin=0] [zetaMax=6] [trackerScaleMin=0] [trackerScaleMax=20]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int Option = argc > 3 ? atoi(argv[3]) : 0;
    double energyLossScale = 2;
    double trackerEnergyLossScale = 1;
    if (Option != 4)
    {
        energyLossScale = argc > 4 ? atof(argv[4]) : 2;
        trackerEnergyLossScale = argc > 5 ? atof(argv[5]) : 1;
    }
    const double globalBetaMax = argc > 4 ? atof(argv[4]) : 0.9;
    const double globalZetaMin = Option == 4 && argc > 7 ? atof(argv[7]) : 0.0;
    const double globalZetaMax = Option == 4 && argc > 8 ? atof(argv[8]) : 6.0;
    const double trackerScaleMin = Option == 4 && argc > 9 ? atof(argv[9]) : 0.0;
    const double trackerScaleMax = Option == 4 && argc > 10 ? atof(argv[10]) : 20.0;
    const std::string energyLossScaleModeName =
        Option == 4 ? (argc > 5 ? argv[5] : "all") : (argc > 6 ? argv[6] : "all");
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

    const std::string referencePointName =
        Option == 4 ? (argc > 6 ? argv[6] : "center") : (argc > 7 ? argv[7] : "center");
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
        info += "Energy Loss Scale: " + std::to_string(energyLossScale) +
                "\nTracker Energy Loss Scale: " + std::to_string(trackerEnergyLossScale) + "\n";
    else if (Option == 4)
        info += "Global Beta Maximum: " + std::to_string(globalBetaMax) +
                "\nGlobal Zeta Range: [" + std::to_string(globalZetaMin) +
                ", " + std::to_string(globalZetaMax) + "]" +
                "\nGlobal Tracker Scale Range: [" + std::to_string(trackerScaleMin) +
                ", " + std::to_string(trackerScaleMax) + "]\n";

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
                              referencePoint)
                   ? 0
                   : 1;
    case 1:
        return Util::saveEnergyLossScale(inputFile, outputFile, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 2:
        return Util::benchmarkBetaNL(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 3:
        return Util::saveBetaDiff(inputFile, outputFile, energyLossScale, energyLossScaleMode, referencePoint) ? 0 : 1;
    case 4:
        return Util::saveGlobalEnergyLossScale(inputFile, outputFile, globalBetaMax,
                                               globalZetaMin, globalZetaMax,
                                               energyLossScaleMode, referencePoint,
                                               trackerScaleMin, trackerScaleMax)
                   ? 0
                   : 1;
    default:
        Util::test();
        return 1;
    }
}
