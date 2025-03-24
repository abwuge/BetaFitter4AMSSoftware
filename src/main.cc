#include <iostream>

#include "Util.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <inputFile.root> <outputFile.root> [<Option>] [<Energy Loss Scale>]" << std::endl;
        std::cout << "Option: " << std::endl;
        std::cout << "  -2: Save energy loss information to ROOT file" << std::endl;
        std::cout << "  -1: Save magnetic field information to ROOT file" << std::endl;
        std::cout << "   0: Save beta information to ROOT file" << std::endl;
        std::cout << "   1: Save energy loss scacle information to ROOT file" << std::endl;
        std::cout << "   2: Benchmark BetaNL::Beta() function execution time" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int Option = argc > 3 ? atoi(argv[3]) : 0;
    double energyLossScale = argc > 4 ? atof(argv[4]) : 2;

    std::string info = "\nInput file: " + inputFile + "\nOutput file: " + outputFile + "\nOption: " + std::to_string(Option) + "\n";

    if (Option == 0)
        info += "Energy Loss Scale: " + std::to_string(energyLossScale) + "\n";

    std::cout << info << std::endl;

    switch (Option)
    {
    case -2:
        return Util::saveEnergyLoss(inputFile, outputFile) ? 0 : 1;
    case -1:
        return Util::saveMagneticField(outputFile) ? 0 : 1;
    case 0:
        return Util::saveBeta(inputFile, outputFile, energyLossScale) ? 0 : 1;
    case 1:
        return Util::saveEnergyLossScale(inputFile, outputFile) ? 0 : 1;
    case 2:
        return Util::benchmarkBetaNL(inputFile, outputFile, energyLossScale) ? 0 : 1;
    default:
        std::cerr << "Error: Invalid option " << Option << std::endl;
        return 1;
    }
}