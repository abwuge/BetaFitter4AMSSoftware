#include <iostream>
#include <string>

#include "amschain.h"

#include "Util.hh"
#include "DataProcessor.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <inputFile> <outputFile.root> <maxEventsNumber>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int maxEventsNumber = argc > 3 ? std::stoi(argv[3]) : -1;

    AMSChain chain;
    if (!Util::addInputFile(inputFile, chain))
    {
        std::cerr << "Error: could not add input file" << std::endl;
        return 1;
    }

    DataProcessor processor(outputFile);
    processor.processEvents(chain, maxEventsNumber);

    return 0;
}