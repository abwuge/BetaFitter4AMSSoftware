#include "Util.hh"
#include "amschain.h"
#include <fstream>

bool Util::addInputFile(std::string &inputFile, AMSChain &chain)
{
    try
    {
        if (inputFile.substr(inputFile.find_last_of(".") + 1) == "root")
        {
            chain.Add(inputFile.c_str());
            return true;
        }

        std::ifstream file(inputFile);
        if (!file.is_open())
            return false;

        std::string line;
        while (std::getline(file, line))
        {
            if (!line.empty())
            {
                chain.Add(line.c_str());
            }
        }

        return true;
    }
    catch (...)
    {
        return false;
    }
}