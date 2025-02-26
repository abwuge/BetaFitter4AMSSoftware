#ifndef __UTIL_HH__
#define __UTIL_HH__

#include <string>

#include "amschain.h"

/**
 * @brief Utility functions
 */
namespace Util
{
    /**
     * @brief Add input file to the chain
     *
     * @param inputFile Input file
     * @param chain Chain
     * @return true if the file was added successfully
     * @return false if the file was not added successfully
     */
    bool addInputFile(std::string &inputFile, AMSChain &chain);
}

#endif // __UTIL_HH__