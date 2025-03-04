#ifndef __UTIL_HH__
#define __UTIL_HH__

#include <string>
#include <vector>
#include "amschain.h"
#include "ParticleData.hh"
#include <TCanvas.h>
#include <TGraph.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TAxis.h>

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

    /**
     * @brief Load particle data from ROOT file
     * 
     * @param inputFile Path to the ROOT file containing particle data
     * @return std::vector<ParticleData> Vector of particle data, empty if loading fails
     */
    std::vector<ParticleData> loadParticleData(const std::string &inputFile);

    /**
     * @brief Draw particle trajectory and TOF hits in 3D
     * 
     * @param particle Particle data containing MC truth and TOF hits
     * @param outputPath Path to save the output plot (PDF format)
     * @param nPoints Number of points to plot on the trajectory (default: 100)
     * @return true if drawing was successful
     * @return false if drawing failed
     */
    bool drawTrajectory(const ParticleData &particle, const std::string &outputPath, int nPoints = 100);

    /**
     * @brief Load particle data from amstreea in ROOT file
     * 
     * @param inputFile Path to the ROOT file containing MC particle data
     * @return std::vector<ParticleData> Vector of particle data, empty if loading fails
     */
    std::vector<ParticleData> loadParticleData2(const std::string &inputFile);
}

#endif // __UTIL_HH__