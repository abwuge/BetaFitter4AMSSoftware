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
     * @brief Load particle data from amstreea in ROOT file
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @return std::vector<ParticleData> Vector of particle data, empty if loading fails
     */
    std::vector<ParticleData> loadParticleData(const std::string &inputFile);

    /**
     * @brief Get particle mass from PDG ID
     *
     * Supports Geant3 particle IDs and Geant4 PDG IDs if no Geant3 counterpart is found
     *
     * @param pdgId Geant3 particle ID or Geant4 PDG ID
     * @return float Mass in GeV/c^2, 0 if particle not recognized
     */
    float getMassFromPDG(int pdgId, double charge);
}

#endif // __UTIL_HH__