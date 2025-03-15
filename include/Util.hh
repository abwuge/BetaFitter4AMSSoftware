#ifndef __UTIL_HH__
#define __UTIL_HH__

#include <string>
#include <vector>
#include "amschain.h"
#include "ParticleData.hh"
#include "BetaNL.hh"
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
     * @brief Get particle mass from PDG ID
     *
     * Supports Geant3 particle IDs and Geant4 PDG IDs if no Geant3 counterpart is found
     *
     * @param pdgId Geant3 particle ID or Geant4 PDG ID
     * @return float Mass in GeV/c^2, 0 if particle not recognized
     */
    float getMass(int pdgId, double charge);

    /**
     * @brief Load particle data from amstreea in ROOT file
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @return std::vector<ParticleData> Vector of particle data, empty if loading fails
     */
    std::vector<ParticleData> loadParticleData(const std::string &inputFile);

    /**
     * @brief Save magnetic field information to ROOT file
     *
     * @param outputFile Path to save the ROOT file containing magnetic field data
     * @return bool True if saving succeeds, false otherwise
     */
    bool saveMagneticField(const std::string &outputFile);

    /**
     * @brief Calculate TOF energy loss for a particle
     *
     * @param particle Particle data
     * @param beta Particle beta
     * @param energyLoss Array to store energy loss in each TOF layer [output]
     * @return bool True if calculation succeeds, false otherwise
     */
    bool calculateTOFEnergyLoss(const ParticleData &particle, double beta, double energyLoss[4]);

    /**
     * @brief Calculate energy loss information
     * 
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing energy loss information
     * @return bool True if calculation succeeds, false otherwise
     */
    bool saveEnergyLoss(const std::string &inputFile, const std::string &outputFile);

    /**
     * @brief Convert ParticleData to BetaNLPars
     * 
     * @param particle Particle data to convert
     * @return BetaNLPars Parameters for beta non-linear reconstruction
     */
    BetaNLPars convertToBetaNLPars(const ParticleData &particle);
}

#endif // __UTIL_HH__