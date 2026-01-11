#ifndef __UTIL_HH__
#define __UTIL_HH__

#include <string>
#include <vector>
#include "amschain.h"
#include "ParticleData.hh"
#include "BetaNL.h"
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
     * @brief Save beta reconstruction results to ROOT file
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing fit results
     * @return bool True if fit succeeds, false otherwise
     */
    bool saveBeta(const std::string &inputFile, const std::string &outputFile, double energyLossScale = 2);

    /**
     * @brief Save magnetic field information to ROOT file
     *
     * @param outputFile Path to save the ROOT file containing magnetic field data
     * @return bool True if saving succeeds, false otherwise
     */
    bool saveMagneticField(const std::string &outputFile);

    /**
     * @brief Calculate energy loss information
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing energy loss information
     * @return bool True if calculation succeeds, false otherwise
     */
    bool saveEnergyLoss(const std::string &inputFile, const std::string &outputFile);

    /**
     * @brief Calculate energy loss scale factor
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing energy loss scale factor
     * @return bool True if calculation succeeds, false otherwise
     */
    bool saveEnergyLossScale(const std::string &inputFile, const std::string &outputFile);

    /**
     * @brief Benchmark BetaNL::Beta() function average CPU time
     * 
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the benchmark results
     * @param energyLossScale Energy loss scale factor
     * @return bool True if benchmark succeeds, false otherwise
     */
    bool benchmarkBetaNL(const std::string &inputFile, const std::string &outputFile, double energyLossScale = 2);

    bool saveBetaDiff(const std::string &inputFile, const std::string &outputFile, double energyLossScale = 2);

    void test();
}

#endif // __UTIL_HH__