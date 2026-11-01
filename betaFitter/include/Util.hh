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
     * @brief Get the isotope-abundance-weighted mass for a reconstructed charge
     *
     * @param charge Reconstructed nuclear charge Z
     * @return float Average mass in GeV/c^2, 0 for unsupported charges
     */
    float getAverageMass(int charge);

    /**
     * @brief Load particle data from amstreea in ROOT file
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param referencePoint Reference point used to derive MC beta
     * @return std::vector<ParticleData> Vector of particle data, empty if loading fails
     */
    std::vector<ParticleData> loadParticleData(
        const std::string &inputFile,
        BetaReferencePoint referencePoint = BetaReferencePoint::AMSCenter);

    /**
     * @brief Save beta reconstruction results to ROOT file
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing fit results
     * @param energyLossScale Scale factor applied to selected TOF energy losses
     * @param trackerEnergyLossScale Scale factor applied to tracker energy losses
     * @param scaleConfig Optional charge-dependent scale configuration file
     * @return bool True if fit succeeds, false otherwise
     */
    bool saveBeta(const std::string &inputFile,
                  const std::string &outputFile,
                  double energyLossScale = 2,
                  double trackerEnergyLossScale = 1,
                  EnergyLossScaleMode energyLossScaleMode = EnergyLossScaleMode::All,
                  BetaReferencePoint referencePoint = BetaReferencePoint::AMSCenter,
                  const std::string &scaleConfig = "none",
                  const std::string &hybridConfig = "none",
                  double materialStep = 0.0);

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
     * @brief Fit per-event TOF and tracker energy loss scale factors
     *
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the ROOT file containing energy loss scale factor
     * @return bool True if the output is written, false otherwise
     */
    bool saveEnergyLossScale(const std::string &inputFile,
                             const std::string &outputFile,
                             EnergyLossScaleMode energyLossScaleMode = EnergyLossScaleMode::All,
                             BetaReferencePoint referencePoint = BetaReferencePoint::AMSCenter);

    /**
     * @brief Export MC samples used to train the hybrid energy-loss model
     *
     * The output contains, for four physically defined detector regions, the
     * charge-aware Bethe/material prediction, MC true energy loss, and matched
     * TOF/tracker deposits. Region boundaries follow the hit TOF-bar surfaces.
     * All energy quantities are stored in GeV.
     *
     * @param inputFile Path to an MC amstreea ROOT file
     * @param outputFile Path to the hybrid training ROOT file
     * @return bool True if at least one valid training sample is written
     */
    bool saveHybridEnergyLossTraining(const std::string &inputFile,
                                      const std::string &outputFile,
                                      double materialStep = 0.0);

    /**
     * @brief Fit hybrid energy-loss calibration parameters from a training tree
     *
     * @param inputFile ROOT file containing a hybridTraining tree
     * @param outputFile Text calibration file consumed by saveBeta
     * @return bool True if all available charge/segment bins are fitted
     */
    bool trainHybridEnergyLoss(const std::string &inputFile,
                               const std::string &outputFile);

    /**
     * @brief Benchmark BetaNL::Beta() function average CPU time
     * 
     * @param inputFile Path to the ROOT file containing particle data
     * @param outputFile Path to save the benchmark results
     * @param energyLossScale Energy loss scale factor
     * @return bool True if benchmark succeeds, false otherwise
     */
    bool benchmarkBetaNL(const std::string &inputFile,
                         const std::string &outputFile,
                         double energyLossScale = 2,
                         EnergyLossScaleMode energyLossScaleMode = EnergyLossScaleMode::All,
                         BetaReferencePoint referencePoint = BetaReferencePoint::AMSCenter);

    bool saveBetaDiff(const std::string &inputFile,
                      const std::string &outputFile,
                      double energyLossScale = 2,
                      EnergyLossScaleMode energyLossScaleMode = EnergyLossScaleMode::All,
                      BetaReferencePoint referencePoint = BetaReferencePoint::AMSCenter);

    void test();
}

#endif // __UTIL_HH__
