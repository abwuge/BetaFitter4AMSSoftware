#ifndef __DATAPROCESSOR_HH__
#define __DATAPROCESSOR_HH__

#include <string>
#include <memory>
#include <TFile.h>
#include <TTree.h>

#include "amschain.h"
#include "ParticleData.hh"
#include "BetaFitter.hh"
#include "ParticlePropagator.hh"

class DataProcessor {
public:
    DataProcessor(const std::string& outputFileName);
    ~DataProcessor();
    
    bool processEvents(AMSChain& chain, int maxEvents = -1);

private:
    void setupTree();
    bool processParticle(ParticleR* particle);
    bool selectMainParticle(AMSEventR* event, int& selectedIndex);
    
    // Member variables
    std::unique_ptr<TFile> outputFile;
    TTree* tree;
    ParticleData particleData;
    
    // Variables for particle selection
    int ibetah{-1};     // Selected BetaH index
    int itrtrack{-1};   // Selected track index
};

#endif // __DATAPROCESSOR_HH__
