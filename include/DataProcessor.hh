#ifndef __DATAPROCESSOR_HH__
#define __DATAPROCESSOR_HH__

#include <string>
#include <memory>
#include <TFile.h>
#include <TTree.h>

#include "amschain.h"
#include "ParticleData.hh"

class DataProcessor {
public:
    DataProcessor(const std::string& outputFileName);
    ~DataProcessor();
    
    bool processEvents(AMSChain& chain, int maxEvents = -1);

private:
    void setupTree();
    bool processParticle(ParticleR* particle);

    std::unique_ptr<TFile> outputFile;
    TTree* tree;
    ParticleData particleData;
};

#endif // __DATAPROCESSOR_HH__
