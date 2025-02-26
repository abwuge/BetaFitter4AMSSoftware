#include "DataProcessor.hh"
#include <iostream>

DataProcessor::DataProcessor(const std::string &outputFileName)
{
    std::string fileName = outputFileName;
    if (fileName.substr(fileName.find_last_of(".") + 1) != "root")
    {
        fileName += ".root";
    }
    outputFile.reset(new TFile(fileName.c_str(), "RECREATE"));
    setupTree();
}

DataProcessor::~DataProcessor()
{
    if (outputFile)
    {
        outputFile->Write();
        outputFile->Close();
    }
}

void DataProcessor::setupTree()
{
    tree = new TTree("particle", "Particle Information");
    tree->Branch("charge", &particleData.charge, "charge/F");

    // Use vectors for branches instead of arrays
    tree->Branch("tof_x", &particleData.tof_x);
    tree->Branch("tof_y", &particleData.tof_y);
    tree->Branch("tof_z", &particleData.tof_z);
}

bool DataProcessor::processParticle(ParticleR *particle)
{
    if (!particle)
        return false;

    particleData.charge = particle->Charge;
    for (int tof = 0; tof < 4; tof++)
    {
        particleData.tof_x[tof] = particle->TOFCoo[tof][0];
        particleData.tof_y[tof] = particle->TOFCoo[tof][1];
        particleData.tof_z[tof] = particle->TOFCoo[tof][2];
    }
    tree->Fill();
    return true;
}

bool DataProcessor::processEvents(AMSChain &chain, int maxEvents)
{
    int numEvents = maxEvents > 0 ? maxEvents : chain.GetEntries();

    for (int eventNumber = 0; eventNumber < numEvents; eventNumber++)
    {
        AMSEventR *event = (AMSEventR *)chain.GetEvent(eventNumber);
        if (!event)
            continue;

        int nParticles = event->nParticle();
        for (int i = 0; i < nParticles; i++)
            processParticle(event->pParticle(i));
    }
    return true;
}
