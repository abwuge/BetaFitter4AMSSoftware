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

    // Original particle properties
    tree->Branch("mass", &particleData.mass, "mass/F");
    tree->Branch("charge", &particleData.charge, "charge/F");
    tree->Branch("momentum", &particleData.momentum, "momentum/F");
    tree->Branch("beta", &particleData.beta, "beta/F");
    tree->Branch("Theta", &particleData.Theta, "Theta/F");
    tree->Branch("Phi", &particleData.Phi, "Phi/F");

    // Hit information - Using fixed size arrays with MAX_HITS
    tree->Branch("hitX", particleData.hitX, Form("hitX[%d]/F", ParticleData::MAX_HITS));
    tree->Branch("hitY", particleData.hitY, Form("hitY[%d]/F", ParticleData::MAX_HITS));
    tree->Branch("hitZ", particleData.hitZ, Form("hitZ[%d]/F", ParticleData::MAX_HITS));
    tree->Branch("hitTime", particleData.hitTime, Form("hitTime[%d]/F", ParticleData::MAX_HITS));
    tree->Branch("hitTimeError", particleData.hitTimeError, Form("hitTimeError[%d]/F", ParticleData::MAX_HITS));

    // MC truth information
    tree->Branch("mcBeta", &particleData.mcBeta, "mcBeta/F");
    tree->Branch("mcMomentum", &particleData.mcMomentum, "mcMomentum/F");
    tree->Branch("mcMass", &particleData.mcMass, "mcMass/F");
    tree->Branch("mcCharge", &particleData.mcCharge, "mcCharge/I");  // 改为I表示int类型
    tree->Branch("mcPdgId", &particleData.mcPdgId, "mcPdgId/I");
    tree->Branch("isMC", &particleData.isMC, "isMC/O");
}

bool DataProcessor::processParticle(ParticleR *particle)
{
    if (!particle)
        return false;

    // Process reconstructed information
    particleData.mass = particle->Mass;
    particleData.charge = particle->Charge;
    particleData.momentum = particle->Momentum;
    particleData.beta = particle->Beta;
    particleData.Theta = particle->Theta;
    particleData.Phi = particle->Phi;

    for (int tof = 0; tof < ParticleData::MAX_HITS; tof++)
    {
        particleData.hitX[tof] = particle->TOFCoo[tof][0];
        particleData.hitY[tof] = particle->TOFCoo[tof][1];
        particleData.hitZ[tof] = particle->TOFCoo[tof][2];
    }

    BetaHR *beta = particle->pBetaH();
    if (!beta)
        return false;

    for (int i = 0; i < ParticleData::MAX_HITS; i++)
    {
        particleData.hitTime[i] = beta->GetTime(i);
        particleData.hitTimeError[i] = beta->GetETime(i);
    }

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

        // Get MC information if available
        if (event->nMCEventg() > 0)
        {
            MCEventgR *mcEvent = event->GetPrimaryMC();
            particleData.isMC = true;
            particleData.mcPdgId = mcEvent->Particle;
            particleData.mcCharge = mcEvent->Charge;
            particleData.mcMass = mcEvent->Mass;
            particleData.mcMomentum = mcEvent->Momentum;
            // Calculate beta from MC momentum and mass
            float energy = sqrt(particleData.mcMomentum * particleData.mcMomentum +
                                particleData.mcMass * particleData.mcMass);
            particleData.mcBeta = particleData.mcMomentum / energy;
        }
        else
        {
            particleData.isMC = false;
            particleData.mcPdgId = 0;
            particleData.mcCharge = 0.0f;
            particleData.mcMass = 0.0f;
            particleData.mcMomentum = 0.0f;
            particleData.mcBeta = 0.0f;
        }
#ifdef false
        int mainParticleIdx = -1;
        if (!selectMainParticle(event, mainParticleIdx))
            continue;

        if (!processParticle(event->pParticle(mainParticleIdx)))
            continue;
#endif
        if (!processParticle(event->GetPrimaryParticle()))
            continue;

        tree->Fill();
    }
    return true;
}

bool DataProcessor::selectMainParticle(AMSEventR *event, int &selectedIndex)
{
    if (!event)
        return false;

#ifdef _GOODTKTOFPAR_
    if (event->nTrTrack() < 1)
        return false;
#endif

    // Arrays to store maximum charge values
    float tkiqm[3] = {0}; // Inner tracker charge
    selectedIndex = -1;
    ibetah = -1;
    itrtrack = -1;

    // Loop over all tracks
    for (int iTrk = 0; iTrk < event->nTrTrack(); iTrk++)
    {
        TrTrackR *trk = event->pTrTrack(iTrk);
        if (!trk)
            continue;

        // Check track quality and pattern
        bool hasL1 = (trk->GetBitPatternJ() & 1);
        bool hasL9 = (trk->GetBitPatternJ() & (1 << 8));
        float ntkq = 0;

#if defined(_HEINNERPRESCALE_) || defined(_IONINNERPRESCALE_)
        ntkq = trk->GetInnerQ();
        if (ntkq > tkiqm[0])
        {
            tkiqm[0] = ntkq; // Inner
        }
#endif

#ifdef _SAVENEGSCALE_
        int nhiti = 0;
        for (int ilay = 1; ilay < 8; ilay++)
        {
            if ((trk->GetBitPatternJ() & (1 << ilay)) > 0)
                nhiti++;
        }
        if (nhiti >= 5 && ntkq > 0.6)
        {
            // Loop over BetaH measurements
            for (int iBeta = 0; iBeta < event->NBetaH(); iBeta++)
            {
                BetaHR *beta = event->pBetaH(iBeta);
                if (!beta)
                    continue;
                if (beta->iTrTrack() != iTrk)
                    continue;

                float obeta = beta->GetBeta();
                float tkrig = trk->GetRigidity();

                if (fabs(obeta) > 0.4 && fabs(tkrig) > 0.7 && tkrig / obeta < 0)
                {
                    // Store best negative particle candidate
                    if (ibetah < 0 || itrtrack < 0)
                    {
                        ibetah = iBeta;
                        itrtrack = iTrk;
                        selectedIndex = iBeta;
                    }
                }
            }
        }
#endif

        // Process BetaH measurements for track
        for (int iBeta = 0; iBeta < event->NBetaH(); iBeta++)
        {
            BetaHR *beta = event->pBetaH(iBeta);
            if (!beta)
                continue;
            if (beta->iTrTrack() != iTrk)
                continue;

            // Beta quality checks
            if (beta->GetBeta() <= 0)
                continue;

            // Select best beta measurement
            if (ibetah < 0 || itrtrack < 0)
            {
                ibetah = iBeta;
                itrtrack = iTrk;
                selectedIndex = iBeta;
            }
            else
            {
                BetaHR *curBeta = event->pBetaH(ibetah);

                int betaTofClusters = 0;
                int curBetaTofClusters = 0;

                for (int i = 0; i < 4; i++)
                {
                    if (beta->GetClusterHL(i))
                        betaTofClusters++;
                    if (curBeta->GetClusterHL(i))
                        curBetaTofClusters++;
                }

                float betaChi2 = beta->GetNormChi2T();
                float curBetaChi2 = curBeta->GetNormChi2T();

                if (betaTofClusters > curBetaTofClusters ||
                    (betaTofClusters == curBetaTofClusters && betaChi2 < curBetaChi2))
                {
                    ibetah = iBeta;
                    itrtrack = iTrk;
                    selectedIndex = iBeta;
                }
            }
        }
    }

#ifdef _GOODTKTOFPAR_
    if (ibetah < 0 || itrtrack < 0)
        return false;
#endif

    for (unsigned int i = 0; i < event->NParticle(); i++)
    {
        ParticleR *p = event->pParticle(i);
        if (p && p->iBetaH() == ibetah)
        {
            selectedIndex = i;
            break;
        }
    }

    if (selectedIndex < 0)
        return false;

    return true;
}
