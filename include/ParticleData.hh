#ifndef __PARTICLEDATA_HH__
#define __PARTICLEDATA_HH__

struct ParticleData
{
    // Particle properties
    float mass;          // Mass
    float charge;        // Charge
    float momentum;      // Momentum
    float betaLinear;    // Reconstructed beta value using linear approximation
    float innerRigidity; // Inner tracker rigidity
    float betaRigidity;  // Reconstructed beta value using rigidity

    // Particle hit information (TOF)
    float initDir[3];
    float initCoo[3];
    static constexpr int TOF_MAX_HITS = 4;
    float TOF_hitZ[TOF_MAX_HITS];
    float TOF_hitTime[TOF_MAX_HITS];
    float TOF_hitTimeError[TOF_MAX_HITS];
    float TOF_hitEdep[TOF_MAX_HITS];
    float TOF_length[TOF_MAX_HITS];

    // Particle hit information (Tracker)
    static constexpr int TRACKER_MAX_HITS = 9;
    float TRACKER_hitX[TRACKER_MAX_HITS];
    float TRACKER_hitY[TRACKER_MAX_HITS];
    float TRACKER_hitZ[TRACKER_MAX_HITS];
    float TRACKER_hitError[TRACKER_MAX_HITS];

    // MC truth information
    bool isMC;          // Flag to indicate if it's MC particle
    int mcGeantId;      // MC particle Geant 3/4 ID
    float mcMass;       // MC truth mass
    int mcCharge;       // MC truth charge
    float mcMomentum;   // MC truth momentum
    float mcBeta;       // MC truth beta
    float mcInitCoo[3]; // MC truth initial position (cm)
    float mcInitDir[3]; // MC truth initial direction cosines

    // Constructor to initialize arrays and particle properties with appropriate default values
    ParticleData() : mass(0.0f),
                     charge(0),
                     momentum(0.0f),
                     betaLinear(0.0f),
                     innerRigidity(0.0f),
                     betaRigidity(0.0f),
                     initDir{},
                     initCoo{},
                     TOF_hitZ{},
                     TOF_hitTime{},
                     TOF_hitTimeError{},
                     TOF_hitEdep{},
                     TOF_length{},
                     TRACKER_hitX{},
                     TRACKER_hitY{},
                     TRACKER_hitZ{},
                     TRACKER_hitError{},
                     isMC(false),
                     mcGeantId(0),
                     mcMass(0.0f),
                     mcCharge(0),
                     mcMomentum(0.0f),
                     mcBeta(0.0f),
                     mcInitCoo{},
                     mcInitDir{}
    {
    }
};

#endif // __PARTICLEDATA_HH__