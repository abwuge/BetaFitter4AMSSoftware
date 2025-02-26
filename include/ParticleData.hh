#ifndef __PARTICLEDATA_HH__
#define __PARTICLEDATA_HH__

#include <vector>

struct ParticleData
{
    std::vector<float> tof_x;
    std::vector<float> tof_y;
    std::vector<float> tof_z;

    float charge;

    // Constructor to initialize vectors with appropriate size
    ParticleData() : 
        tof_x(4, 0.0f), 
        tof_y(4, 0.0f), 
        tof_z(4, 0.0f), 
        charge(0.0f) {}
};

#endif // __PARTICLEDATA_HH__