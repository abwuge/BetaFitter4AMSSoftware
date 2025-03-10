#include "ParticlePropagator.hh"
#include <cmath>
#include <stdexcept>

ParticlePropagator::ParticlePropagator(const AMSPoint &pos,
                                       const AMSDir &dir,
                                       double momentum,
                                       double mass,
                                       double charge)
    : TrProp(pos, dir), _initPos(pos), _initDir(dir), _energyLossScale(1.0)
{
    if (charge == 0)
    {
        printf("ParticlePropagator: Due to energy loss calculations, "
               "charge cannot be zero. Setting charge to 1.\n");
        charge = 1;
    }

    _chrgSign = (momentum > 0) ? 1 : -1;

    _rigidity = momentum / charge;
    _momentum = momentum;
    _energy = sqrt(mass * mass + momentum * momentum);
    SetMassChrg(mass, charge);
}

ParticlePropagator::ParticlePropagator(const ParticleData &data)
    // : ParticlePropagator(AMSPoint(data.mcCoo[0], data.mcCoo[1], data.mcCoo[2]),
    //                      AMSDir(data.mcDir[0], data.mcDir[1], data.mcDir[2]),
    //                      data.mcMomentum, data.mcMass, data.mcCharge)
    : ParticlePropagator(AMSPoint(data.TRACKER_hitX[0], data.TRACKER_hitY[0], data.TRACKER_hitZ[0]),
                         AMSDir(data.TRACKER_dir[0], data.TRACKER_dir[1], data.TRACKER_dir[2]),
                         data.momentum, data.mass, data.charge)
{
    for (int i = 0; i < 4; ++i)
    {
        tof_z[i] = data.TOF_hitZ[i];
        energyLoss[i] = data.TOF_hitEdep[i];
    }

    for (int i = 0; i < 9; ++i)
        tracker_z[i] = data.TRACKER_hitZ[i];
}

double ParticlePropagator::GetBeta() const
{
    return sqrt(1.0 - (_mass * _mass) / (_energy * _energy));
}

double ParticlePropagator::GetMomentum() const
{
    return _momentum;
}

void ParticlePropagator::UpdateWithEnergyLoss(int i)
{
    // Update energy ensuring it stays above rest mass
    _energy = std::max(_mass, _energy - energyLoss[i] * _energyLossScale);

    // Update momentum and rigidity
    _momentum = sqrt(_energy * _energy - _mass * _mass);
    _rigidity = _momentum / _chrg;
}

bool ParticlePropagator::PropagateToTOF(double trackerHitX[ParticleData::TRACKER_MAX_HITS],
                                        double trackerHitY[ParticleData::TRACKER_MAX_HITS],
                                        double TOFHitTime[ParticleData::TOF_MAX_HITS])
{
    double total_time = 0;
    int i = 0, j = 1; // i for TOF, j for tracker

    // Propagate according to z position (descending order)
    while (i < ParticleData::TOF_MAX_HITS || j < ParticleData::TRACKER_MAX_HITS)
    {
        // Get next target z position
        double tof_next = (i < ParticleData::TOF_MAX_HITS) ? tof_z[i] : -999999;
        double tracker_next = (j < ParticleData::TRACKER_MAX_HITS) ? tracker_z[j] : -999999;

        // Check if we're done
        if (tof_next < -99999 && tracker_next < -99999)
            break;

        // Determine which layer to propagate to
        bool is_tof = tof_next > tracker_next;
        double z_target = is_tof ? tof_next : tracker_next;

        // Propagate to the layer
        double layer_length = TrProp::Propagate(z_target);
        if (layer_length < 0)
            return false;

        // Update for TOF or Tracker
        double layer_time = layer_length / (GetBeta() * SPEED_OF_LIGHT);
        total_time += layer_time;
        
        if (is_tof)
        {
            TOFHitTime[i] = total_time;
            UpdateWithEnergyLoss(i);
            i++;
        }
        else
        {
            trackerHitX[j] = _p0x;
            trackerHitY[j] = _p0y;
            j++;
        }
    }

    return true;
}

bool ParticlePropagator::resetPropagator(double beta)
{
    _momentum = _mass * beta / sqrt(1 - beta * beta);
    _momentum *= _chrgSign;
    _energy = sqrt(_mass * _mass + _momentum * _momentum);
    _rigidity = _momentum / _chrg;

    _p0x = _initPos.x();
    _p0y = _initPos.y();
    _p0z = _initPos.z();
    _dxdz = _initDir.x() / _initDir.z();
    _dydz = _initDir.y() / _initDir.z();

    return true;
}
