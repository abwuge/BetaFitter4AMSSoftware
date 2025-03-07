#include "ParticlePropagator.hh"
#include <cmath>
#include <stdexcept>

ParticlePropagator::ParticlePropagator(const AMSPoint &pos,
                                       const AMSDir &dir,
                                       double momentum,
                                       double mass,
                                       double charge)
    : TrProp(pos, dir), _initPos(pos), _initDir(dir)
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
    _energy = std::max(_mass, _energy - energyLoss[i]);

    // Update momentum and rigidity
    _momentum = sqrt(_energy * _energy - _mass * _mass);
    _rigidity = _momentum / _chrg;
}

bool ParticlePropagator::PropagateToTOF(double hitX[4], double hitY[4],
                                        double hitTime[4], double pathLength[4])
{
    double total_length = 0;
    double total_time = 0;

    // Propagate directly to each TOF layer
    for (int i = 0; i < 4; ++i)
    {
        double z_target = tof_z[i];

        // Calculate beta before propagation
        double current_beta = GetBeta();
        if (current_beta <= 0)
            return false;

        // Propagate to the TOF layer
        double layer_length = TrProp::Propagate(z_target);
        if (layer_length < 0)
            return false;

        // Record hit coordinates
        hitX[i] = _p0x;
        hitY[i] = _p0y;

        // Calculate time to reach this layer
        double layer_time = layer_length / (current_beta * SPEED_OF_LIGHT);

        // Update cumulative path length and time
        total_length += layer_length;
        total_time += layer_time;

        pathLength[i] = total_length;
        hitTime[i] = total_time;

        UpdateWithEnergyLoss(i);
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
