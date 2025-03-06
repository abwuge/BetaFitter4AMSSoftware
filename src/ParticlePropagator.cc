#include "ParticlePropagator.hh"
#include <cmath>
#include <stdexcept>

ParticlePropagator::ParticlePropagator(const AMSPoint &pos,
                                       const AMSDir &dir,
                                       double momentum,
                                       double mass,
                                       int charge)
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
        tof_z[i] = data.TOF_hitZ[i];
}

double ParticlePropagator::GetBeta() const
{
    return sqrt(1.0 - (_mass * _mass) / (_energy * _energy));
}

double ParticlePropagator::GetGamma() const
{
    return _energy / _mass;
}

double ParticlePropagator::GetMomentum() const
{
    return _momentum;
}

void ParticlePropagator::UpdateWithEnergyLoss(const AMSPoint &start_point,
                                              const AMSDir &direction,
                                              double z_target)
{
    // Calculate energy loss using Betalhd
    double energy_loss = _rigidity ? Betalhd::CalculateEnergyLoss(start_point, direction, _rigidity,
                                                                  start_point.z(), z_target,
                                                                  _mass, _chrg, 0)
                                   : 0;

    // Update energy ensuring it stays above rest mass
    _energy = std::max(_mass, _energy - energy_loss);

    // Update momentum and rigidity
    _momentum = sqrt(_energy * _energy - _mass * _mass);
    _rigidity = _momentum / _chrg;
}

double ParticlePropagator::PropagateToZ(double z_target)
{
    // Save current state
    AMSPoint start_point = GetP0();
    AMSDir direction = GetDir();

    // Propagate to target plane
    double len = TrProp::Propagate(z_target);
    if (len < 0)
        return -1;

    // Update kinematics with energy loss
    // UpdateWithEnergyLoss(start_point, direction, z_target);

    return GetBeta();
}

bool ParticlePropagator::PropagateToTOF(double hitX[4], double hitY[4],
                                        double hitTime[4], double pathLength[4])
{
    double total_length = 0;
    const int steps_per_layer = 10; // Number of sub-steps between TOF layers

    for (int i = 0; i < 4; ++i)
    {
        double z_start = _p0z;
        double z_target = tof_z[i];
        double dz = (z_target - z_start) / steps_per_layer;
        double layer_length = 0;
        double layer_time = 0;
        double current_z_target = z_start;

        // Propagate in small steps
        for (int step = 0; step < steps_per_layer; ++step)
        {
            double current_beta = GetBeta();
            if (current_beta <= 0)
                return false;

            AMSPoint start_point = GetP0();
            AMSDir start_dir = GetDir();

            current_z_target += dz;

            double len = TrProp::Propagate(current_z_target);
            if (len < 0)
                return false;

            layer_length += len;
            layer_time += len / (current_beta * SPEED_OF_LIGHT);

            // Update energy loss after each step
            UpdateWithEnergyLoss(start_point, start_dir, current_z_target);
        }

        // Store hit position and direction for this layer
        _hitPoints[i] = GetP0();
        _hitDirs[i] = GetDir();

        // Record hit coordinates
        hitX[i] = _hitPoints[i].x();
        hitY[i] = _hitPoints[i].y();

        // Calculate path length and time
        total_length += layer_length;
        pathLength[i] = total_length;
        hitTime[i] = (i ? hitTime[i - 1] : 0) + layer_time;
    }

    return true;
}

bool ParticlePropagator::resetPropagator(double beta)
{
    if (beta <= 0 || beta >= 1)
    {
        printf("ParticlePropagator: Invalid beta value. Must be between 0 and 1.\n");
        return false;
    }

    _momentum = _mass * beta / sqrt(1 - beta * beta);
    _momentum *= _chrgSign;
    _energy = sqrt(_mass * _mass + _momentum * _momentum);
    _rigidity = _momentum / _chrg;

    _p0x = _initPos.x();
    _p0y = _initPos.y();
    _p0z = _initPos.z();
    _dxdz = (_initDir.z() != 0) ? _initDir.x() / _initDir.z() : 0;
    _dydz = (_initDir.z() != 0) ? _initDir.y() / _initDir.z() : 0;

    return true;
}
