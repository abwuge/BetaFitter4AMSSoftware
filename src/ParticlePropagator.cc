#include "ParticlePropagator.hh"
#include <cmath>
#include <stdexcept>

ParticlePropagator::ParticlePropagator(const AMSPoint &pos,
                                       const AMSDir &dir,
                                       double momentum,
                                       double mass,
                                       int charge)
{
    if (!charge)
    {
        printf("ParticlePropagator: Due to energy loss calculations, charge cannot be zero. Setting charge to 1.\n");
        charge = 1;
    }

    TrProp::TrProp(pos, dir, momentum / charge);
    SetMassChrg(mass, charge);
    _momentum = momentum;
    _energy = sqrt(mass * mass + momentum * momentum);
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
    double energy_loss = Betalhd::CalculateEnergyLoss(start_point, direction, _rigidity,
                                                      start_point.z(), z_target,
                                                      _mass, _chrg, 0);

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
    UpdateWithEnergyLoss(start_point, direction, z_target);

    return GetBeta();
}

bool ParticlePropagator::PropagateToTOF(double hitX[4], double hitY[4],
                                        double hitTime[4], double pathLength[4])
{
    double total_length = 0;

    for (int i = 0; i < 4; ++i)
    {
        // Save current state before propagation
        AMSPoint start_point = GetP0();
        AMSDir start_dir = GetDir();

        // Propagate to TOF layer
        double len = TrProp::Propagate(TOF_Z[i]);
        if (len < 0)
            return false;

        // Store hit position and direction
        _hitPoints[i] = GetP0();
        _hitDirs[i] = GetDir();

        // Record hit coordinates
        hitX[i] = _hitPoints[i].x();
        hitY[i] = _hitPoints[i].y();

        // Calculate path length and time
        total_length += len;
        pathLength[i] = total_length;
        hitTime[i] = total_length / (GetBeta() * SPEED_OF_LIGHT);

        // Update kinematics for next layer (except for last layer)
        if (i < 3)
            UpdateWithEnergyLoss(start_point, start_dir, TOF_Z[i]);
    }

    return true;
}

bool ParticlePropagator::resetPropagator(const AMSPoint &pos,
                                         const AMSDir &dir,
                                         double momentum)
{
    if (momentum <= 0)
        return false;

    _p0x = pos.x();
    _p0y = pos.y();
    _p0z = pos.z();
    _dxdz = (dir.z() != 0) ? dir.x() / dir.z() : 0;
    _dydz = (dir.z() != 0) ? dir.y() / dir.z() : 0;

    _momentum = momentum;
    _energy = sqrt(_mass * _mass + momentum * momentum);
    _rigidity = momentum / _chrg;

    return true;
}