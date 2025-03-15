#ifndef __BETANL_HH__
#define __BETANL_HH__

#include <TMath.h>

#include "TrFit.h"

class BetaNLPars
{
public:
    static constexpr size_t nTOF = 4;                    // Number of TOF hits
    static constexpr double SPEED_OF_LIGHT = 29.9792458; // Speed of light in cm/ns

public:
    // Constructors & Destructors
    // ---------------------------------------------------------------------------

    /**
     * Default constructor
     * @note Default values:
     * - Position: (0, 0, 0)
     * - Direction: (0, 0, 1)
     * - Beta: 0.8
     * - Mass: 0.938 GeV/c^2
     * - Charge: 1
     */
    BetaNLPars() {};

    /**
     * Constructor with parameters
     * @param pos Initial position in cm of the particle
     * @param dir Initial direction of the particle
     * @param beta Initial beta of the particle
     * @param mass Mass in GeV/c^2 of the particle
     * @param charge Charge in e of the particle (must not be 0)
     * @param zTOF Z positions in cm of TOF hits
     * @param energyDeposited Energy deposited in GeV at TOF hits
     * @param hitTime Hit times in ns at TOF hits
     * @param hitTimeError Hit time errors in ns at TOF hits
     */
    BetaNLPars(
        AMSPoint pos,
        AMSDir dir,
        double beta,
        double mass,
        int charge,
        std::vector<double> zTOF,
        std::vector<double> energyDeposited,
        std::vector<double> hitTime,
        std::vector<double> hitTimeError);

    /**
     * Constructor with parameters
     * @param pos Initial position in cm of the particle
     * @param dir Initial direction of the particle
     * @param beta Initial beta of the particle
     * @param mass Mass in GeV/c^2 of the particle
     * @param charge Charge in e of the particle (must not be 0)
     * @param zTOF Z positions in cm of TOF hits
     * @param energyDeposited Energy deposited in GeV at TOF hits
     * @param hitTime Hit times in ns at TOF hits
     * @param hitTimeError Hit time errors in ns at TOF hits
     */
    BetaNLPars(
        const AMSPoint pos,
        const AMSDir dir,
        const double beta,
        const double mass,
        const int charge,
        const double zTOF[nTOF],
        const double energyDeposited[nTOF],
        const double hitTime[nTOF],
        const double hitTimeError[nTOF]);

    /**
     * Constructor with parameters
     * @param pos Initial position in cm of the particle
     * @param dir Initial direction of the particle
     * @param beta Initial beta of the particle
     * @param mass Mass in GeV/c^2 of the particle
     * @param charge Charge in e of the particle (must not be 0)
     * @param zTOF Z positions in cm of TOF hits
     * @param energyDeposited Energy deposited in GeV at TOF hits
     * @param hitTime Hit times in ns at TOF hits
     * @param hitTimeError Hit time errors in ns at TOF hits
     */
    BetaNLPars(
        const AMSPoint pos,
        const AMSDir dir,
        const double beta,
        const double mass,
        const int charge,
        const float zTOF[nTOF],
        const float energyDeposited[nTOF],
        const float hitTime[nTOF],
        const float hitTimeError[nTOF]);

    /**
     * Destructor
     */
    virtual ~BetaNLPars() {};

    // Getters
    // ---------------------------------------------------------------------------

    /**
     * Get the beta (v/c) of the particle
     * @return Beta value
     */
    double Beta() const { return _beta; }

    /**
     * Get the mass in GeV/c^2 of the particle
     * @return Mass (GeV/c^2)
     */
    double Mass() const { return _mass; }

    /**
     * Get the charge in e of the particle
     * @note Charge always > 0
     * @note Refer to BetaNL::Momentum() for sign
     * @return Charge
     */
    int Charge() const { return _charge; }

    /**
     * Get the position in cm of the particle
     * @return Position
     */
    AMSPoint Pos() const { return _pos; }

    /**
     * Get the direction of the particle
     * @return Direction
     */
    AMSDir Dir() const { return _dir; }

    // Setters
    // ---------------------------------------------------------------------------

    /**
     * DON'T ALLOW SETTERS FOR NOW
     */

    // Functions
    // ---------------------------------------------------------------------------

    /**
     * Get the momentum in GeV/c of the particle
     * @note > 0 for positive charge, < 0 for negative charge
     * @return Momentum (GeV/c)
     */
    double Momentum() const { return _mass * _beta / sqrt(1 - _beta * _beta); }

    /**
     * Get the rigidity in GeV/c of the particle
     * @note > 0 for positive charge, < 0 for negative charge
     * @return Rigidity (GeV/c)
     */
    double Rigidity() const { return Momentum() / _charge; }

private:
    BetaNLPars(AMSPoint pos, AMSDir dir, double beta, double mass, int charge);

private:
    // Particle Information
    // ---------------------------------------------------------------------------
    AMSPoint _pos = AMSPoint(0, 0, 0); // Position in cm of the particle
    AMSDir _dir = AMSDir(0, 0, 1);     // Direction of the particle
    double _beta = 0.8;                // Beta of the particle
    double _mass = 0.938;              // Mass in GeV/c^2 of the particle
    int _charge = 1;                   // Charge in e of the particle (must not be 0)

    // Hit Information
    // ---------------------------------------------------------------------------
    std::vector<double> _zTOF;            // Z positions in cm of TOF hits
    std::vector<double> _energyDeposited; // Energy deposited in GeV at TOF hits
    std::vector<double> _hitTime;         // Hit times in ns at TOF hits
    std::vector<double> _hitTimeError;    // Hit time errors in ns at TOF hits

    friend class BetaNL;
};

class BetaNL
{
public:
    // Constructors & Destructors
    // ---------------------------------------------------------------------------

    /**
     * Constructor with BetaNLPars
     * @param pars Parameters for the beta non-linear reconstruction
     */
    BetaNL(BetaNLPars pars, double energyLossScale = 2)
        : _pars(std::make_shared<BetaNLPars>(pars)),
          _energyLossScale(energyLossScale) {};

    /**
     * Destructor
     */
    virtual ~BetaNL() {};

    // Getters
    // ---------------------------------------------------------------------------

    /**
     * Get the reconstructed beta value by the non-linear method
     * @return Reconstructed beta value
     */
    double Beta() { return 1 / InvBeta(); }

    /**
     * Get the reconstructed 1/beta value by the non-linear method
     * @return Reconstructed 1/beta value
     */
    double InvBeta() { return reconstruct(); }

    // Functions
    // ---------------------------------------------------------------------------

    /**
     * Use Monte Carlo beta to calculate energy loss scale factor
     * @param mcBeta Monte Carlo beta value
     * @return Energy loss scale factor
     */
    double EnergyLossScale(double mcBeta);

private:
    TrProp Propagator(const double beta) const;                 // Get the particle propagator with given beta
    std::vector<double> propagate(const double beta) const;     // Propagate the particle with given beta
    double betaChi2(const double *invBeta) const;               // Calculate chi-square for beta reconstruction
    double scaleChi2(const double *scale, const double mcBeta); // Calculate chi-square for energy loss scale
    double reconstruct();                                       // Reconstruct the 1/beta value

private:
    std::shared_ptr<BetaNLPars> _pars;          // Parameters for the beta non-linear reconstruction
    std::shared_ptr<double> _invBeta = nullptr; // Reconstructed 1/beta value
    double _energyLossScale = 2;                // Energy loss scale factor
};

#endif // __BETANL_HH__