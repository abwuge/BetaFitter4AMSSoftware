"""Frozen weak-beta reconstruction with a one-way material transport to RICH.

The learned component is the existing anchored-energy beta model. Transport
uses the native elemental material map and the Bethe mean stopping power.
There is no fit to MC momenta or to the RICH velocity. Losses are expected
interval losses, including passive material, rather than observed deposits.
"""
import json
import sys
from pathlib import Path
import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/c_isotope_proportion_closure'
sys.path.insert(0, str(BASE))
import run_closure as rc
import run_energy_residualized as er

MODEL = ROOT / 'results/c_joint_beta_rigidity_closure/expanded_validation/energy_comparison/energy_frozen_model.json'
MAP = '/cvmfs/ams.cern.ch/Offline/AMSDataDir/v5.00/LAPP/dEdxPDF/g4mscan.root'
Z_PLANES = np.array([1.698, -2.318, -25.212, -29.228, -62.1, -65.2, -69.975])
PLANE_NAMES = ['L5', 'L6', 'L7', 'L8', 'TOF_S3_mean_plane', 'TOF_S4_mean_plane', 'post_S4', 'RICH_emission']
LOSS_NAMES = [a + '_to_' + b for a, b in zip(PLANE_NAMES[:-1], PLANE_NAMES[1:])]
A = np.array([1.0079, 12.0107, 14.0067, 15.9994, 18.9984, 22.9897, 26.9815, 28.0855, 207.2])
ZA = np.array([.99212, .49954, .49976, .50002, .47372, .47847, .48181, .49848, .39575])
I_MEV = np.array([19.2e-6, 78e-6, 82e-6, 95e-6, 115e-6, 149e-6, 166e-6, 173e-6, 823e-6])


def arrays(value):
    if isinstance(value, dict):
        return {k: arrays(v) for k, v in value.items()}
    if isinstance(value, list):
        return np.asarray(value)
    return value


def center_features(a):
    """The measured feature schema of the frozen carbon model."""
    t = np.asarray(a['tof_tl'], float)
    x = np.asarray(a['tof_leng'], float) / rc.C_LIGHT
    cp = x - x.mean(axis=1, keepdims=True)
    ct = t - t.mean(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.sum(cp * ct, axis=1) / np.sum(cp * cp, axis=1)
        res = ct - u[:, None] * cp
        span = np.ptp(x, axis=1)
        p = np.asarray(a['tof_pos'], float)
        dz = p[:, 3, 2] - p[:, 0, 2]
        continuous = np.column_stack([
            p[:, :, 0] / 120., p[:, :, 1] / 120.,
            (p[:, 3, 0] - p[:, 0, 0]) / dz,
            (p[:, 3, 1] - p[:, 0, 1]) / dz, cp, res / span[:, None],
            np.log(np.asarray(a['tof_edep'], float) / 36.),
            2 * np.log(np.asarray(a['tk_ql'][:, 1:8], float) / 6.),
            np.log(span / 4.)])
    rigidity = np.asarray(a['tk_rigidityk'][:, 1, 1], float)
    return {'u': u, 'R': rigidity,
            'bin': np.searchsorted(rc.R_EDGES, rigidity, side='right') - 1,
            'continuous': continuous, 'bars': a['tof_barid']}


class CenterModel:
    def __init__(self, path=MODEL):
        self.state = arrays(json.loads(Path(path).read_text()))

    def predict(self, d):
        if np.any((d['bin'] < 0) | (d['bin'] >= 4)):
            raise ValueError('The frozen beta model supports 2.5 <= R < 6 GV.')
        s = self.state
        b = rc.apply_feature_state(d, s['base_state'])
        z = er.apply_standardizer(er.energy_raw(d, 'fitted', s['response']), s['standardizer'])
        der = np.clip((np.column_stack([z, z*z]) - s['derived_center']) / s['derived_scale'], -5, 5)
        x = np.column_stack([b, der])[:, 1:]
        anchor = s['anchor']
        p = (x - anchor['bin_centers'][d['bin']]).dot(anchor['null_projection']) / anchor['scale']
        correction = p.dot(s['parameters'][1:p.shape[1]+1])
        return d['u'] - correction


def core_selection(a):
    d = center_features(a)
    q = np.where(a['tk_qin'][:, 0, 2] < 2.5, a['tk_q'][:, 1], a['tk_qin'][:, 0, 2])
    ok = np.isfinite(q) & (q >= 5.5) & (q < 6.5) & (a['tof_bpatt'] == 4444)
    ok &= np.isfinite(d['u']) & (d['u'] > 0) & np.all(np.isfinite(d['continuous']), axis=1)
    ok &= (d['bin'] >= 0) & (d['bin'] < 4)
    ok &= np.all((a['tof_barid'] >= 0) & (a['tof_barid'] < 10), axis=1)
    ok &= np.all(np.isfinite(a['tk_pos'][:, 1:8]), axis=(1, 2))
    ok &= (a['tk_pos'][:, 4, 2] > 0) & (a['tk_pos'][:, 5, 2] < 0)
    rel = np.abs(a['tk_erigidity1'][:, 2, 1] * d['R'])
    chi = a['tk_chis1'][:, 2, 2, 1, :]
    ok &= np.isfinite(rel) & (rel > 0) & (rel < .20)
    ok &= np.all(np.isfinite(chi) & (chi > 0) & (chi <= 10), axis=1)
    return ok


def rich_quality(a, strict=False):
    """Reference-only cuts; no TOF/RICH beta-agreement or mass-window cut."""
    b = np.asarray(a['rich_beta'][:, 0], float)
    naf = a['rich_NaF'].astype(bool)
    xy = np.max(np.abs(a['rich_pos'][:, :2]), axis=1)
    rad = np.sum(a['rich_pos'][:, :2] ** 2, axis=1)
    npe = a['rich_npe']
    bc, pc = a['rich_BetaConsistency'], a['rich_PMTChargeConsistency']
    ok = (a['irich'] >= 0) & (a['itrtrack'] >= 0) & (a['rich_itrtrack'] == a['itrtrack'])
    ok &= a['rich_good'].astype(bool) & a['rich_clean'].astype(bool)
    ok &= (a['rich_pb'] >= (.05 if strict else .01))
    ok &= (a['rich_pmt'] >= (4 if strict else 3)) & (a['rich_used'] >= (6 if strict else 4))
    ok &= np.isfinite(b) & (b >= .78) & (b <= 1.02) & np.where(naf, b >= .78, b >= .955)
    ok &= np.isfinite(a['rich_ebeta1']) & (a['rich_ebeta1'] > 0) & (a['rich_ebeta1'] < (.003 if strict else .005))
    ok &= np.isfinite(bc) & (bc >= 0) & (bc <= np.where(naf, .01, .005))
    ok &= np.isfinite(pc) & (pc >= 0) & (pc <= 15)
    ok &= np.all(np.isfinite(npe), axis=1) & (npe[:, 2] > 0)
    ok &= (npe[:, 0] >= .4*npe[:, 2]) & (npe[:, 1] >= np.where(naf, 1., 2.))
    ok &= np.where(naf, xy <= 17, (xy >= 19) & (rad <= 3500) & ~np.isin(a['rich_tile'], [3,7,12,20,87,100,108]))
    return ok & np.all(np.isfinite(a['rich_pos']), axis=1)


def track_at_z(pos, direction, z):
    """Cubic Hermite interpolation of fitted Tracker layer states L5--L9."""
    p = np.asarray(pos[:, 4:9], float)
    v = np.asarray(direction[:, 4:9], float)
    z = np.broadcast_to(z, (len(p),))
    j = np.clip(np.sum(z[:, None] < p[:, 1:, 2], axis=1), 0, 3)
    row = np.arange(len(p))
    pa, pb = p[row, j], p[row, j+1]
    va, vb = v[row, j], v[row, j+1]
    h = pb[:, 2] - pa[:, 2]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (z-pa[:, 2])/h
        da, db = va[:, :2]/va[:, 2, None], vb[:, :2]/vb[:, 2, None]
        xy = ((2*t**3-3*t*t+1)[:,None]*pa[:,:2] + (t**3-2*t*t+t)[:,None]*h[:,None]*da
              + (-2*t**3+3*t*t)[:,None]*pb[:,:2] + (t**3-t*t)[:,None]*h[:,None]*db)
        derivative = ((6*t*t-6*t)[:,None]/h[:,None]*pa[:,:2] + (3*t*t-4*t+1)[:,None]*da
                      + (-6*t*t+6*t)[:,None]/h[:,None]*pb[:,:2] + (3*t*t-2*t)[:,None]*db)
    return xy, np.sqrt(1 + np.sum(derivative**2, axis=1))


def stopping_power(beta, mass, charge=6):
    """Bethe mean in GeV cm^2/g, with exact maximum electron energy transfer.

Density, shell, Barkas and restricted-loss corrections are not fitted here.
The declared carbon interval has mostly beta below 0.96.
"""
    beta = np.clip(beta, .05, 1-1e-10)
    g2 = 1/(1-beta*beta)
    me = .51099895
    m = np.maximum(1000*mass, me)
    tmax = 2*me*beta*beta*g2/(1+2*np.sqrt(g2)*me/m+(me/m)**2)
    bracket = .5*np.log(2*me*beta[:,None]**2*g2[:,None]*tmax[:,None]/I_MEV[None,:]**2) - beta[:,None]**2
    return .000307075*charge**2*ZA[None,:]/beta[:,None]**2*np.maximum(bracket, 0)


class MaterialMap:
    def __init__(self, path=MAP):
        self.path = path
        pieces = []
        with uproot.open(path, object_cache=None, array_cache=None) as f:
            for j in range(1, 10):
                hist = f['hist%d' % j]
                if j == 1:
                    self.edges = [ax.edges() for ax in hist.axes]
                    iz0 = np.searchsorted(self.edges[2], -78., side='right')-1
                    iz1 = np.searchsorted(self.edges[2], 4., side='left')
                    self.z_edges = self.edges[2][iz0:iz1+1]
                pieces.append(hist.values()[:, :, iz0:iz1].astype(np.float32, copy=True))
                del hist
        self.values = np.stack(pieces, axis=-1)
        del pieces

    def column_density(self, pos, direction, zh, zl, samples=5):
        n = len(pos)
        out = np.zeros((n, 9), float)
        inside = np.ones(n, bool)
        length_z = np.maximum(zh-zl, 0)
        for k in range(samples):
            z = zh - length_z*(k+.5)/samples
            xy, metric = track_at_z(pos, direction, z)
            ix = np.searchsorted(self.edges[0], xy[:,0], side='right')-1
            iy = np.searchsorted(self.edges[1], xy[:,1], side='right')-1
            iz = np.searchsorted(self.z_edges, z, side='right')-1
            valid = np.isfinite(xy).all(axis=1) & np.isfinite(metric)
            valid &= (ix>=0)&(ix<self.values.shape[0])&(iy>=0)&(iy<self.values.shape[1])&(iz>=0)&(iz<self.values.shape[2])
            active = length_z > 0
            inside &= ~active | valid
            rows = np.flatnonzero(active & valid)
            out[rows] += self.values[ix[rows],iy[rows],iz[rows]].astype(float)*(length_z[rows]*metric[rows]/samples)[:,None]
        return out*A[None,:], inside

    def build_columns(self, pos, direction):
        n = len(pos)
        xy, _ = track_at_z(pos, direction, np.full(n, -75.1))
        naf = np.max(np.abs(xy), axis=1) < 17.5
        target = np.where(naf, -75.552, -74.705)
        faces = np.r_[Z_PLANES, np.arange(-76., 2.), -74.705, -75.552]
        faces = np.sort(np.unique(faces))[::-1]
        faces = faces[(faces <= Z_PLANES[0]) & (faces >= -75.552)]
        column = np.zeros((n, len(faces)-1, 9), np.float32)
        valid = np.isfinite(pos[:,4:9]).all(axis=(1,2)) & (np.abs(pos[:,4:9,:2])<1000).all(axis=(1,2))
        valid &= np.isfinite(direction[:,4:9]).all(axis=(1,2)) & (np.abs(direction[:,4:9,2])>.05).all(axis=1)
        valid &= (np.diff(pos[:,4:9,2],axis=1)<0).all(axis=1)
        seg = np.empty(len(faces)-1, np.int8)
        for k, (high, low) in enumerate(zip(faces[:-1], faces[1:])):
            zh = np.maximum(high, target)
            zl = np.maximum(low, target)
            column[:,k], inside = self.column_density(pos, direction, zh, zl)
            valid &= inside
            seg[k] = min(int(np.sum(.5*(high+low) < Z_PLANES[1:])), 6)
        return {'column':column, 'segment':seg, 'valid':valid, 'target_z':target, 'NaF_predicted':naf, 'faces':faces}


def propagate(u, rigidity, geometry, material_scale=1., mass_override=None):
    """Return beta at each plane and positive expected losses in GeV.

With no mass override, infer the continuous mass from measured p=Z R and
center beta. No discrete isotope assignment or true particle mass is used.
Unphysical input beta remains flagged rather than being clipped into a peak.
"""
    u = np.asarray(u, float)
    p = 6*np.asarray(rigidity, float)
    valid = geometry['valid'] & np.isfinite(u) & (u > 1) & np.isfinite(p) & (p > 0)
    up = np.where(valid, u, 1.1)
    if mass_override is None:
        mass = p*np.sqrt(up*up-1)
        energy = p*up
    else:
        mass = np.full(len(u), float(mass_override))
        energy = mass/np.sqrt(1-1/(up*up))
    initial_energy = energy.copy()
    losses = np.zeros((len(u), 7))
    beta_planes = np.full((len(u), 8), np.nan)
    beta_planes[:,0] = 1/u
    for k, segment in enumerate(geometry['segment']):
        beta = np.sqrt(np.maximum(1-(mass/energy)**2, 1e-8))
        stopping = stopping_power(beta, mass)
        loss = material_scale*np.sum(stopping*geometry['column'][:,k], axis=1)
        # Midpoint integration uses the updated stopping power, not a fitted scale.
        mid_energy = np.maximum(energy-.5*loss, mass*(1+1e-8))
        mid_beta = np.sqrt(np.maximum(1-(mass/mid_energy)**2, 1e-8))
        loss = material_scale*np.sum(stopping_power(mid_beta,mass)*geometry['column'][:,k], axis=1)
        energy -= loss
        valid &= np.isfinite(energy) & (energy > mass)
        losses[:,segment] += loss
        beta_planes[:,segment+1] = np.sqrt(np.maximum(1-(mass/np.maximum(energy,mass))**2,0))
    beta_planes[~valid] = np.nan
    losses[~valid] = np.nan
    return {'beta_planes':beta_planes, 'loss_GeV':losses, 'valid':valid,
            'mass_GeV':np.where(valid,mass,np.nan), 'initial_energy_GeV':np.where(valid,initial_energy,np.nan),
            'final_energy_GeV':np.where(valid,energy,np.nan)}
