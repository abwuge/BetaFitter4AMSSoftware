"""Measured-input batch API for the linked L5/RICH beta estimator."""
from pathlib import Path
from typing import Dict, Mapping
import numpy as np
import rich_transport as rt


class RichBetaEstimator:
    """Return center beta, RICH beta and seven expected material losses.

    Arrays use the amsd69n branch layout. The center fit is restricted to
    reconstructed carbon and 2.5 <= R < 6 GV. Charge/track-quality selection
    is the caller's responsibility; unsupported rigidity is rejected here.
    MC truth and RICH measurements are neither required nor read.
    """
    REQUIRED = ('tof_tl','tof_leng','tof_pos','tof_edep','tof_barid',
                'tk_ql','tk_rigidityk','tk_pos','tk_dir')

    def __init__(self, model_path: str = None, material_path: str = rt.MAP):
        default = rt.ROOT / 'results/rich_transport/center_model.json'
        self.center = rt.CenterModel(Path(model_path) if model_path else default)
        self.material_path = material_path
        self._material = None

    def predict(self, measured: Mapping[str, np.ndarray], batch_size: int = 16000) -> Dict[str, np.ndarray]:
        if batch_size < 1:
            raise ValueError('batch_size must be positive')
        absent = [k for k in self.REQUIRED if k not in measured]
        if absent:
            raise ValueError('Missing measured branches: ' + ', '.join(absent))
        a = {k:np.asarray(measured[k]) for k in self.REQUIRED}
        n = len(a['tof_tl'])
        if any(len(v) != n for v in a.values()):
            raise ValueError('Measured arrays have different event counts')
        d = rt.center_features(a)
        if not np.isfinite(d['continuous']).all() or not np.isfinite(d['u']).all() or np.any(d['u']<=0):
            raise ValueError('Non-finite or invalid TOF/charge feature input')
        u = self.center.predict(d)
        result = {'beta_center':1/u, 'u_center':u, 'beta_rich_pred':np.full(n,np.nan),
                  'beta_planes':np.full((n,8),np.nan), 'loss_GeV':np.full((n,7),np.nan),
                  'transport_valid':np.zeros(n,bool), 'inferred_mass_GeV':np.full(n,np.nan),
                  'rich_target_z_cm':np.full(n,np.nan), 'NaF_predicted':np.zeros(n,bool)}
        if n and self._material is None:
            self._material = rt.MaterialMap(self.material_path)
        for start in range(0,n,batch_size):
            sl = slice(start,min(start+batch_size,n))
            geo = self._material.build_columns(a['tk_pos'][sl],a['tk_dir'][sl])
            pred = rt.propagate(u[sl],d['R'][sl],geo)
            result['beta_rich_pred'][sl] = pred['beta_planes'][:,-1]
            result['beta_planes'][sl] = pred['beta_planes']
            result['loss_GeV'][sl] = pred['loss_GeV']
            result['transport_valid'][sl] = pred['valid']
            result['inferred_mass_GeV'][sl] = pred['mass_GeV']
            result['rich_target_z_cm'][sl] = geo['target_z']
            result['NaF_predicted'][sl] = geo['NaF_predicted']
        return result
