import sys,uproot,numpy as np
with uproot.open(sys.argv[1]) as f:
 t=f['amstreea'];a=t.arrays(['tk_adc','tk_qlqs','tk_stripx'],library='np')
 assert t.num_entries>0
 for k in ['tk_adc','tk_qlqs']:assert np.count_nonzero(a[k])>0
 print('VALIDATED',t.num_entries,'events',len(t.keys()),'branches')
