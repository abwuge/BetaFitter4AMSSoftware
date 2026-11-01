#!/usr/bin/env python3
"""Package this study, without raw detector data or temporary input caches."""
import json
from pathlib import Path
import zipfile
import rich_transport as rt


def main():
    root=rt.ROOT;dest=root/'results/rich_transport'
    output=root/'results/rich_transport_validation.zip'
    names=['rich_transport.py','rich_beta_estimator.py','evaluate_rich_transport.py',
           'audit_rich_reference.py','prepare_rich_validation.py','evaluate_rich_charge_reference.py',
           'plot_rich_validation.py','verify_rich_validation.py','package_rich_validation.py']
    files=[root/'analysis'/name for name in names]
    files += [root/'tests/test_rich_transport.py',rt.MODEL,rt.BASE/'run_closure.py',rt.BASE/'run_energy_residualized.py']
    files += list(p for p in dest.rglob('*') if p.is_file())
    with zipfile.ZipFile(str(output),'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(set(files)):
            z.write(str(p),'beta_rich_validation/'+str(p.relative_to(root)))
        z.writestr('beta_rich_validation/PACKAGE.txt',
                   'This archive contains the new analysis code, frozen model, numerical results,\n'
                   'event-level predictions and figures. Original ROOT ntuples, the CVMFS material\n'
                   'map and pre-existing training/test caches are not included. Full inference\n'
                   'and source-level reproduction require the existing sdiat project environment.\n'
                   'The saved prediction tables permit independent metric/figure review.\n')
    with zipfile.ZipFile(str(output)) as z:
        if z.testzip() is not None:raise RuntimeError('Archive CRC failure')
        print(json.dumps({'path':str(output),'bytes':output.stat().st_size,'members':len(z.namelist())}),flush=True)

if __name__=='__main__':main()
