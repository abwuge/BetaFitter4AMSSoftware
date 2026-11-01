#!/usr/bin/env python3
"""Package the research handoff, code, and diagnostic outputs."""
import hashlib
import json
import platform
from pathlib import Path
import sys
import zipfile
import numpy as np
import uproot

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research'
CODE=['flight_calibration_core.py','research_flight_calibration_inventory.py',
      'research_flight_calibration_geometry.py','research_flight_calibration_closure.py',
      'research_flight_calibration_anchor.py','research_flight_calibration_fixed_timing.py',
      'research_flight_calibration_ends.py','verify_flight_calibration_research.py',
      'package_flight_calibration_research.py']


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    verification=json.loads((OUT/'verification.json').read_text())
    if not verification['passed']:raise RuntimeError('Verification has not passed')
    runtime=dict(python=sys.version,numpy=np.__version__,uproot=uproot.__version__,platform=platform.platform())
    (OUT/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    files=[ROOT/'analysis'/name for name in CODE]
    files += [ROOT/'tests/test_flight_calibration_core.py',ROOT/'docs/flight-calibration-handoff.md']
    files += sorted(OUT.glob('*.json'))+sorted(OUT.glob('*.log'))
    files=[p for p in files if p.name!='package_manifest.json']
    manifest={str(p.relative_to(ROOT)):dict(sha256=digest(p),bytes=p.stat().st_size) for p in files}
    if len(manifest)!=len(files):raise RuntimeError('Duplicate package member')
    manifest_text=json.dumps(dict(scope='Research handoff; not a deployable flight beta model',files=manifest),indent=2)+'\n'
    (OUT/'package_manifest.json').write_text(manifest_text)
    path=OUT/'flight_calibration_research.zip'
    temporary=path.with_suffix('.partial.zip')
    with zipfile.ZipFile(str(temporary),'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for p in files:archive.write(str(p),str(p.relative_to(ROOT)))
        archive.writestr('MANIFEST.json',manifest_text)
    with zipfile.ZipFile(str(temporary)) as archive:
        assert archive.testzip() is None
        for name,item in manifest.items():
            assert hashlib.sha256(archive.read(name)).hexdigest()==item['sha256'],name
    temporary.replace(path)
    print(json.dumps(dict(path=str(path),bytes=path.stat().st_size,sha256=digest(path),members=len(manifest)+1,
                          handoff_sha256=digest(ROOT/'docs/flight-calibration-handoff.md')),indent=2),flush=True)

if __name__=='__main__':main()
