"""Replay an archived manifest using verified inputs in a new output directory."""
import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(65536),b''):h.update(block)
    return h.hexdigest()


def replay(manifest, output, binary=None, allow_different=False):
    manifest=Path(manifest).resolve();data=json.loads(manifest.read_text());config=data['effective_config'].copy()
    binary=Path(binary or data['executable']).resolve()
    if sha(binary)!=data['binary_sha256'] and not allow_different:raise ValueError('Binary SHA-256 mismatch')
    for key in ('precompute_snapshot','final_snapshot','final_cell_snapshot'):
        value=config.get(key,'');path=Path(value)
        if value and (path.is_absolute() or '..' in path.parts):config[key]=key+'_'+path.name
    original_runtime=config.pop('runtime_control','')
    config['runtime_control']=''
    units=None
    references=[]
    for item in data['inputs']:
        path=manifest.parent/item['path']
        if sha(path)!=item['sha256']:raise ValueError('Input SHA-256 mismatch: '+str(path))
        if item['role']=='geometry':config['geometry_file']=str(path)
        if item['role']=='reference_data':references.append(str(path))
        if item['role']=='units':units=json.loads(path.read_text())
    with tempfile.TemporaryDirectory(prefix='mc3d-replay-') as tmp:
        tmp=Path(tmp)
        if units:
            units.pop('source_file',None)
            if 'runtime_control' in units['normalized_parameters']:units['normalized_parameters']['runtime_control']=''
            if 'geometry_file' in units['normalized_parameters']:units['normalized_parameters']['geometry_file']=config['geometry_file']
            units['provenance']={'replay_manifest':str(manifest),'original_provenance':units.get('provenance'),'allow_different_binary':allow_different}
            (tmp/'units.json').write_text(json.dumps(units));config['units_file']=str(tmp/'units.json')
        provenance={'replay_manifest':str(manifest),'original_binary_sha256':data['binary_sha256'],'allow_different_binary':allow_different,'files':references,'runtime_control_disabled':bool(original_runtime),'original_runtime_policy':data.get('runtime_policy')}
        (tmp/'provenance.json').write_text(json.dumps(provenance))
        config['provenance_file']=str(tmp/'provenance.json')
        (tmp/'solver.json').write_text(json.dumps(config,allow_nan=False,ensure_ascii=False))
        return subprocess.call([str(binary),'--config',str(tmp/'solver.json'),'--output-dir',str(Path(output).resolve())])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('manifest');p.add_argument('--output-dir',required=True);p.add_argument('--binary');p.add_argument('--allow-different-binary',action='store_true');a=p.parse_args()
    raise SystemExit(replay(a.manifest,a.output_dir,a.binary,a.allow_different_binary))
