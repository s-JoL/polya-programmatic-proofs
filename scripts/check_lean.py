from pathlib import Path
import subprocess,re,json,hashlib
R=Path(__file__).resolve().parent.parent; B=R/'build';B.mkdir(exist_ok=True)
records=[]
for name in ['FullTorus','K3Interfaces']:
 p=R/(name+'.lean');s=p.read_text()
 assert not re.search(r'\b(sorry|admit|native_decide|unsafe)\b',s)
 assert not re.search(r'^\s*(axiom|opaque)\s',s,re.M)
 out=subprocess.run(['lake','env','lean','-j1','-s8192','-DwarningAsError=true',p.name],cwd=R,capture_output=True,text=True)
 log=out.stdout+out.stderr;(B/(name+'.log')).write_text(log)
 assert out.returncode==0 and 'sorryAx' not in log,log
 allowed={'propext','Classical.choice','Quot.sound'}
 for names in re.findall(r'depends on axioms:\s*\[([^\]]*)\]',log):
  assert {n.strip() for n in names.split(',')} <= allowed, names
 records.append({'module':name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'axioms':[l for l in log.splitlines() if 'depends on axioms' in l or 'does not depend' in l]})
# The actual first minor is 240 at (-1,1); the false lower bound must fail.
first=next(l for l in (R/'FullTorus.lean').read_text().splitlines() if l.startswith('def minor1 '))
bad=B/'NegativeControl.lean';bad.write_text('import Mathlib.Basic.Real.Basic\nimport Mathlib.Tactic.NormNum\n'+first+'\nexample : (241:ℝ) ≤ minor1 (-1) 1 := by norm_num [minor1]\n')
o=subprocess.run(['lake','env','lean','-j1','-s8192','-R',str(B),str(bad)],cwd=R,capture_output=True,text=True)
(B/'NegativeControl.log').write_text(o.stdout+o.stderr)
assert o.returncode!=0 and 'unsolved goals' in o.stdout+o.stderr
(B/'lean-receipt.json').write_text(json.dumps({'status':'PASS','fresh_kernel_checks':True,'records':records,'negative_control_rejected':True,'full_pde_formalization':False},indent=2))
print('PASS: fresh Lean kernel checks; false corner bound rejected')
