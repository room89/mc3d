"""Plot a saved acceptance audit (optional matplotlib dependency)."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('report',type=Path);p.add_argument('audit',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();r=json.loads(a.report.read_text());d=json.loads(a.audit.read_text())
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
left=r['summaries']['C'];right=r['sparta_summaries']['C'];x=left['grid']
y=[100*(u-v) for u,v in zip(left['density'],right['density'])]
u=[100*(u+v) for u,v in zip(left['density_uncertainty'],right['density_uncertainty'])]
ax=axes[0,0];ax.plot(x,y,label='Mean difference');ax.fill_between(x,[v-e for v,e in zip(y,u)],[v+e for v,e in zip(y,u)],alpha=.2,label='Current summed uncertainty')
for bound in (-1,1):ax.axhline(bound,color='red',ls='--',lw=1)
ax.set(title='C: MC3D minus SPARTA density',xlabel='Aligned x / upstream mean free path',ylabel='% of density jump');ax.legend(fontsize=8)
names=list('ABCDEFG');summaries=[r['summaries'][n] for n in names]
ax=axes[0,1];xx=list(range(8));vals=[s['width_10_90'] for s in summaries]+[right['width_10_90']]
ax.errorbar(xx,[v['mean'] for v in vals],yerr=[v['uncertainty'] for v in vals],fmt='o',capsize=4)
ax.set(xticks=xx,xticklabels=names+['SPARTA C'],title='10–90% width: current means and uncertainties',ylabel='Width / upstream mean free path')
ax=axes[1,0];xx=list(range(7))
ax.errorbar([x-.08 for x in xx],[s['gradient_thickness']['mean'] for s in summaries],yerr=[s['gradient_thickness']['uncertainty'] for s in summaries],fmt='o',capsize=3,label='Mean of window thicknesses (statistical U)')
ax.plot([x+.08 for x in xx],[d['mc3d'][n]['aligned_mean_profile_diagnostic']['all_seeds']['gradient_thickness'] for n in names],'s',label='Thickness of mean profile (exploratory)')
ax.set(xticks=xx,xticklabels=names,title='Gradient estimator depends on averaging order',ylabel='Thickness / upstream mean free path');ax.legend(fontsize=8)
ax=axes[1,1];counts=[sum(s['plateau_failed_windows'] for s in d['mc3d'][n]['seeds']) for n in names]
ax.bar(names,counts);ax.set(title='Plateau deviations > 1% in individual windows',ylabel='Windows out of 192')
for i,v in enumerate(counts):ax.text(i,v+.4,str(v),ha='center')
for ax in axes.flat:ax.grid(alpha=.2)
fig.suptitle('Completed runs; physical acceptance remains incomplete',fontsize=15)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
