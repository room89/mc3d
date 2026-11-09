"""Plot estimator calibration and provisional precision gates; requires matplotlib."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser(description=__doc__);p.add_argument('calibration',type=Path);p.add_argument('reanalysis',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
c=json.loads(a.calibration.read_text());r=json.loads(a.reanalysis.read_text());cases=c['noisy'];x=np.arange(4)
labels=[f"{v['shape']}\nAR={v['phi']}" for v in cases]
fig,ax=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
ax[0].bar(x-.18,[100*v['mean_relative_bias'][1] for v in cases],.36,label='Estimated 50% front')
ax[0].bar(x+.18,[100*v['oracle_alignment_relative_bias'][1] for v in cases],.36,label='Known true front (control)')
ax[0].axhline(0,color='black',lw=.5);ax[0].set(xticks=x,xticklabels=labels,ylabel='Mean relative bias (%)',title='Gradient thickness: alignment bias');ax[0].legend(fontsize=8)
ax[1].bar(x-.18,[100*v['statistical_coverage'][1] for v in cases],.36,label='Statistical radius')
ax[1].bar(x+.18,[100*v['with_processing_coverage'][1] for v in cases],.36,label='+ processing sensitivity')
ax[1].axhline(95,color='red',ls='--',lw=1);ax[1].set(xticks=x,xticklabels=labels,ylim=(0,112),ylabel='Truth inside interval (%)',title='Coverage: only 16 experiments per case');ax[1].legend(fontsize=8,loc='lower right')
names=['width_10_90','gradient_thickness','density'];g=r['comparisons']['C vs sparta-C'];ratios=[g[n]['uncertainty']/g[n]['precision_limit'] for n in names]
ax[2].bar(['10–90% width','Gradient','Density'],ratios,color='tab:orange');ax[2].axhline(1,color='red',ls='--',lw=1)
for i,v in enumerate(ratios):ax[2].text(i,v+.15,f'{v:.2f}×',ha='center')
ax[2].set(ylabel='Estimated radius / allowed radius',title='MC3D C vs SPARTA C: precision unmet')
for a0 in ax:a0.grid(axis='y',alpha=.2)
fig.suptitle('Experimental mean-profile analysis — not an acceptance result',fontsize=15)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
