#!/usr/bin/env python3
"""Extract JSON from dashboard.html and save two matplotlib PNGs."""
import json
import re
import pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HTML = pathlib.Path('thermoembo_run_stable2/dashboard.html')
text = HTML.read_text()
# Support both old API (const d = {...}) and new API (const entries = [...])
m = re.search(r'const entries\s*=\s*(\[.*?\]);', text, re.DOTALL)
if m:
    entries = json.loads(m.group(1))
    d = entries[0]['series']
else:
    m = re.search(r'const d\s*=\s*(\{.*?\});', text, re.DOTALL)
    d = json.loads(m.group(1))

t = d['times']


def save_fig(fname, ylabel, mean, lo, hi, color):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.fill_between(t, lo, hi, color=color, alpha=0.15, label='min–max')
    ax.plot(t, mean, color=color, lw=2, marker='o', ms=5, label='mean')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(60))
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {fname}')


save_fig('temperature_tumor.png',
         'Temperature (hK)',
         d['temp_mean'], d['temp_min'], d['temp_max'],
         '#DC3232')

save_fig('concentration_tumor.png',
         'Embolic fluid concentration',
         d['conc_mean'], d['conc_min'], d['conc_max'],
         '#3264DC')
