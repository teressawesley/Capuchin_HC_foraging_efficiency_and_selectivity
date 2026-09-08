import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

root = Path(__file__).resolve().parents[1]
cols = ['man_hands_duration_s', 'bite_pull_duration_s', 'bite_shell_duration_s', 'hit_surface_duration_s', 'pound_stone_duration_s', 'roll_scrub_duration_s']
counts = {n: 0 for n in range(7)}
with (root / 'generated_data/eff_seq_single_proc_s.csv').open(newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        if row['success'] == '1':
            n = sum(row[c] not in ('', 'NA', 'NaN') and float(row[c]) != 0 for c in cols)
            counts[n] += 1
assert counts[0] == 0
fig, ax = plt.subplots(figsize=(8, 5), dpi=180)
x = list(range(1, 7))
y = [counts[n] for n in x]
bars = ax.bar(x, y, width=.66, color='#377D91')
for bar, value in zip(bars, y):
    ax.text(bar.get_x()+bar.get_width()/2, value-2 if value else .7,
            f'n={value}', ha='center', va='top' if value else 'bottom',
            color='white' if value else '#374151', fontsize=12, fontweight='bold')
ax.set(title=f'Successful sequences (success = 1; N = {sum(y)})',
       xlabel='Number of techniques occurring', ylabel='Number of sequences',
       xticks=x, ylim=(0, 47))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.spines[['top','right']].set_visible(False)
ax.set_axisbelow(True)
ax.yaxis.grid(True, alpha=.18)
fig.text(.5, .025, 'Nonzero, non-NA durations across six processing techniques, including roll/scrub.',
         ha='center', fontsize=9, color='#4B5563')
fig.tight_layout(rect=(0,.06,1,1))
fig.savefig(root / 'outputs/success_technique_counts.png', bbox_inches='tight')
