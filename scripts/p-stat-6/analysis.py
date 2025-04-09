#imports
import pandas as pd
import scipy.stats as stats
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd

sys.path.append('C:/Users/mient/Desktop/biology/masters---AD-signaling-pathway-using-Hsp90-inhibitor/scripts/')

from beta_actin_data import beta_actin_data
from pstat6_data import pstat6_data  

# Preparing data
control_beta_actin = beta_actin_data["control"][0]
stimulation_beta_actin = beta_actin_data["stimulation"][0]
treatment_01_beta_actin = beta_actin_data["0.1MM_treatment"][0]
treatment_1_beta_actin = beta_actin_data["1MM_treatment"][0]

# Individual normalized values
normalized_control = [control_beta_actin / p for p in pstat6_data["control"]]
normalized_stimulation = [stimulation_beta_actin / p for p in pstat6_data["stimulation"]]
normalized_treatment_01 = [treatment_01_beta_actin / p for p in pstat6_data["0.1MM_treatment"]]
normalized_treatment_1 = [treatment_1_beta_actin / p for p in pstat6_data["1MM_treatment"]]

# Combine the data for plotting
data = normalized_control + normalized_stimulation + normalized_treatment_01 + normalized_treatment_1
group = ['Control']*len(normalized_control) + ['Stimulation']*len(normalized_stimulation) + ['0.1MM Treatment']*len(normalized_treatment_01) + ['1MM Treatment']*len(normalized_treatment_1)
df = pd.DataFrame({'value': data, 'group': group})

# Add significance markers
# Statistical test
p_values_shapiro = []
for group_data, name in zip([normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1],
                             ["Control", "Stimulation", "0.1MM Treatment", "1MM Treatment"]):
    stat, p = stats.shapiro(group_data)
    p_values_shapiro.append(p)
    print(f'Shapiro-Wilk Test for {name}: Statistic={stat}, p-value={p}')
    if p > 0.05:
        print(f'{name}: Data is normally distributed')
    else:
        print(f'{name}: Data is NOT normally distributed')

# Homogeneity of variances test (Levene's test)
stat, p_levene = stats.levene(normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1)
print("\nLevene’s Test - Statistic:", stat, "p-value:", p_levene)
if p_levene > 0.05:
    print('Variances are equal (can use ANOVA)')
else:
    print('Variances are NOT equal (use Welch ANOVA or non-parametric test)')

# Choose correct statistical test
if all(p > 0.05 for p in p_values_shapiro) and p_levene > 0.05:
    stat, p_anova = stats.f_oneway(normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1)
    print("\nANOVA - Statistic:", stat, "p-value:", p_anova)
    significant_groups = pairwise_tukeyhsd(df['value'], df['group'], alpha=0.05)
elif all(p > 0.05 for p in p_values_shapiro) and p_levene <= 0.05:
    stat, p_welch = stats.ttest_ind(normalized_stimulation, normalized_treatment_01, equal_var=False)
    print('Welch’s t-test - Statistic:', stat, 'p-value:', p_welch)
    significant_groups = pairwise_tukeyhsd(df['value'], df['group'], alpha=0.05)
else:
    stat, p_kw = stats.kruskal(normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1)
    print('Kruskal-Wallis Test - Statistic:', stat, 'p-value:', p_kw)
    significant_groups = pairwise_tukeyhsd(df['value'], df['group'], alpha=0.05)

# Tukey's Test for pairwise comparisons
print("\nTukey's Test Results:")
print(significant_groups)


#----plots----
group_order = ['Control', 'Stimulation', '0.1MM Treatment', '1MM Treatment']
palette = ['#BB8FCE', '#A569BD', '#8E44AD', '#6C3483']  


sns.set_theme(style="whitegrid")

# plotting
plt.figure(figsize=(10, 6))
ax = sns.boxplot(x='group', y='value', data=df, order=group_order, palette=palette, linewidth=1.3)


means = df.groupby('group')['value'].mean().reindex(group_order)
for i, mean in enumerate(means):
    ax.plot(i, mean, marker='o', color='white', markersize=8, markeredgewidth=2, markeredgecolor='#4A235A')


ax.set_title('Phosphorylation of STAT6', fontsize=18, fontweight='bold', color='#4A235A')  # Update title
ax.set_ylabel('P-STAT6 / Beta Actin', fontsize=14, color='#4A235A')  # Update ylabel
ax.set_xlabel('')
ax.tick_params(axis='x', labelsize=12, colors='#4A235A')
ax.tick_params(axis='y', labelsize=12, colors='#4A235A')
sns.despine()

# adding significance markers
significant_pairs = significant_groups.summary().data[1:]
y_max = df['value'].max() + 0.3
offset = 0.05
for pair in significant_pairs:
    group1, group2, p_val = pair[0], pair[1], pair[2]
    if group1 not in group_order or group2 not in group_order:
        continue
    x1 = group_order.index(group1)
    x2 = group_order.index(group2)
    height = y_max + offset
    ax.plot([x1, x1, x2, x2], [height, height + 0.02, height + 0.02, height], lw=1.2, c='#4A235A')
    ax.text((x1 + x2) / 2, height + 0.03, '*' if p_val < 0.05 else 'ns',
            ha='center', va='bottom', fontsize=12, color='#4A235A')
    offset += 0.15

plt.tight_layout()
plt.show()

