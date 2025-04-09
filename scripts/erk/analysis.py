#imports
import pandas as pd
import scipy.stats as stats
import sys
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import seaborn as sns

sys.path.append('C:/Users/mient/Desktop/biology/masters---AD-signaling-pathway-using-Hsp90-inhibitor/scripts/')


from beta_actin_data import beta_actin_data
from erk_data import erk_data

#preparing data
control_beta_actin = beta_actin_data["control"][0]
stimulation_beta_actin = beta_actin_data["stimulation"][0]
treatment_01_beta_actin = beta_actin_data["0.1MM_treatment"][0]
treatment_1_beta_actin = beta_actin_data["1MM_treatment"][0]

# Individual normalized values
normalized_control = [control_beta_actin / p for p in erk_data["control"]]
normalized_stimulation = [stimulation_beta_actin / p for p in erk_data["stimulation"]]
normalized_treatment_01 = [treatment_01_beta_actin / p for p in erk_data["0.1MM_treatment"]]
normalized_treatment_1 = [treatment_1_beta_actin / p for p in erk_data["1MM_treatment"]]

# Normality test (Shapiro-Wilk) for each group
p_values_shapiro = []
for group, name in zip(
    [normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1],
    ["Control", "Stimulation", "0.1MM Treatment", "1MM Treatment"]
):
    stat, p = stats.shapiro(group)
    p_values_shapiro.append(p)
    print("\n"f'Shapiro-Wilk Test for {name}: Statistic={stat}, p-value={p}')
    if p > 0.05:
        print(f'{name}: Data is normally distributed')
    else:
        print(f'{name}: Data is NOT normally distributed')

# Homogeneity of variances test (Levene's test)
stat, p_levene = stats.levene(
    normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1
)
print("\n"f'Levene’s Test - Statistic: {stat}, p-value: {p_levene}')
if p_levene > 0.05:
    print('Variances are equal (can use ANOVA)')
else:
    print('Variances are NOT equal (use Welch ANOVA or non-parametric test)')

# Choosing the correct statistical test
if all(p > 0.05 for p in p_values_shapiro) and p_levene > 0.05:
    # Parametric test (ANOVA)
    stat, p_anova = stats.f_oneway(
        normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1
    )
    print("\n"f'ANOVA - Statistic: {stat}, p-value: {p_anova}')
elif all(p > 0.05 for p in p_values_shapiro) and p_levene <= 0.05:
    # Welch ANOVA (for unequal variances)
    stat, p_welch = stats.ttest_ind(normalized_stimulation, normalized_treatment_01, equal_var=False)
    print(f'Welch’s t-test - Statistic: {stat}, p-value: {p_welch}')
else:
    # Non-parametric test (Kruskal-Wallis)
    stat, p_kw = stats.kruskal(
        normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1
    )
    print(f'Kruskal-Wallis Test - Statistic: {stat}, p-value: {p_kw}')

# Tukey's test
data = normalized_control + normalized_stimulation + normalized_treatment_01 + normalized_treatment_1
group = ['Control']*len(normalized_control) + ['Stimulation']*len(normalized_stimulation) + ['0.1MM Treatment']*len(normalized_treatment_01) + ['1MM Treatment']*len(normalized_treatment_1)
df = pd.DataFrame({'value': data, 'group': group})

tukey_result = pairwise_tukeyhsd(df['value'], df['group'], alpha=0.05)
print(tukey_result)

# plotting
group_order = ['Control', 'Stimulation', '0.1MM Treatment', '1MM Treatment']
palette = ['#FF0099', '#CC0099', '#990066', '#660066']

sns.set_theme(style="whitegrid")

plt.figure(figsize=(10, 6))
ax = sns.boxplot(x='group', y='value', data=df, order=group_order, palette=palette, linewidth=1.3)

# Add means to the plot
means = df.groupby('group')['value'].mean().reindex(group_order)
for i, mean in enumerate(means):
    ax.plot(i, mean, marker='o', color='white', markersize=8, markeredgewidth=2, markeredgecolor='#4A235A')

ax.set_title('Phosphorylation of ERK', fontsize=18, fontweight='bold', color='#4A235A')
ax.set_ylabel('ERK / Beta Actin', fontsize=14, color='#4A235A')
ax.set_xlabel('')
ax.tick_params(axis='x', labelsize=12, colors='#4A235A')
ax.tick_params(axis='y', labelsize=12, colors='#4A235A')
sns.despine()

# Add significance markers
significant_pairs = tukey_result.summary().data[1:]
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