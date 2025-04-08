#imports
import sys
sys.path.append('/Users/vermenea/Documents/biologia/masters/analysis/scripts')

from beta_actin_data import beta_actin_data
from pstat3_data import pstat3_data
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

#preparing data
control_beta_actin = beta_actin_data["control"][0]
stimulation_beta_actin = beta_actin_data["stimulation"][0]
treatment_01_beta_actin = beta_actin_data["0.1MM_treatment"][0]
treatment_1_beta_actin = beta_actin_data["1MM_treatment"][0]

control_ps3 = pstat3_data["control"][0]
stimulation_ps3 = pstat3_data["stimulation"][0]
treatment_01_ps3 = pstat3_data["0.1MM_treatment"][0]
treatment_1_ps3 = pstat3_data["1MM_treatment"][0]

normalized_control = [control_beta_actin / control_ps3]
normalized_stimulation = [stimulation_beta_actin / stimulation_ps3]
normalized_treatment_01 = [treatment_01_beta_actin / treatment_01_ps3]
normalized_treatment_1 = [treatment_1_beta_actin / treatment_1_ps3]

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

#plotting
labels = ['Control', 'Stimulation', '0.1MM Treatment', '1MM Treatment']
normalized_values = [normalized_control[0], normalized_stimulation[0], normalized_treatment_01[0], normalized_treatment_1[0]]

plt.figure(figsize=(10, 6))
plt.bar(labels, normalized_values, color=['#FF0099', '#CC0099', '#990066', '#660066'])
plt.xlabel('Normalized pSTAT3 Levels with Beta-actin')
plt.ylabel('Normalized pSTAT3')
plt.title('Phosphorylation of STAT3')
plt.show()