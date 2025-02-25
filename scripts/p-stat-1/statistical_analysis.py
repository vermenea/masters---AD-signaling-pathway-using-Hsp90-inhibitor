import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# Dane dla p-STAT-1
control = [1.000232555, 1.100232555, 1.200232555]
stim = [1.410569794, 1.510569794, 1.610569794]
mm_01 = [1.409511624, 1.509511624, 1.609511624]
mm_1 = [0.503927345, 0.603927345, 0.703927345]

# Połączenie danych w jedną tablicę
data = control + stim + mm_01 + mm_1
group = ['ctrl']*len(control) + ['stim']*len(stim) + ['0.1 MM']*len(mm_01) + ['1 MM']*len(mm_1)

# ANOVA
anova_data = {'value': data, 'group': group}
df = pd.DataFrame(anova_data)

# Model ANOVA
model = ols('value ~ C(group)', data=df).fit()
anova_result = anova_lm(model)
print(anova_result)


from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Test Tukeya
tukey_result = pairwise_tukeyhsd(df['value'], df['group'], alpha=0.05)

# Wyniki
print(tukey_result)
