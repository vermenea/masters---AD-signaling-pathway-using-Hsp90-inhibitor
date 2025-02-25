#imports
import sys
sys.path.append('/Users/vermenea/Documents/biologia/masters/analysis/scripts')

from beta_actin_data import beta_actin_data
from erk_data import data as erk_data
import matplotlib.pyplot as plt

#preparing data
control_beta_actin = beta_actin_data["control"][0]
stimulation_beta_actin = beta_actin_data["stimulation"][0]
treatment_01_beta_actin = beta_actin_data["0.1MM_treatment"][0]
treatment_1_beta_actin = beta_actin_data["1MM_treatment"][0]

control_erk = erk_data["control"][0]
stimulation_erk = erk_data["stimulation"][0]
treatment_01_erk = erk_data["0.1MM_treatment"][0]
treatment_1_erk = erk_data["1MM_treatment"][0]

normalized_control = control_beta_actin / control_erk
normalized_stimulation = stimulation_beta_actin / stimulation_erk
normalized_treatment_01 = treatment_01_beta_actin / treatment_01_erk
normalized_treatment_1 = treatment_1_beta_actin / treatment_1_erk

#plotting
labels = ['Control', 'Stimulation', '0.1MM Treatment', '1MM Treatment']
normalized_values = [normalized_control, normalized_stimulation, normalized_treatment_01, normalized_treatment_1]

plt.figure(figsize=(10, 6))
plt.bar(labels, normalized_values, color=['#FF0099', '#CC0099', '#990066', '#660066'])
plt.xlabel('Normalized ERK Levels with Beta-actin')
plt.ylabel('Normalized ERK')
plt.title('Phosphorylation of ERK')
plt.show()