import pandas as pd

#data after densynometry and before normalization with beta-actin
pstat6_data = pd.DataFrame({
  "control": [331431, 331431, 331431],
    "stimulation": [156147, 156147, 156147],
    "0.1MM_treatment": [156147, 156147, 156147],
    "1MM_treatment": [211121, 211121, 211121]
})

# Calculate means
pstat6_means = pstat6_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(pstat6_means, columns=['mean']).transpose()

print("\nPSTAT6 Data:")
print(pstat6_data)
print("\nMeans DataFrame:")
print(means_df)