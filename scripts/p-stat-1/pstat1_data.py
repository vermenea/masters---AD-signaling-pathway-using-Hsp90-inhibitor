import pandas as pd

#data after densynometry and before normalization with beta-actin
pstat1_data = pd.DataFrame({
   "control": [221431, 220000, 223000],
    "stimulation": [105090, 106000, 104500],
    "0.1MM_treatment": [105127, 106500, 104800],
    "1MM_treatment": [211038, 212000, 210500]
})

# Calculate means
pstat1_means = pstat1_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(pstat1_means, columns=['mean']).transpose()

print("\nPSTAT1 Data:")
print(pstat1_data)
print("\nMeans DataFrame:")
print(means_df)