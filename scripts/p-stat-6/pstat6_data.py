import pandas as pd

#data after densynometry and before normalization with beta-actin
pstat6_data = pd.DataFrame({
  "control": [331431, 330000, 333000],  
  "stimulation": [156147, 157000, 155500],
  "0.1MM_treatment": [156127, 157500, 155800],
  "1MM_treatment": [211121, 212000, 210500]
})

# Calculate means
pstat6_means = pstat6_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(pstat6_means, columns=['mean']).transpose()

print("\nPSTAT6 Data:")
print(pstat6_data)
print("\nMeans DataFrame:")
print(means_df)