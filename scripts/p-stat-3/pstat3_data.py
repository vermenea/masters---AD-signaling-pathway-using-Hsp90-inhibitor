data = {
   
}

import pandas as pd

#data after densynometry and before normalization with beta-actin
pstat3_data = pd.DataFrame({
   "control": [214904, 215000, 214950],
   "stimulation": [204668, 204700, 204680],
   "0.1MM_treatment": [204672, 204690, 204675],
   "1MM_treatment": [204641, 204660, 204650]
})

# Calculate means
pstat3_means = pstat3_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(pstat3_means, columns=['mean']).transpose()

print("\nPSTAT3 Data:")
print(pstat3_data)
print("\nMeans DataFrame:")
print(means_df)

# Export pstat3_data and means_df for use in other scripts
__all__ = ['pstat3_data', 'means_df']