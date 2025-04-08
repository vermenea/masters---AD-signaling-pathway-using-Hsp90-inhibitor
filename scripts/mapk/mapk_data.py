import pandas as pd

data = {
    "control": [221469, 221500, 221480],
    "stimulation": [211094, 211120, 211110],
    "0.1MM_treatment": [165197, 165210, 165205],
    "1MM_treatment": [118330, 118350, 118340]
}

# Convert data dictionary to DataFrame
mapk_data = pd.DataFrame(data)

# Calculate means
mapk_means = mapk_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(mapk_means, columns=['mean']).transpose()

print("\nMAPK Data:")
print(mapk_data)
print("\nMeans DataFrame:")
print(means_df)