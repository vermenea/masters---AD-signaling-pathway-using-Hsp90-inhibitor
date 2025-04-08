import pandas as pd

data = {
    "control": [221469, 221500, 221480],
    "stimulation": [211094, 211120, 211110],
    "0.1MM_treatment": [165197, 165210, 165205],
    "1MM_treatment": [118330, 118350, 118340]
}

# Convert data dictionary to DataFrame
erk_data = pd.DataFrame(data)

# Calculate means
erk_means = erk_data.mean()

# Create DataFrame from means
means_df = pd.DataFrame(erk_means, columns=['mean']).transpose()

print("\nERK Data:")
print(erk_data)
print("\nMeans DataFrame:")
print(means_df)