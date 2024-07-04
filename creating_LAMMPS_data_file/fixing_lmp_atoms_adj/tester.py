# Function to calculate the sum of the fourth column and count names with cumulative charges
def analyze_data(lines):
    fourth_column_sum = 0.0
    name_counts = {}
    name_cumulative_charge = {}

    for line in lines:
        parts = line.split()
        charge = float(parts[3])
        name = parts[-1]

        # Sum the charges in the fourth column
        fourth_column_sum += charge

        # Count occurrences of each name
        if name in name_counts:
            name_counts[name] += 1
            name_cumulative_charge[name] += charge
        else:
            name_counts[name] = 1
            name_cumulative_charge[name] = charge

    return fourth_column_sum, name_counts, name_cumulative_charge

# Read the text file
file_path = 'original.txt'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Analyze the data
fourth_column_sum, name_counts, name_cumulative_charge = analyze_data(lines)

# Print the results
print(f"Sum of the fourth column: {fourth_column_sum}")
print("\nName Counts and Cumulative Charges:")
for name in name_counts:
    print(f"{name}: Count = {name_counts[name]}, Cumulative Charge = {name_cumulative_charge[name]}")

