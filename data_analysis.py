import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

# Load your CSV file
data = pd.read_csv("output.csv")

# 1. Cluster Distribution with Progress Bar
cluster_counts = data['cluster'].value_counts()

plt.figure(figsize=(10, 6))
for index, value in tqdm(zip(cluster_counts.index, cluster_counts.values), total=len(cluster_counts), desc="Plotting Cluster Distribution Bars"):
    plt.bar(index, value)
plt.title('Cluster Distribution')
plt.xlabel('Cluster')
plt.ylabel('Number of Data Points')
plt.show()

# 2. Distribution of Representative Points by Cluster
representative_counts = data[data['is_representative'] == True]['cluster'].value_counts()

plt.figure(figsize=(10, 6))
for index, value in tqdm(zip(representative_counts.index, representative_counts.values), total=len(representative_counts), desc="Plotting Representative Points Bars"):
    plt.bar(index, value)
plt.title('Distribution of Representative Points by Cluster')
plt.xlabel('Cluster')
plt.ylabel('Number of Representative Points')
plt.show()

# 3. Average Retention Time by Cluster
avg_retention_time = data.groupby('cluster')['retention_time'].mean()

plt.figure(figsize=(10, 6))
for index, value in tqdm(zip(avg_retention_time.index, avg_retention_time.values), total=len(avg_retention_time), desc="Plotting Retention Time by Cluster"):
    plt.plot(index, value, marker='o')
plt.title('Average Retention Time by Cluster')
plt.xlabel('Cluster')
plt.ylabel('Average Retention Time')
plt.show()

# 4. Distribution of Precursor m/z
mz_distribution = data['precursor_mz']

plt.figure(figsize=(10, 6))
for i in tqdm(range(30), desc="Plotting Precursor m/z Distribution"):
    plt.hist(mz_distribution, bins=30, edgecolor='black')
plt.title('Distribution of Precursor m/z')
plt.xlabel('Precursor m/z')
plt.ylabel('Frequency')
plt.show()

# 5. Distribution of Precursor Charges
precursor_charge_distribution = data['Precursor_charge'].value_counts()

plt.figure(figsize=(10, 6))
for index, value in tqdm(zip(precursor_charge_distribution.index, precursor_charge_distribution.values), total=len(precursor_charge_distribution), desc="Plotting Precursor Charge Bars"):
    plt.bar(index, value)
plt.title('Distribution of Precursor Charges')
plt.xlabel('Precursor Charge')
plt.ylabel('Frequency')
plt.show()
