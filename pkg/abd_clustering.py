from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

def plot_kmeans(values, k_start=1, k_end=20):
    scaler = StandardScaler()
    scaler.fit(values)
    scaled_values = scaler.transform(values).T
    inertia = []
    for k in range(k_start,k_end+1):
        km = KMeans(n_clusters=k)
        km.fit(scaled_values)
        inertia.append(km.inertia_)
    plt.plot(range(k_start,k_end+1), inertia,'o-')
    plt.xlabel('Clusters K')
    plt.ylabel('Distances to Cluster Centers')
    plt.title('K-Means Fit')
    plt.xticks(range(k_start,k_end+1))
    plt.tight_layout()
    plt.show()
    
def kmeans_clustering(values, k):
    scaler = StandardScaler()
    scaler.fit(values)
    scaled_values = scaler.transform(values).T
    km = KMeans(n_clusters=k)
    km.fit(scaled_values)
    labels = km.labels_
    return [i+1 for i in labels]

def abd_label(dataset, k):
    dataset.abd_labels = kmeans_clustering(dataset.values, k)
    return dataset