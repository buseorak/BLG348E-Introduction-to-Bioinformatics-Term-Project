import pysam
import seaborn as sn 
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting

VCF_folder_path = "/home/buse/Desktop/ITU/Courses/Bioinformatics/term_project/BLG348E-Introduction-to-Bioinformatics-Term-Project/VCF/"

pipeline_types = ["norecal-bowtie-mutect",
                  "norecal-bowtie-somaticsniper",
                  "norecal-bowtie-strelka",
                  "norecal-bwa-mutect",
                  "norecal-bwa-somaticsniper",
                  "norecal-bwa-strelka",
                  "recal-bowtie-mutect",
                  "recal-bowtie-somaticsniper",
                  "recal-bowtie-strelka",
                  "recal-bwa-mutect",
                  "recal-bwa-somaticsniper",
                  "recal-bwa-strelka"]

files = []
for pipeline_type in pipeline_types:
    files.append(VCF_folder_path + pipeline_type + "/output3_filtered.recode.vcf")

all_ids = []
for file_name in files:
    with pysam.VariantFile(file_name, "r") as file:
        for variant in file:
            variant_id = str(variant).split("\t")[1]
            if variant_id not in all_ids:
                all_ids.append(variant_id)

data = np.zeros((len(files), len(all_ids)))

for file_name_i in range(len(files)):
    file_name = files[file_name_i]
    with pysam.VariantFile(file_name, "r") as file:
        for variant in file:
            variant_id = str(variant).split("\t")[1]
            if variant_id in all_ids:
                variant_index = all_ids.index(variant_id)
                data[file_name_i, variant_index] = 1

scaler = StandardScaler()
data_normalized = scaler.fit_transform(data)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(data_normalized)

# Plot the results without class labels
plt.figure(figsize=(8, 6))
plt.scatter(X_pca[:, 0], X_pca[:, 1])

plt.title('PCA on Variants for Pipelines')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')

for i in range(X_pca.shape[0]):
    plt.text(X_pca[i, 0], X_pca[i, 1], pipeline_types[i], color='black', fontsize=8, ha='right', va='bottom')

plt.savefig("PCA/pca_2_components.png")

plt.show()