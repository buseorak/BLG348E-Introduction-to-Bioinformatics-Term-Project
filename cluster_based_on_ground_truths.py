import pysam
import seaborn as sn 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

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

ground_truth_file_path = "/home/buse/Desktop/ITU/Courses/Bioinformatics/term_project/BLG348E-Introduction-to-Bioinformatics-Term-Project/hc_bed_filtered.recode.vcf"

files = [ground_truth_file_path]
for pipeline_type in pipeline_types:
    files.append(VCF_folder_path + pipeline_type + "/output3_filtered.recode.vcf")

ground_truth_ids = []
with pysam.VariantFile(ground_truth_file_path, "r") as ground_truth_file:
    for variant in ground_truth_file:
        ground_truth_ids.append(str(variant).split("\t")[1])

heatmap_data = np.zeros((len(ground_truth_ids), len(files)))

for file_name_i in range(len(files)):
    file_name = files[file_name_i]
    with pysam.VariantFile(file_name, "r") as file:
        for variant in file:
            variant_id = str(variant).split("\t")[1]
            if variant_id in ground_truth_ids:
                variant_index = ground_truth_ids.index(variant_id)
                heatmap_data[variant_index, file_name_i] = 1

colors = ['#82BAED', '#13599C']
custom_cmap = ListedColormap(colors)
sn.clustermap(data=heatmap_data, method='average', cmap=custom_cmap, fmt=".4f", figsize=(15, 8), xticklabels=["ground-truth", *pipeline_types])
plt.title('Clustering Pipelines Based on Ground Truths')

plt.savefig("Heatmaps/clustering_based_on_ground_truths.png")
plt.show()
