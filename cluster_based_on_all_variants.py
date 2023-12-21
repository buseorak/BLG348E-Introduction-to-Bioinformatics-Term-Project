import pysam
import seaborn as sn 
import matplotlib.pyplot as plt
import numpy as np

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

all_ids = []
for file_name in files:
    with pysam.VariantFile(file_name, "r") as file:
        for variant in file:
            variant_id = str(variant).split("\t")[1]
            if variant_id not in all_ids:
                all_ids.append(variant_id)

heatmap_data = np.zeros((len(all_ids), len(files)))

for file_name_i in range(len(files)):
    file_name = files[file_name_i]
    with pysam.VariantFile(file_name, "r") as file:
        for variant in file:
            variant_id = str(variant).split("\t")[1]
            if variant_id in all_ids:
                variant_index = all_ids.index(variant_id)
                heatmap_data[variant_index, file_name_i] = 1

sn.clustermap(data=heatmap_data, method='average', cmap='coolwarm', fmt=".4f", figsize=(15, 8), xticklabels=["ground-truth", *pipeline_types])
plt.title('Clustering Pipelines Based on All Variants')

plt.savefig("Heatmaps/clustering_based_on_all_variants.png")

plt.show()
