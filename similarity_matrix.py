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

file_names = [ground_truth_file_path]
for pipeline_type in pipeline_types:
    file_names.append(VCF_folder_path + pipeline_type + "/output3_filtered.recode.vcf")

similarity_matrix = np.zeros((len(file_names), len(file_names)))

for i in range(len(file_names)):
    for j in range(i, len(file_names)):
        filtered1_ids = []
        with pysam.VariantFile(file_names[i], "r") as filtered1_vcf:
            for variant in filtered1_vcf:
                filtered1_ids.append(str(variant).split("\t")[1])

        filtered2_ids = []
        with pysam.VariantFile(file_names[j], "r") as filtered2_vcf:
            for variant in filtered2_vcf:
                filtered2_ids.append(str(variant).split("\t")[1])
        
        intersection_size = len(set(filtered1_ids) & set(filtered2_ids))
        union_size = len(set(filtered1_ids) | set(filtered2_ids))
        similarity_matrix[i, j] = intersection_size / union_size
        similarity_matrix[j, i] = similarity_matrix[i, j]

sn.clustermap(data=similarity_matrix, method='average', annot=True, cmap='coolwarm', fmt=".2f", figsize=(10, 5), xticklabels=["ground-truth", *pipeline_types], yticklabels=["ground-truth", *pipeline_types])
plt.title('Clustermap of Pipelines')

plt.savefig("Heatmaps/similarity_matrix.png")

plt.show()
