import pysam
import seaborn as sn 
import matplotlib.pyplot as plt
import seaborn as sn 

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
ground_truth_ids = []
with pysam.VariantFile(ground_truth_file_path, "r") as ground_truth_vcf:
    for variant in ground_truth_vcf:
        ground_truth_ids.append(str(variant).split("\t")[1])

# 4 Metrics (Accuracy, Recall, Precision, F1-Score) on columns, 12 pipelines on rows.
performance_metrics = []

for pipeline_type in pipeline_types:
    metrics = []

    filtered_file_path = VCF_folder_path + pipeline_type + "/output3_filtered.recode.vcf"
    unfiltered_file_path = VCF_folder_path + pipeline_type + "/output3.recode.vcf"

    filtered_ids = []
    with pysam.VariantFile(filtered_file_path, "r") as filtered_vcf:
        for variant in filtered_vcf:
            filtered_ids.append(str(variant).split("\t")[1])

    counter_tn = 0
    with pysam.VariantFile(unfiltered_file_path, "r") as unfiltered_vcf:
        for variant in unfiltered_vcf:
            variant_id = str(variant).split("\t")[1]
            if (variant_id not in filtered_ids) and (variant_id not in ground_truth_ids):
                counter_tn += 1

    counter_tp = 0
    counter_fp = 0
    for filtered_id in filtered_ids:
        if filtered_id in ground_truth_ids:
            counter_tp += 1
        else:
            counter_fp += 1

    counter_fn = 0
    for ground_truth_id in ground_truth_ids:
        if ground_truth_id not in filtered_ids:
            counter_fn += 1

    accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
    recall = counter_tp / (counter_tp + counter_fn)
    precision = counter_tp / (counter_tp + counter_fp)
    f1_score = 2 * precision * recall / (precision + recall)
    metrics.append(accuracy)
    metrics.append(recall)
    metrics.append(precision)
    metrics.append(f1_score)
    performance_metrics.append(metrics)

    confusion_matrix = [[counter_tp, counter_fn], [counter_fp, counter_tn]]
    #sn.heatmap(data=confusion_matrix, annot=True, fmt=".0f")
    #plt.title("Confusion Matrix for " + pipeline_type)
    #plt.xlabel("Predicted Value")
    #plt.ylabel("Actual Value")
    
    #plt.savefig("Performance_Metrics/" + pipeline_type + "-confusion-matrix.png")

    #plt.show()

fig, ax = plt.subplots(figsize=(20, 10))
sn.heatmap(data=performance_metrics, annot=True, fmt=".6f", annot_kws={"size": 16}, xticklabels=["Accuracy", "Recall", "Precision", "F1-Score"], yticklabels=pipeline_types)
plt.title("Table of Performance Metrics")
plt.xlabel("Performance Metrics")
plt.ylabel("Pipelines")
plt.savefig("Performance_Metrics/performance_metrics.png")
plt.show()
