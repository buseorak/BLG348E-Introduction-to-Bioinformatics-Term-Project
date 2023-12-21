import pysam

pipeline_types = ["norecal-bowtie-mutect",
                  "norecal-bowtie-strelka",
                  "norecal-bwa-mutect",
                  "norecal-bwa-strelka",
                  "recal-bowtie-mutect",
                  "recal-bowtie-strelka",
                  "recal-bwa-mutect",
                  "recal-bwa-strelka"]

for pipeline_type in pipeline_types:
    file_path = "/home/buse/Desktop/ITU/Courses/Bioinformatics/term_project/BLG348E-Introduction-to-Bioinformatics-Term-Project/VCF/" + pipeline_type + "/output3.recode.vcf"
    output_file_path = "/home/buse/Desktop/ITU/Courses/Bioinformatics/term_project/BLG348E-Introduction-to-Bioinformatics-Term-Project/VCF/" + pipeline_type + "/output3_filtered.recode.vcf"

    with pysam.VariantFile(file_path, "r") as vcf_in:
        with pysam.VariantFile(output_file_path, "w", header=vcf_in.header) as vcf_out:
            for variant in vcf_in:
                lst = str(variant).split("\t")
                if lst[6] == "PASS" or lst[6] == ".":
                    vcf_out.write(variant)
