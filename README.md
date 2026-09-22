**Variant Segregation Workflow**

This Snakemake workflow processes genomic data (PLINK2 format) to identify segregating variants in families across multiple inheritance models (de novo, recessive, dominant, and compound heterozygous). 

🛠 Prerequisites and Dependencies\
The pipeline relies on several containerized tools via Docker/Singularity, as well as local binaries. It requires snakemake version >=9
  
📂 Expected Input Directory Structure\
meta/segregation_report_ped.tsv: TSV file containing the PED files with the header.\
meta/segregation_report_fam.ped: Standard PED file for slivar inheritance logic.\
meta/samples_to_keep_family.txt: List of family samples to extract from the larger dataset.\
input/pfiles/: Directory containing input PLINK2 genotype files (.pgen, .pvar, .psam) per chromosome.

📊 Expected Outputs\
Upon successful completion, the workflow generates the following finalized files:\
Family Segregation Reports: families/reports/{family_id}_segregation_variants.csv\
The aggregated reports for each family containing segregating variants within the family filtered by gnomAD genome frequency, sorted by impact and inheritance mode.

🚀 Usage\
To execute the workflow, run Snakemake from the directory.\
snakemake --sdm conda apptainer env-modules --snakefile gp2_family_variant_segregation.smk
