## creating an epigentic clock in daphnia magna

data provided by [Liu *et al.* (2025)](https://doi.org/10.1186/s13072-025-00580-y) and [Hearn *et al.* (2021)](https://doi.org/10.1186/s13072-020-00379-z)

## preprocessing genome data

using ALICE3: University of Leicester HPC

#### quality control 

[multiqc.sh](docs/multiqc.sh) : trimming and [MultiQC](https://seqera.io/multiqc/) to assess read quality

#### alignment

preparing reference genomes with [Bismark](https://github.com/FelixKrueger/Bismark)

`/bin/Bismark-0.22.3/bismark_genome_prepatation genome_folder`


workflow using [snakemake](https://snakemake.github.io/)

`conda create --name snakemake_env`

`conda activate snakemake_env `

`conda install bioconda::snakemake`

`conda install bioconda::bowtie2`

`conda install bioconda::samtools`

[config.yaml](docs/config.yaml) : example of required configuration file

[Snakefile](docs/Snakefile) : alignment workflow

[snakemake.slm](docs/snakemake.slm) : slurm script to execute Snakefile

## differential methylation

specific to [Hearn *et al.* (2021)](https://doi.org/10.1186/s13072-020-00379-z) data

[Hearn_pca.R](docs/Hearn_pca) : PCA of methylation patterns

[Hearn_DSS_inputs.R](docs/Hearn_DSS_inputs) : processing data in preparation for DSS

specific to C32 strain 

[C32_DSS.R](docs/C32_DSS.R) : DSS analysis to identify significant age-associated CpGs and methylation percentages

## epigenetic age / predicted age

[C32_epiclock.R](docs/C32_epiclock.R) : epigenetic clock and extraction of clock loci (edited from code provided by [Eamonn Mallon](https://github.com/EamonnMallon/))

[C32_loci_age_pred.R](docs/C32_loci_age_pred.R) : predicting age using methylation levels and Bham2 clock loci
