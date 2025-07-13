## creating a robust epigentic clock in daphnia magna

data provided by [Liu *et al.*](https://doi.org/10.1186/s13072-025-00580-y) and [Hearn *et al.*](https://doi.org/10.1186/s13072-020-00379-z)

## preprocessing genome data
using ALICE3: University of Leicester HPC


### qualtiy control 
[multiqc.sh](docs/multiqc.sh) -- trimming and [MultiQC](https://seqera.io/multiqc/) to assess read quality

### read alignment
preparing reference genomes with [Bismark](https://github.com/FelixKrueger/Bismark)

`/bin/Bismark-0.22.3/bismark_genome_prepatation genome_folder`

`/bin/Bismark-0.22.3/bismark_genome_prepatation lambda_genome`

[snakemake](https://snakemake.github.io/) workflow

`conda create --name snakemake_env`

`conda activate snakemake_env `


[Snakefile](docs/Snakefile) -- alignment workflow
[snakemake.slm](docs/snakemake.slm) -- slurm script to execute Snakefile
