It is considered good practice to follow the principals of cohesion in coding,
i.e., put like things with like things and each block has a purpose. While 
someone philosophical, there are some complexity arguments that support
this approach. With that in mind, the following directory outlines what each
files purpose is in this project

| Filename                    | Purpose                                   |
| functions.R                 | ACML fitting functions                    |
| generate-data.R             | Simulate data for a given model           |
| impute.R                    | Imputation fitting functions              |
| missing.sh                  | Script to check for missing runs on ACCRE |
| ods4lda.slurm               | Slurm script to run on ACCRE              |
| run.R                       | The main loop / work horse                |
| setup.R                     | Parameters for the runs                   |
| sim-accre.R                 | Simulate using ACCRE (sapply)             |
| sim-local.R                 | Simulate using local CPU (mcapply)        |
| two_phase.R                 | Ran Tao's 2 phase fitting function        |
| report/performance-acml.Rmd | Generate a summary of results             |
| report/figure.R             | Generate a summary figure                 |