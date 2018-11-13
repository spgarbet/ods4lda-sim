It is considered good practice to follow the principals of cohesion in coding,
i.e., put like things with like things and each block has a purpose. While 
someone philosophical, there are some complexity arguments that support
this approach. With that in mind, the following directory outlines what each
files purpose is in this project

| Filename                    | Purpose                                   |
| --------------------------- | ----------------------------------------- |
| acml.R                      | ACML fitting functions                    |
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


The main loop function `simulation(run, count, save.raw=FALSE)` runs the model
setting the random seed to the `run` variable and `count` to the indexed parameter
setup. The `save.raw` flag will save the raw model fits if desired. Output goes 
into the `output` directory in the form of `output/run-<run>-<count>.RData` and
if raw was requested `output/raw-<run>-<count>.RData` likewise.

All output to the terminal is controlled by the `progress` function at the top
of the `run.R` code. This is helpful if one decides to eliminate all output, 
this function can just be set to do nothing without having to edit all occurances
of output in the code. This is a good example of functional programming's benefits.

When run on ACCRE the output of each job goes into the `status` directory as
`status/job<run>.out`. Also the ACCRE SBATCH command controls the run numbers
and is the most commonly edited line. For example, it's commonly set to execute
2000 runs from 1-2000.