README
================

In this document, we describe how to produce simulation results included
in our manuscript for supporting our findings.

# Data

We provided an R program `Simulations_data.R` to generate data used in
our simulation studies. Those simulated data will be saved in a
sub-folder named `data` under the working directory.

# R scripts

- `FUNs.R` (self-defined functions)
- `FUNs_cpp.cpp`
- `Simulations_data.R`
- `Simulations.R`
- `Sum_res.R`

# Bash shell scripts

- `sim_data.sh`
- `sim.sh`
- `sum_res.sh`

# Intermediate results

In the simulation study, a total of 10 simulation scenarios were
considered, and within each scenario 100 replicates were generated.
Since one simulation takes about 5 to 6 minutes to run on a computer
running Apple M2 pro chip with 32 Gb memory, the R code provided here
requires to save intermediate results in the sub-folder named
`inter_res` under the working directory. Those saved intermediate
results will be used to generate the Table 1 included in the manuscript.

# Instructions for use

There are two options to run those provided R scripts: (1) RStudio or
(2) Bash Shell scripts. When running R codes in RStudio, please create a
new R project and make sure R scripts are stored at the same folder
where the created R project `*.Rproj` is; when running R codes using
provided bash shell scripts, please set the working directory to the
folder where R scripts and bash shell scripts are.

Run-time was approximated using a computer with Apple M2 pro chip and 32
Gb memory.

Steps to produce results presented in Table 1 and discussed in Section 3
are given below.

**Step 1**. Generate data

Run `Simulations_data.R` in RStudio; or execute the bash script
`sim_data.sh` using the command `bash sim_data.sh` in terminal.

`dat_inde_ns_st.rda`, `dat_pos_s_ns_st.rda`, `dat_pos_m_ns_st.rda`,
`dat_nega_s_ns_st.rda`, `dat_nega_m_ns_st.rda` will be saved to the
sub-folder named `data` under the working directory.

**Step 2**. Fit different models to simulated data

Run `Simulations.R` in RStudio (arguments `isim`, `sceindex`, and
`j.index` have to be specified); or execute the bash script `sim.sh`.
One simulation takes about 5-6 minutes to run, we have a total of 1000
simulations (10 simulation scenarios, and 100 simulations were conducted
within each scenario). Those intermediate results will be saved to the
sub-directory named `inter_res` under the working directory.

**Step 3**. Summarize simulation results to generate the Table 1 in the
manuscript.

Run `Sum_res.R` in RStudio; or execute the bash script `sum_res.sh`.
`Table1.xlsx` will be saved to the sub-folder named `table`.
