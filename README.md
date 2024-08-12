
# Lattice Boltzmann Method (LBM) Optimization Project
- zh_CN [简体中文](/README.zh_CN.md)
## Introduction
The Lattice Boltzmann Method (LBM) is a powerful and flexible method for simulating complex fluid flows. However, to fully utilize its capabilities, optimizing the underlying code to enhance performance is crucial. This project details my comprehensive optimization process, which covers several stages:

## Optimization Stages

### Serial Optimization
I began with serial optimization techniques to lay a solid foundation for performance improvement. I implemented loop fusion and pointer swapping techniques to streamline the execution process and reduce overhead.

### Vectorization Optimization
During the optimization process, I analyzed the impact of memory access patterns on overall performance. To enhance data access speed and facilitate effective compiler vectorization, I first adjusted the data structures used in the code to better suit vector operations. I then applied memory alignment techniques, which optimized the data loading process. To further improve vectorization efficiency, I experimented with various compilers and their different versions, applying specific compiler flags. By analyzing the optimization reports provided by each compiler, I was able to select the most suitable compilation strategy for our project, thereby maximizing the execution efficiency and CPU parallel processing capability of the code.

### Parallel Acceleration Using OpenMP
My optimization efforts ultimately utilized OpenMP, a parallel programming model, to accelerate computation. By expanding the code from a single core to 28 cores, I achieved significant improvements in performance and efficiency, making my LBM simulations more suitable for complex scenarios.

### Distributed Memory Parallelization Using MPI
This part of the project extended the optimization of LBM code using the Message Passing Interface (MPI). The optimized code runs on the BlueCrystal supercomputer, utilizing 4 nodes, each equipped with 28 cores. Unlike the shared memory model used in OpenMP, MPI adopts a distributed memory model suitable for multi-node computing environments. Each node has its own independent physical memory, and data exchanges are explicitly conducted through networks or other communication methods. I implemented data allocation and initialization for each process and adopted load balancing strategies. Additionally, through the Halo Exchange strategy, I optimized inter-node data communication, optimizing memory storage while ensuring data accuracy.

## Results
<img src="/pic/1.png" alt="scalability" width="450" height="350">
<img src="/pic/2.png" alt="speedup" width="500" height="350">
## Compiling and running

### Compiling and running the vectorized or OpenMP version of the code

Type `make`. Editing the values for `CC` and `CFLAGS` in the Makefile can be used to enable different compiler options or use a different compiler. These can also be passed on the command line:

    $ make CFLAGS="-O3 -fopenmp -DDEBUG"

Input parameter and obstacle files are all specified on the command line of the `d2q9-bgk` executable.

Usage:

    $ ./d2q9-bgk <paramfile> <obstaclefile>
eg:

    $ ./d2q9-bgk input_256x256.params obstacles_256x256.dat

### Compiling and running MPI code locally

To compile the MPI program, use the following command:

    $ mpicc -o result d2q9-bgk.c -lm

To run the MPI program with 4 processes, use this command:

    $ mpirun -np 4 ./result input_128x128.params obstacles_128x128.dat

In this example, `4` represents the number of processes to be used.

## Checking results

An automated result checking function is provided that requires you to load a particular Python module (`module load languages/anaconda2/5.0.1`). Running `make check` will check the output file (average velocities and final state) against some reference results. By default, it should look something like this:

    $ make check
    python check/check.py --ref-av-vels-file=check/128x128.av_vels.dat --ref-final-state-file=check/128x128.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    Total difference in av_vels : 5.270812566515E-11
    Biggest difference (at step 1219) : 1.000241556248E-14
      1.595203170657E-02 vs. 1.595203170658E-02 = 6.3e-11%

    Total difference in final_state : 5.962977334129E-11
    Biggest difference (at coord (6,2)) : 1.000588500943E-14
      3.329122639178E-02 vs. 3.329122639179E-02 = 3e-11%

    Both tests passed!

This script takes both the reference results and the results to check (both average velocities and final state). This is also specified in the makefile and can be changed like the other options:

    $ make check REF_AV_VELS_FILE=check/128x256.av_vels.dat REF_FINAL_STATE_FILE=check/128x256.final_state.dat
    python check/check.py --ref-av-vels-file=check/128x256.av_vels.dat --ref-final-state-file=check/128x256.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    ...

All the options for this script can be examined by passing the --help flag to it.

    $ python check/check.py --help
    usage: check.py [-h] [--tolerance TOLERANCE] --ref-av-vels-file
                    REF_AV_VELS_FILE --ref-final-state-file REF_FINAL_STATE_FILE
    ...


## Running on BlueCrystal Phase 4

When you wish to submit a job to the queuing system on BlueCrystal, you should use the job submission script provided.

    $ sbatch job_submit_d2q9-bgk

This will dispatch a job to the queue, which you can monitor using the
`squeue` command:

    $ squeue -u $USER

When finished, the output from your job will be in a file called
`d2q9-bgk.out`:

    $ less d2q9-bgk.out

If you wish to run a different set of input parameters, you should
modify `job_submit_d2q9-bgk` to update the value assigned to `options`.

