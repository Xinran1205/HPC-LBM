
# Introduction

The Lattice Boltzmann Method (LBM) is a powerful and flexible approach to simulate complex fluid flows. However, to fully exploit its capabilities, optimizing the underlying code for performance is essential. This project details my comprehensive optimization process, which spans several stages:

### Serial Optimization

I began with serial optimization techniques to lay a solid foundation for performance improvements. Techniques such as loop fusion and pointer swapping were implemented to streamline the execution flow and reduce overhead.

### Vectorization

Recognizing the importance of memory access patterns for performance, I revised the data structures used within the code. Additionally, memory alignment techniques were applied to enhance data access speed and facilitate vectorization by the compiler.

### Parallel Acceleration with OpenMP

The culmination of my optimization efforts involved leveraging OpenMP, a parallel programming model, to accelerate computation. By scaling the code from a single core to utilizing 28 cores, I achieved significant gains in performance and efficiency, making my LBM simulations more practical for complex scenarios.

### Distributed memory parallelism with MPI

This part extends the optimization of the LBM code by employing the Message Passing Interface (MPI). Running on the BlueCrystal supercomputer, the optimized code utilizes four nodes, each equipped with 28 cores. This segment builds upon the previously optimized serial code. In contrast to the shared memory model used in OpenMP, MPI employs a distributed memory model suitable for multi-node computing environments. Each node has its own independent physical memory, and data exchange is explicitly conducted through network or other communication methods. In this project phase, MPI alone was utilized to manage inter-node communications, achieving performance levels that significantly enhance the practicality of LBM simulations for even more complex scenarios.

## Compiling and running

To compile type `make`. Editing the values for `CC` and `CFLAGS` in the Makefile can be used to enable different compiler options or use a different compiler. These can also be passed on the command line:

    $ make CFLAGS="-O3 -fopenmp -DDEBUG"

Input parameter and obstacle files are all specified on the command line of the `d2q9-bgk` executable.

Usage:

    $ ./d2q9-bgk <paramfile> <obstaclefile>
eg:

    $ ./d2q9-bgk input_256x256.params obstacles_256x256.dat

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


# Visualisation

You can view the final state of the simulation by creating a .png image file using a provided Gnuplot script:

    $ gnuplot final_state.plt
