#!/bin/bash
#BSUB -J poissonmulti[2-3]
#BSUB -o logs/%J_mpi%I_multi.out.txt
#BSUB -e logs/%J_mpi%I_multi.err.txt
#BSUB -q hpc
#BSUB -n 48
#BSUB -R "rusage[mem=2G]"
#BSUB -R "span[ptile=24]"
#BSUB -W 3:00
##BSUB -env MPIVER
#BSUB -R "select[model == XeonE5_2650v4]"

unset LSB_AFFINITY_HOSTFILE

module load studio
module load "mpi/${MPIVER:-3.1.5-gcc-9.2.0}"

# BASH_LIB=$(pwd)
# ROOT_LIB=$BASH_LIB/..
# SOURCE_LIB=$ROOT_LIB/Source

# H=(100 200 500 1000)
# W=(100 200 500)
# N=(3 4 8 16 24)
H=(200)
W=(200)
N=(48)

iter_max=2000
tolerance=0
start_T=15
# output_type=0
output_type=4

INSTANCE="mpi${LSB_JOBINDEX}"

mkdir -p "logs"

echo "$INSTANCE"
cd "$INSTANCE"

#make clean
make -B
for ((h=0;h<${#H[@]};h++));
do
    for ((w=0;w<${#W[@]};w++));
    do
        for ((n=0;n<${#N[@]};n++));
        do
            #printf "H: %d, W: %d, N: %d\n" "${H[$h]}" "${W[$w]}" "${N[$n]}"
            echo "mpirun --map-by node --bind-to socket --rank-by socket --report-bindings -n ${N[$n]} ./poisson ${H[$h]} ${W[$w]} $iter_max $tolerance $start_T $output_type"
            mpirun --map-by node --bind-to socket --rank-by socket --report-bindings -n "${N[$n]}" ./poisson "${H[$h]}" "${W[$w]}" $iter_max $tolerance $start_T $output_type
        done
    done
done
