#!/bin/bash
#BSUB -J poissonomp
#BSUB -o logs/%J_omp-twonodes.out.txt
#BSUB -e logs/%J_omp-twonodes.err.txt
#BSUB -q hpc
#BSUB -n 48
#BSUB -R "rusage[mem=4G]"
#BSUB -R "span[block=24]"
#BSUB -W 3:00
##BSUB -env MPIVER
#BSUB -R "select[model == XeonE5_2650v4]"

unset LSB_AFFINITY_HOSTFILE

module load studio
module load "mpi/${MPIVER:-3.1.5-gcc-9.2.0}"

# BASH_LIB=$(pwd)
# ROOT_LIB=$BASH_LIB/..
# SOURCE_LIB=$ROOT_LIB/Source

H=(200)
W=(500)
N=(4)
# H=(500)
# W=(200)
# N=(8)

iter_max=2000
tolerance=0
start_T=15
output_type=0

INSTANCE="mpi-omp"

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
            # export OMP_NUM_THREADS="${N[$n]}"
            echo "OMP_NUM_THREADS=${N[$n]} mpirun --map-by ppr:2:node ./poisson ${H[$h]} ${W[$w]} $iter_max $tolerance $start_T $output_type"
            mpirun -n 4 --map-by node --bind-to socket --rank-by socket -x "OMP_NUM_THREADS=${N[$n]}" ./poisson "${H[$h]}" "${W[$w]}" $iter_max $tolerance $start_T $output_type
            # mpirun -n 4 -x "OMP_NUM_THREADS=${N[$n]}" ./poisson "${H[$h]}" "${W[$w]}" $iter_max $tolerance $start_T $output_type
        done
    done
done
