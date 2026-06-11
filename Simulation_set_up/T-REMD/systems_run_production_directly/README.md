TPR files contain all essential runtime parameters, topology and initial structures. Navigate to the directory of each system, and you can directly run Replica Exchange Molecular Dynamics (REMD).

Command: mpirun --oversubscribe -np 41 gmx_mpi mdrun -multidir replica_{1..41} -replex 1000 -deffnm md
