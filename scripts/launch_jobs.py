import sys

ks = eval(sys.argv[1])
ss = eval(sys.argv[2])
ms = eval(sys.argv[3])

for k in ks:
    for s in ss:
        for m in ms:
            with open(f"k_{k}_s{s}_m{m}.sbatch", "w") as outfile:
                outfile.write(f"""#!/bin/bash --login                                               
########## SBATCH Lines for Resource Request ##########                         

#SBATCH --time=3:59:00
#SBATCH --job-name m{m}_k{k}_mult{s} 
#SBATCH --mem-per-cpu=500M                                             
####SBATCH --array=1-10                                                            
########## Command Lines to Run ##########                                      

mkdir /mnt/scratch/dolsonem/genomic_stability/data/MUTPROB_{m}_K_{k}_FITNESSMULT_{s}
cd /mnt/scratch/dolsonem/genomic_stability/data/MUTPROB_{m}_K_{k}_FITNESSMULT_{s}

#mkdir $SLURM_ARRAY_TASK_ID
#cd $SLURM_ARRAY_TASK_ID
for i in {{1..10}}
do
    mkdir $i
    cd $i
    ../../../genomic_stability_model -TIME_STEPS 1000000 -INITIAL_FITNESS 1 -GAMMA_K {k} -MUT_PROB {m} -FITNESS_MULT {s} -INIT_POP_SIZE 500 > run.log
    cd ..
done
""")