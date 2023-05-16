#!/bin/bash

# Get first job id

jid=$(sbatch create_slurm_template_205.sh | cut -d ' ' -f4)

# Remainder jobs
for k in {210..300..5};
    do temp="${k}"
        jid=$(sbatch --dependency=afterok:${jid} create_slurm_template_${k}.sh | cut -d ' ' -f4)
    done
