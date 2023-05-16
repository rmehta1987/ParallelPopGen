#!/usr/bin/env python
import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


def make_scripts(sample_location: str, sample_size: int, slurm_prefix: str, script_name: str):
    samples_directory = os.listdir(sample_location)  # Get names of sample directories in the sample folder
    samples_directory.sort()
    #sratch = os.environ['SCRATCH']
    #data_dir = os.path.join(scratch, '/project/LizardLips')

    # Make top level directories
    #mkdir_p(job_directory)
    #mkdir_p(data_dir)

    #lizards=["LizardA","LizardB"]

    for j, a_folder in enumerate(samples_directory):

        folder_num = int(a_folder.split('_')[0]) # hacky for now
        job_file = f'{slurm_prefix}_{folder_num}.sh'
        path_name = f'{sample_location}{a_folder}'

        with open(job_file, 'w') as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%d_samp.job\n" % folder_num)
            fh.writelines("#SBATCH --output=%d.out\n" % folder_num)
            fh.writelines("#SBATCH --time=12:00:00\n")
            fh.writelines("#SBATCH --partition=gpu\n")
            fh.writelines("#SBATCH --gres=gpu:1\n")
            fh.writelines("#SBATCH --account=pi-jjberg\n")
            fh.writelines("#SBATCH --ntasks-per-node=1\n")
            fh.writelines("#SBATCH --cpus-per-task=2\n\n")

            fh.writelines("module load cuda\n")
            fh.writelines("bash {} {} {}\n".format(script_name, path_name, sample_size))

            #os.system("sbatch %s" %job_file)
            #break # for testing


make_scripts('samples/', 170, 'create_slurm_template', 'run_lots_of_sims_template.sh')
