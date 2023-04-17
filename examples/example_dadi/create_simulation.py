import numpy as np
import subprocess # for running bcftools/plink
import time



'''def unrollVEP(path_to_files: str, path_to_bcf: str):
    """Unrolls the VEP field into separate columns

    Args:
        path_to_files (str): where the filtered vcf/file is 
        path_to_vcf (str): need the path to original vcf file so we can see the header of the vcf to extract VEP field 
        (bcftools view -h gnomad.exomes.r2.1.1.sites.1.vcf.bgz | grep -i 'ID=vep')
        path_to_bcf (str): location of bcftools
    """    

    #location of plink /usr/bin/bcftools
    def execute_command(command):
        process = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        out, errs = process.communicate()
        check_returncode = process.poll()
        if check_returncode == 0:
            print("Found header for vep")
        else:
            raise Exception("Something went wrong: {}".format(errs))
        return out
    
    out = execute_command("{} view -h {} | grep -i \'ID=vep\'".format(path_to_bcf, path_to_files)).decode('utf-8')
    # Get everything after keyword format in field and before ending double quote example: (##INFO..,Description="....Format: ....|Lof_info") -- note the space
    # Regex will get "....|Lof_info
    regex = r"(?<=Format: ).*?(?=\">)"
    match = re.findall(regex, out) 
    assert len(match) == 1, print("There should only be one large VEP field, multiple found {} {}".format(len(match), match))

    # Begin to unroll
    s = match[0].split('|')

    return s
'''

def run_sim(path: str, low: float, high: float, size: int):
    """Runs simulation based on parameters

    Args:
        path (str): path to bash script that runs the simulation
        low (float): low endpoint of linspace
        high (float): high endpoint of linspace
        size (int): how many samples from the linspace
    """

    def execute_command(command):
        print("Executing simulation for {}".format(command))
        process = subprocess.Popen(command.split(), shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
        out, errs = process.communicate()
        print(out)
        check_returncode = process.poll()
        if check_returncode == 0:
            print("Ran simulation")
        else:
            raise Exception("Something went wrong: {}".format(errs))
        return out
    samples = np.linspace(low, high, size)
    for sample in samples:
        start = time.time()
        sample = -1*np.abs(sample)
        out = execute_command("bash {} {}".format(path, sample))
        end = time.time()
        minutes = (end-start)/60
        print("This sample {} took {} minutes to calculate".format(sample, minutes))


    
run_sim('examples/example_dadi/run_single_sel_coefficient.sh', -.990, -1.0, 100)
