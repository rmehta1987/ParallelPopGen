import numpy as np
from itertools import islice

import re
import os



def find_not_sim_samples(data_path_1: str, data_path_2: str, file_name: str):

        # First filter out files to match specific file; ie sfs...{}.txt
    lsdirs = os.listdir(data_path_1)
    regex_remove_last = re.compile(r'((?!last).)*$')
    sel_list_files = list(filter(regex_remove_last.match,lsdirs))

    regex_capture_selection_files = re.compile(r'^.*[{|}].*$')
    sel_list_files = list(filter(regex_capture_selection_files.match, sel_list_files))
    regex_capture_selection_coef = re.compile(r'(?<={)(.*)(?=})')
    sel_coef = []
    for a_file in sel_list_files:
        sel_coef.append(regex_capture_selection_coef.search(a_file)[1]) # gets the actual selection coefficient)
    #print(sel_coef)

    sel_coef = np.asarray(sel_coef,dtype=float)
    print(sel_coef)
    lsdirs = os.listdir(data_path_2)
    samples = []
    for a_file in lsdirs:
        if 'sample_' in a_file:
            print(a_file)
            samples.append(np.loadtxt(f'{data_path_2}/{a_file}', dtype=float))

    all_samps = np.concatenate(samples)
    #print(all_samps)
    #print(all_samps.shape)
    diff = np.setdiff1d(all_samps, sel_coef)
    print("Found this number of un-simulated samples: {}".format(diff.shape))
    #print(diff)
    diffs = np.array_split(diff, diff.shape[0]%40+1)
    for j, a in enumerate(diffs):
        np.savetxt(f'sample_{j}.txt',a,fmt='%.6e')
    #np.savetxt(file_name, diff, fmt='%.6e')
    print('Finished finding not simulated samples')

def main():

    data_path = '/project/jjberg/mehta5/ParallelPopGen/examples/example_robertson/chr10_genome_wide_100k_sample_size'
    data_path2 = '/project/jjberg/mehta5/ParallelPopGen/examples/example_robertson/robertson_samples/'
    file_name = 'needtoknow2.txt'
    find_not_sim_samples(data_path, data_path2, file_name)


if __name__ == "__main__":
    main()


