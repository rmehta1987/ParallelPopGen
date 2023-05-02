import numpy as np
import os
import re


regex = r"{(.*?)}"
files = os.listdir()
path='/project/jjberg/mehta5/ParallelPopGen/examples/example_dadi/linspace_samples/'
sample_files = os.listdir('/project/jjberg/mehta5/ParallelPopGen/examples/example_dadi/linspace_samples/')
the_gammas = [np.loadtxt(path+b,dtype=object) for b in sample_files]
the_gammas = list(np.concatenate(the_gammas))
print("length of gamma samples {}.".format(len(the_gammas)))

need = []
for a_file in files:
    matches = re.search(regex, a_file)
    if matches is not None:
        gamma = matches.group()[1:-1] # get everything between curly brace besides the curly brace
        if gamma in the_gammas:
            the_gammas.remove(gamma)
the_gammas = np.array_split(np.asarray(the_gammas,dtype=float),50)
for i, a_file in enumerate(the_gammas):
    np.savetxt('sample_{}.txt'.format(i),a_file,fmt='%.6e')





