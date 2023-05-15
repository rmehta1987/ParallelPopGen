import numpy as np

# From 0 to .0001
a = np.linspace(0,.0001,200)

# From .0001 to .001
b = np.linspace(0.0001,.001,200)

# From .001 to .01
c = np.linspace(0.001,.01,200)

# From .01 to .1
d = np.linspace(0.01,0.1,100)

e = np.linspace(0.1,0.25,25)

samps = np.array_split(np.unique(np.concatenate((a, b, c,d,e ))),100)

for j, a_samp in enumerate(samps):
    np.savetxt('sample_{}.txt'.format(j),a_samp,fmt='%.6e')

