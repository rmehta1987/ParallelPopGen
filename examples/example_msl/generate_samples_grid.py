import numpy as np
import os

def generate(low: float, high: float, num: int, path: str):

    samples_uniform = np.linspace(low, high, num)
    how_many_batches = int(num/100)
    batches = np.split(samples_uniform, how_many_batches)

    if not (os.path.isdir(path)):
        try:
            os.mkdir(path)
        except OSError:
            print("Error in making directory")


    for j, batch in enumerate(batches):
        size_of_batch = batch.shape[0]
        if j % 5 == 0:
            dir_path = f'{path}{j}_set_of_samples/'
            if not(os.path.isdir(dir_path)):
                os.mkdir(dir_path)
        path2 = f'{dir_path}{j}_{size_of_batch}_samples.txt'
        np.savetxt(path2, batch, fmt='%.9e')



generate(-.990, -1e-5, 90000, 'samples/')
