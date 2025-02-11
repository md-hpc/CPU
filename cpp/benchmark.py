import time
import os

for N in range(3,20):
    print(N, end=' ')
    
    start = time.time()
    os.system(f"./cells --timesteps 10 --universe-size {N} --particles {80 * N ** 3}")
    stop = time.time()
    print(stop - start, end=' ')
    
    start = time.time()
    os.system(f"./list-thread --timesteps 10 --universe-size {N} --particles {80 * N ** 3}")
    stop = time.time()
    print(stop - start)
