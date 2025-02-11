import time
import os

for N in range(3,20):
    fp = open("benchmark.out","a")
    print(N, end=' ', file=fp)
    
    argv = f"--timesteps 10 --universe-size {N} --particles {80 * N ** 3}"

    start = time.time()
    os.system(f"./cells {argv}")
    stop = time.time()
    print(stop - start, end=' ', file=fp)
    
    start = time.time()
    os.system(f"./list {argv}")
    stop = time.time()
    print(stop - start, end=' ', file=fp)
    
    start = time.time()
    os.system(f"./list-thread {argv}")
    stop = time.time()
    print(stop - start, file=fp)
    
    fp.close()
