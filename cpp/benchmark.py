import time
import os


for N in range(3,30):
    fp = open("benchmark.out","a")
    print(N, end=' ', file=fp)
    
    ppc_c = 80
    ppc_l = int(80 * 1.2 ** 3)

    argv_c = f"--timesteps 17 --universe-size {N} --particles {ppc_c * N ** 3}"
    argv_l = f"--timesteps 17 --universe-size {N} --particles {ppc_l * N ** 3}"
    start = time.time()
    os.system(f"./cells-thread {argv_c}")
    stop = time.time()
    print(stop - start, end=' ', file=fp)
    
    start = time.time()
    os.system(f"./list-thread {argv_l}")
    stop = time.time()
    print(stop - start, end=' ', file=fp)
   
    print(file=fp)

    fp.close()
