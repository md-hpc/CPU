import time
import os

bins = [
    "cells",
    "lists",
]

for t in range(1,16):
    fp = open("benchmark.out","a")
    print(t, end=' ', file=fp)
    
    argv = f"--timesteps 17 --threads {t} --universe-size 20"
    for f in bins:
        start = time.time()
        os.system(f"./{f} {argv}")
        stop = time.time()
        print(stop - start, end=' ', file=fp)
    
    print(file=fp)

    fp.close()
