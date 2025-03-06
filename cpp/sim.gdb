file sim
break simulation.cpp:95 if t == 647
r --particles 5 --timesteps 651
break workers.cpp:71 if hci == 14 || hci == 13
c
