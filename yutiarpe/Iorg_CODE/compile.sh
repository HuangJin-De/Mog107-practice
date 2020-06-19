#/bin/bash
mpifort -free Iorg_idealized.F -heap-arrays 1024 -shared-intel -mcmodel=large -I/home/C.peter50504/.local/include -L//home/C.peter50504/.local/lib -lnetcdff -lnetcdf -o go.exe

