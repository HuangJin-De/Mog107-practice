#/bin/bash
MPI=8
#MPI number should be divisor of NT (time steps in your input data)
for YEAR in 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015
do
mpirun -np $MPI ./go.exe $YEAR
done
