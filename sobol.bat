if exist a.exe del a.exe
gfortran -Wall -Wextra -fcheck=all -fbacktrace kind.f90 stats.f90 sobol.f90 xsobol.f90
if exist a.exe a.exe
