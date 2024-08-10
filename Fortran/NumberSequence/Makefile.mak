F90 = C:\Users\Quake\OneDrive\Dokument\Coding\Fortran\NumberSequence
FFLAGS = -Ofast

# Object Files for build
OBJS = n
FibonacciModule.o n

NumberSequence : $(OBJS)
$F90 -o $@ $(OBJS)
# Object dependencies and compilation
FibonacciModule.o : ./FibonacciModule.f90
$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ ./FibonacciModule.f90