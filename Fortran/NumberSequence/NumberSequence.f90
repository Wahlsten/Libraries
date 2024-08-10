! ======================================================
! Main program
! ======================================================
PROGRAM NumberSequence
   USE FibonacciModule
   USE FareyModule

   IMPLICIT NONE
   INTEGER       		  :: N, a
   INTEGER, DIMENSION(8)  :: number_sequence
   REAL, DIMENSION(24) 	  :: number_sequence2
   
   N = 8
   
   number_sequence = Fibonacci(N)
   number_sequence2 = Farey(N)
   
   Print *, number_sequence
   Print *, number_sequence2

END PROGRAM NumberSequence