! ======================================================
! Main program
! ======================================================
MODULE FibonacciModule
   IMPLICIT NONE

CONTAINS
! ======================================================
! Function computing the Fibonacci number sequence
! ======================================================
	FUNCTION Fibonacci(n) RESULT(fibonacci_v)
	   IMPLICIT NONE
	   INTEGER, INTENT(IN)   :: N
	   INTEGER, DIMENSION(N) :: fibonacci_v
	   INTEGER 				 :: i

	   fibonacci_v(1) = 1
	   fibonacci_v(2) = 1
	   
	   DO i = 3, N
	   
		  fibonacci_v(i) = fibonacci_v(i-1) + fibonacci_v(i-2)

	   END DO
	   
	END FUNCTION Fibonacci
END MODULE FibonacciModule
