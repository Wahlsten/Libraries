! ======================================================
! Main program
! ======================================================
MODULE FareyModule
   IMPLICIT NONE

CONTAINS
! ======================================================
! Function computing the Farey number sequence
! ======================================================
	FUNCTION Farey(n) RESULT(farey_v)
	   IMPLICIT NONE
	   INTEGER, INTENT(IN)     :: N
	   REAL, DIMENSION(3*N)    :: farey_v
	   INTEGER				   :: a, b, c, d, k
	   INTEGER				   :: a_temp, b_temp, l

	   a = 0
	   b = 1
	   c = 1
	   d = N
	   
	   farey_v(1) = 0
	   
	   l = 1
	   DO WHILE (c <= N)
		  l = l + 1
		  k = (N + b) / d
		  a_temp = a
		  b_temp = b
		  a = c
		  b = d
		  c = (k*c-a_temp)
		  d	= (k*d-b_temp)
		  farey_v(l) = REAL(a) / REAL(b)

	   END DO
	   
	END FUNCTION Farey
END MODULE FareyModule
