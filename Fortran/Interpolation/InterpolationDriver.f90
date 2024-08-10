PROGRAM InterpolationDriver
   USE InterpolationRoutines
   IMPLICIT NONE
   INTEGER              :: N = 12
   REAL(KIND=RP), DIMENSION(0:N) :: X, Y, c
   REAL(KIND=RP), DIMENSION(0:2*N) :: new_X, new_Y
   INTEGER                         :: i

   OPEN(13, FILE='poly.dat')
   DO i = 0, N
      X(i) = i * (2.0_RP * pi)/N
      Y(i) = COS(X(i))
      WRITE(13, *) X(i), Y(i)
   END DO
   CLOSE(13)

   CALL ReadInArrays(X, Y, N)
   CALL NewtonCoeff(X, Y, c, N)
   CALL Interpolate(c, X, new_X, new_Y, N)
   CALL WriteData(new_X, new_Y, N)

END PROGRAM InterpolationDriver
