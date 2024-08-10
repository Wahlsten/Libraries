MODULE InterpolationRoutines
   USE Constants
   IMPLICIT NONE

CONTAINS
   SUBROUTINE Interpolate(C, X, X_new, Y_new, N)
      IMPLICIT NONE
      INTEGER                           ,INTENT(IN)  :: N
      REAL(KIND = RP), DIMENSION(0:N)   ,INTENT(IN)  :: C, X
      REAL(KIND = RP), DIMENSION(0:2*N) ,INTENT(OUT) :: X_new, Y_new
      REAL(KIND = RP), EXTERNAL                      :: Poly_interpolant

      ! Local variables
      INTEGER                                        :: i

      DO i = 0,2*N
         X_new(i) = i * (2.0_RP * pi)/N
         Poly_Interpolant(X_new(i), X, C, B, Y_new(i))
      END DO

      RETURN
   END SUBROUTINE Interpolate

!---------------------------------------------------------

   SUBROUTINE NewtonCoeff(X, Y, C, N)
      IMPLICIT NONE
      INTEGER,                       INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: X, Y
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: C

      ! Local variables
      REAL(KIND=RP)                              :: d, u
      INTEGER                                    :: i, j

      C(0) = Y(0)
      DO j = 1,N
         d = X(j) - X(j-1)
         u = C(j-1)
         DO i = j-2, 0, -1
            u = u * (X(j)-X(i)) + C(i)
            d = d * (X(j) - X(i))
         END DO
         C(j) = (Y(j) - u)/d
      END DO

      RETURN
   END SUBROUTINE NewtonCoeff

!---------------------------------------------------------

   SUBROUTINE Poly_Interpolant(x, nodes, coeff, N, y)
      IMPLICIT NONE
      INTEGER,                       INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: nodes, coeff
      REAL(KIND=RP),                 INTENT(IN)  :: x
      REAL(KIND=RP),                 INTENT(OUT) :: y

      ! Local variables
      INTEGER      :: i,j
      REAL(KIND=RP) :: temp

      y = 0.0_RP
      DO i = 0, N
         temp = 1.0_RP
         DO j = 0, i-1
            temp = temp * (x - nodes(j))
         END DO
         y = y + coeff(i) * temp
      END DO

      RETURN
   END SUBTROUTINE Poly_Interpolant

!--------------------------------------------------------

   SUBROUTINE ReadInArrays(X, Y, N)
      IMPLICIT NONE
      INTEGER,                       INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: X, Y
     
      ! Local variables
      INTEGER :: i

      OPEN(12, FILE='poly.dat')
      DO i = 0, N
         READ(12,*) X(i), Y(i)
      END DO

      RETURN
   END SUBROUTINE

!--------------------------------------------------------

   SUBROUTINE WriteData(X, Y, N)
      IMPLICIT NONE
      INTEGER,                       INTENT(IN) :: N
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: X, Y
      ! Local variables
      INTEGER :: i

      OPEN(12, FILE='output.dat')
         WRITE(12, *) X(i), Y(i)
      END DO
      CLOSE(12)

      RETURN
   END SUBROUTINE WriteData

!--------------------------------------------------------

END MODULE InterpolationRoutines
