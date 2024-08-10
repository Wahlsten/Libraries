program LeftRiemann
      IMPLICIT NONE
      INTEGER, PARAMETER :: RP = SELECTED_REAL_KIND(15)
      INTEGER            :: i, N
      REAL(KIND=RP)      :: a, b, x_i, dx, sum

      a = -2.0
      b = 4.0
      N = 25

      dx = (b - a)/N
      sum = 0.0_RP

      DO i = 0, N-1
        x_i = a + i * dx
        sum = sum + dx * x_i**2 / 3.0_RP
      END DO

      WRITE(*,*)sum

end program leftRiemann
