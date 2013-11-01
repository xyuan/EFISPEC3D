SUBROUTINE in_poly_curve(X0, Y0, X, Y, N, L, M)
   !-----------------------------------------------------------------------
   ! GIVEN A POLYGONAL LINE CONNECTING THE VERTICES (X(I),Y(I)) (I = 1,...,N) TAKEN IN THIS ORDER.  
   ! IT IS ASSUMED THAT THE POLYGONAL PATH IS A LOOP, WHERE (X(N),Y(N)) = (X(1),Y(1)) OR THERE IS 
   ! AN ARC FROM (X(N),Y(N)) TO (X(1),Y(1)).  
   !
   ! N.B. THE POLYGON MAY CROSS ITSELF ANY NUMBER OF TIMES.
   !
   ! (X0,Y0) IS AN ARBITRARY POINT AND L AND M ARE VARIABLES.
   ! ON OUTPUT, L AND M ARE ASSIGNED THE FOLLOWING VALUES:
   !    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
   !    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
   !    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH
   !    M =  0   IF (X0,Y0) IS ON OR OUTSIDE THE PATH.  IF (X0,Y0) IS INSIDE THE
   !             PATH THEN M IS THE WINDING NUMBER OF THE PATH AROUND THE POINT (X0,Y0).
   !
   ! FORTRAN 66 VERSION BY A.H. MORRIS
   ! CONVERTED TO ELF90 COMPATIBILITY BY Alan MILLER, 15 FEBRUARY 1997
   !
   IMPLICIT NONE
   !
   REAL            , INTENT(IN)  :: X0, Y0, X(N), Y(N)
   INTEGER         , INTENT(IN)  :: N
   INTEGER         , INTENT(OUT) :: L, M
   
   !LOCAL VARIABLES
   INTEGER                       :: I, N0
   REAL                          :: ANGLE
   REAL                          :: EPS
   REAL                          :: PI
   REAL                          :: PI2
   REAL                          :: SUM
   REAL                          :: THETA
   REAL                          :: THETA1
   REAL                          :: THETAI
   REAL                          :: TOL
   REAL                          :: U
   REAL                          :: V
   
   !     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
   !            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0
   !
   EPS = EPSILON(1.0)
   !
   !-----------------------------------------------------------------------
   N0 = N
   IF (X(1) == X(N) .AND. Y(1) == Y(N)) N0 = N - 1
   PI  = ATAN2(0.0, -1.0)
   PI2 = 2.0*PI
   TOL = 4.0*EPS*PI
   L   = -1
   M   =  0
   !
   U = X(1) - X0
   V = Y(1) - Y0
   IF (U == 0.0 .AND. V == 0.0) GO TO 20
   IF (N0 < 2) RETURN
   THETA1 = ATAN2(V, U)
   !
   SUM   = 0.0
   THETA = THETA1
   DO I = 2, N0
     U = X(I) - X0
     V = Y(I) - Y0
     IF (U == 0.0 .AND. V == 0.0) GO TO 20
     THETAI = ATAN2(V, U)
   !  
     ANGLE = ABS(THETAI - THETA)
     IF (ABS(ANGLE - PI) < TOL) GO TO 20
     IF (ANGLE > PI    ) ANGLE = ANGLE - PI2
     IF (THETA > THETAI) ANGLE = -ANGLE
     SUM   = SUM + ANGLE
     THETA = THETAI
   END DO
   !
   ANGLE = ABS(THETA1 - THETA)
   IF (ABS(ANGLE - PI) < TOL) GO TO 20
   IF (ANGLE > PI    ) ANGLE = ANGLE - PI2
   IF (THETA > THETA1) ANGLE = -ANGLE
   SUM = SUM + ANGLE
   !
   !SUM = 2*PI*M WHERE M IS THE WINDING NUMBER
   M = ABS(SUM)/PI2 + 0.2
   IF (M == 0) RETURN
   L = 1
   IF (SUM < 0.0) M = -M
   RETURN
   !
   !(X0, Y0) IS ON THE BOUNDARY OF THE PATH
   20 L = 0

   RETURN
END SUBROUTINE in_poly_curve
