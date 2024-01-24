C SUBROUTINE DALFK (N,M,CP)
C
C DIMENSION OF           DOUBLE PRECISION CP(N/2 + 1)
C ARGUMENTS
C
C PURPOSE                ROUTINE DALFK COMPUTES DOUBLE PRECISION FOURIER
C                        COEFFICIENTS IN THE TRIGONOMETRIC SERIES
C                        REPRESENTATION OF THE NORMALIZED ASSOCIATED
C                        LEGENDRE FUNCTION PBAR(N,M,THETA) FOR USE BY
C                        ROUTINE DLFPT IN CALCULATING DOUBLE
C                        PRECISION PBAR(N,M,THETA).
C
C                        FIRST DEFINE THE NORMALIZED ASSOCIATED
C                        LEGENDRE FUNCTIONS
C
C                        PBAR(M,N,THETA) = SQRT((2*N+1)*FACTORIAL(N-M)
C                        /(2*FACTORIAL(N+M)))*DSIN(THETA)**M/(2**N*
C                        FACTORIAL(N)) TIMES THE (N+M)TH DERIVATIVE OF
C                        (X**2-1)**N WITH RESPECT TO X=DCOS(THETA)
C
C                        WHERE THETA IS COLATITUDE.
C
C                        THEN SUBROUTINE DALFK COMPUTES THE COEFFICIENTS
C                        CP(K) IN THE FOLLOWING TRIGONOMETRIC
C                        EXPANSION OF PBAR(M,N,THETA).
C
C                        1) FOR N EVEN AND M EVEN, PBAR(M,N,THETA) =
C                           .5*CP(1) PLUS THE SUM FROM K=1 TO K=N/2
C                           OF CP(K)*DCOS(2*K*TH)
C
C                        2) FOR N EVEN AND M ODD, PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=N/2 OF
C                           CP(K)*DSIN(2*K*TH)
C
C                        3) FOR N ODD AND M EVEN, PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=(N+1)/2 OF
C                           CP(K)*DCOS((2*K-1)*TH)
C
C                        4) FOR N ODD AND M ODD,  PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=(N+1)/2 OF
C                           CP(K)*DSIN((2*K-1)*TH)
C
C
C USAGE                  CALL DALFK(N,M,CP)
C
C ARGUMENTS
C
C ON INPUT               N
C                          NONNEGATIVE INTEGER SPECIFYING THE DEGREE OF
C                          PBAR(N,M,THETA)
C
C                        M
C                          IS THE ORDER OF PBAR(N,M,THETA). M CAN BE
C                          ANY INTEGER HOWEVER CP IS COMPUTED SUCH THAT
C                          PBAR(N,M,THETA) = 0 IF ABS(M) IS GREATER
C                          THAN N AND PBAR(N,M,THETA) = (-1)**M*
C                          PBAR(N,-M,THETA) FOR NEGATIVE M.
C
C ON OUTPUT              CP
C                          DOUBLE PRECISION ARRAY OF LENGTH (N/2)+1
C                          WHICH CONTAINS THE FOURIER COEFFICIENTS IN
C                          THE TRIGONOMETRIC SERIES REPRESENTATION OF
C                          PBAR(N,M,THETA)
C
C
C SPECIAL CONDITIONS     NONE
C
C PRECISION              DOUBLE
C
C ALGORITHM              THE HIGHEST ORDER COEFFICIENT IS DETERMINED IN
C                        CLOSED FORM AND THE REMAINIG COEFFICIENTS ARE
C                        DETERMINED AS THE SOLUTION OF A BACKWARD
C                        RECURRENCE RELATION.
C
      SUBROUTINE DALFK (N,M,CP)
      DOUBLE PRECISION 
     1 CP,FNUM,FDEN,FNMH,A1,B1,C1,CP2,FNNP1,FNMSQ,FK,T1,T2,PM1
      DIMENSION       CP(1)
      CP(1) = 0.
      MA = IABS(M)
      IF(MA .GT. N) RETURN
      IF(N-1) 2,3,5
    2 CP(1) = DSQRT(2.D0)
      RETURN
    3 IF(MA .NE. 0) GO TO 4
      CP(1) = DSQRT(1.5D0)
      RETURN
    4 CP(1) = DSQRT(.75D0)
      IF(M .EQ. -1) CP(1) = -CP(1)
      RETURN
    5 IF(MOD(N+MA,2) .NE. 0) GO TO 10
      NMMS2 = (N-MA)/2
      FNUM = N+MA+1
      FNMH = N-MA+1
      PM1 = 1.D0
      GO TO 15
   10 NMMS2 = (N-MA-1)/2
      FNUM = N+MA+2
      FNMH = N-MA+2
      PM1 = -1.D0
   15 T1 = 1.
      T2 = 1.
      IF(NMMS2 .LT. 1) GO TO 20
      FDEN = 2.D0
      DO 18 I=1,NMMS2
      T1 = FNUM*T1/FDEN
      FNUM = FNUM+2.
      FDEN = FDEN+2.
   18 CONTINUE
   20 IF(MA .EQ. 0) GO TO 26
      DO 25 I=1,MA
      T2 = FNMH*T2/(FNMH+PM1)
      FNMH = FNMH+2.
   25 CONTINUE
   26 IF(MOD(MA/2,2) .NE. 0) T1 = -T1
      CP2 = T1*DSQRT((N+.5D0)*T2)/(2.D0**(N-1))
      FNNP1 = N*(N+1)
      FNMSQ = FNNP1-2.D0*MA*MA
      L = (N+1)/2
      IF(MOD(N,2) .EQ. 0 .AND. MOD(MA,2) .EQ. 0) L = L+1
      CP(L) = CP2
      IF(M .GE. 0) GO TO 29
      IF(MOD(MA,2) .NE. 0) CP(L) = -CP(L)
   29 IF(L .LE. 1) RETURN
      FK = N
      A1 = (FK-2.)*(FK-1.)-FNNP1
      B1 = 2.*(FK*FK-FNMSQ)
      CP(L-1) = B1*CP(L)/A1
   30 L = L-1
      IF(L .LE. 1) RETURN
      FK = FK-2.
      A1 = (FK-2.)*(FK-1.)-FNNP1
      B1 = -2.*(FK*FK-FNMSQ)
      C1 = (FK+1.)*(FK+2.)-FNNP1
      CP(L-1) = -(B1*CP(L)+C1*CP(L+1))/A1
      GO TO 30
      END
