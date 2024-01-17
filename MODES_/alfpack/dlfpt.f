C SUBROUTINE DLFPT (N,M,THETA,CP,PB)
C
C DIMENSION OF
C ARGUMENTS
C                        DOUBLE PRECISION CP((N/2)+1)
C
C PURPOSE                ROUTINE DLFPT USES COEFFICIENTS COMPUTED BY
C                        ROUTINE DALFK TO COMPUTE THE DOUBLE PRECISION
C                        NORMALIZED ASSOCIATED LEGENDRE FUNCTION PBAR(N,
C                        M,THETA) AT COLATITUDE THETA.
C
C USAGE                  CALL DLFPT(N,M,THETA,CP,PB)
C
C ARGUMENTS
C
C ON INPUT               N
C                          NONNEGATIVE INTEGER SPECIFYING THE DEGREE OF
C                          PBAR(N,M,THETA)
C                        M
C                          IS THE ORDER OF PBAR(N,M,THETA). M CAN BE
C                          ANY INTEGER. HOWEVER, PBAR(N,M,THETA) = 0
C                          IF ABS(M) IS GREATER THAN N, AND
C                          PBAR(N,M,THETA) = (-1)**M*PBAR(N,-M,THETA)
C                          FOR NEGATIVE M.
C
C                        THETA
C                          DOUBLE PRECISION COLATITUDE IN RADIANS
C
C                        CP
C                          DOUBLE PRECISION ARRAY OF LENGTH (N/2)+1
C                          CONTAINING COEFFICIENTS COMPUTED BY ROUTINE
C                          DALFK
C
C ON OUTPUT              PB
C                          DOUBLE PRECISION VARIABLE CONTAINING
C                          PBAR(N,M,THETA)
C
C SPECIAL CONDITIONS     CALLS TO ROUTINE DLFPT MUST BE PRECEDED BY AN
C                        APPROPRIATE CALL TO ROUTINE DALFK.
C
C PRECISION              DOUBLE
C
C ALGORITHM              THE TRIGONOMETRIC SERIES FORMULA USED BY
C                        ROUTINE DLFPT TO CALCULATE PBAR(N,M,TH) AT
C                        COLATITUDE TH DEPENDS ON M AND N AS FOLLOWS:
C
C                           1) FOR N EVEN AND M EVEN, THE FORMULA IS
C                              .5*CP(1) PLUS THE SUM FROM K=1 TO K=N/2
C                              OF CP(K)*DCOS(2*K*TH)
C                           2) FOR N EVEN AND M ODD. THE FORMULA IS
C                              THE SUM FROM K=1 TO K=N/2 OF
C                              CP(K)*DSIN(2*K*TH)
C                           3) FOR N ODD AND M EVEN, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*DCOS((2*K-1)*TH)
C                           4) FOR N ODD AND M ODD, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*DSIN((2*K-1)*TH)
C
C TIMING                 TIME PER CALL TO ROUTINE DLFPT IS DEPENDENT
C                        ON THE INPUT PARAMETER N.
C
      SUBROUTINE DLFPT (N,M,THETA,CP,PB)
      DOUBLE PRECISION 
     1                CP         ,PB          ,CDT        ,SDT        ,
     2                CT         ,ST         ,SUM        ,CTH     ,THETA
      DIMENSION       CP(1)
      PB = 0.
      MA = IABS(M)
      IF(MA .GT. N) RETURN
      IF (N)  10, 10, 30
   10 IF (MA)  20, 20, 30
   20 PB= DSQRT(5.D-1)
      GO TO 140
   30 NP1 = N+1
      NMOD = MOD(N,2)
      MMOD = MOD(MA,2)
      IF (NMOD)  40, 40, 90
   40 IF (MMOD)  50, 50, 70
   50 KDO = N/2+1
      CDT = DCOS(THETA+THETA)
      SDT = DSIN(THETA+THETA)
      CT = 1.
      ST = 0.
      SUM = .5*CP(1)
      DO  60 KP1=2,KDO
         CTH = CDT*CT-SDT*ST
         ST = SDT*CT+CDT*ST
         CT = CTH
         SUM = SUM+CP(KP1)*CT
   60 CONTINUE
      PB= SUM
      GO TO 140
   70 KDO = N/2
      CDT = DCOS(THETA+THETA)
      SDT = DSIN(THETA+THETA)
      CT = 1.
      ST = 0.
      SUM = 0.
      DO  80 K=1,KDO
         CTH = CDT*CT-SDT*ST
         ST = SDT*CT+CDT*ST
         CT = CTH
         SUM = SUM+CP(K)*ST
   80 CONTINUE
      PB= SUM
      GO TO 140
   90 KDO = (N+1)/2
      IF (MMOD) 100,100,120
  100 CDT = DCOS(THETA+THETA)
      SDT = DSIN(THETA+THETA)
      CT = DCOS(THETA)
      ST = -DSIN(THETA)
      SUM = 0.
      DO 110 K=1,KDO
         CTH = CDT*CT-SDT*ST
         ST = SDT*CT+CDT*ST
         CT = CTH
         SUM = SUM+CP(K)*CT
  110 CONTINUE
      PB= SUM
      GO TO 140
  120 CDT = DCOS(THETA+THETA)
      SDT = DSIN(THETA+THETA)
      CT = DCOS(THETA)
      ST = -DSIN(THETA)
      SUM = 0.
      DO 130 K=1,KDO
         CTH = CDT*CT-SDT*ST
         ST = SDT*CT+CDT*ST
         CT = CTH
         SUM = SUM+CP(K)*ST
  130 CONTINUE
      PB= SUM
  140 RETURN
      END
