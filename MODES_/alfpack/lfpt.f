C SUBROUTINE LFPT (N,M,THETA,CP,PB)
C
C DIMENSION OF
C ARGUMENTS
C                        CP((N/2)+1)
C
C PURPOSE                ROUTINE LFPT USES COEFFICIENTS COMPUTED BY
C                        ROUTINE ALFK TO COMPUTE THE SINGLE PRECISION
C                        NORMALIZED ASSOCIATED LEGENDRE FUNCTION PBAR(N,
C                        M,THETA) AT COLATITUDE THETA.
C
C USAGE                  CALL LFPT(N,M,THETA,CP,PB)
C
C ARGUMENTS
C
C ON INPUT               N
C                          NONNEGATIVE INTEGER SPECIFYING THE DEGREE OF
C                          PBAR(N,M,THETA)
C                        M
C                          IS THE ORDER OF PBAR(N,M,THETA). M CAN BE
C                          ANY INTEGER HOWEVER PBAR(N,M,THETA) = 0
C                          IF ABS(M) IS GREATER THAN N AND
C                          PBAR(N,M,THETA) = (-1)**M*PBAR(N,-M,THETA)
C                          FOR NEGATIVE M.
C
C                        THETA
C                          SINGLE PRECISION COLATITUDE IN RADIANS
C
C                        CP
C                          SINGLE PRECISION ARRAY OF LENGTH (N/2)+1
C                          CONTAINING COEFFICIENTS COMPUTED BY ROUTINE
C                          ALFK
C
C ON OUTPUT              PB
C                          SINGLE PRECISION VARIABLE CONTAINING
C                          PBAR(N,M,THETA)
C
C SPECIAL CONDITIONS     CALLS TO ROUTINE LFPT MUST BE PRECEDED BY AN
C                        APPROPRIATE CALL TO ROUTINE ALFK.
C
C PRECISION              SINGLE
C
C ALGORITHM              THE TRIGONOMETRIC SERIES FORMULA USED BY
C                        ROUTINE LFPT TO CALCULATE PBAR(N,M,TH) AT
C                        COLATITUDE TH DEPENDS ON M AND N AS FOLLOWS:
C
C                           1) FOR N EVEN AND M EVEN, THE FORMULA IS
C                              .5*CP(1) PLUS THE SUM FROM K=1 TO K=N/2
C                              OF CP(K)*COS(2*K*TH)
C                           2) FOR N EVEN AND M ODD. THE FORMULA IS
C                              THE SUM FROM K=1 TO K=N/2 OF
C                              CP(K)*SIN(2*K*TH)
C                           3) FOR N ODD AND M EVEN, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*COS((2*K-1)*TH)
C                           4) FOR N ODD AND M ODD, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*SIN((2*K-1)*TH)
C
C ACCURACY               COMPARISON BETWEEN ROUTINES LFPT AND DOUBLE
C                        PRECISION DLFPT ON THE CRAY1 INDICATES GREATER
C                        ACCURACY FOR GREATER VALUES ON INPUT PARAMETER
C                        N.  AGREEMENT TO 13 PLACES WAS OBTAINED FOR
C                        N=10 AND TO 12 PLACES FOR N=100.
C
C TIMING                 TIME PER CALL TO ROUTINE LFPT IS DEPENDENT ON
C                        THE INPUT PARAMETER N.
C
      SUBROUTINE LFPT (N,M,THETA,CP,PB)
      DIMENSION       CP(1)
C
      PB = 0.
      MA = IABS(M)
      IF(MA .GT. N) RETURN
      IF (N)  10, 10, 30
   10 IF (MA)  20, 20, 30
   20 PB= SQRT(.5)
      GO TO 140
   30 NP1 = N+1
      NMOD = MOD(N,2)
      MMOD = MOD(MA,2)
      IF (NMOD)  40, 40, 90
   40 IF (MMOD)  50, 50, 70
   50 KDO = N/2+1
      CDT = COS(THETA+THETA)
      SDT = SIN(THETA+THETA)
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
      CDT = COS(THETA+THETA)
      SDT = SIN(THETA+THETA)
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
  100 CDT = COS(THETA+THETA)
      SDT = SIN(THETA+THETA)
      CT = COS(THETA)
      ST = -SIN(THETA)
      SUM = 0.
      DO 110 K=1,KDO
         CTH = CDT*CT-SDT*ST
         ST = SDT*CT+CDT*ST
         CT = CTH
         SUM = SUM+CP(K)*CT
  110 CONTINUE
      PB= SUM
      GO TO 140
  120 CDT = COS(THETA+THETA)
      SDT = SIN(THETA+THETA)
      CT = COS(THETA)
      ST = -SIN(THETA)
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
