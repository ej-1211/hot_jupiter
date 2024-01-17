      SUBROUTINE DLFNC1(INIT,M,L,THETA,PB,PZ,P1,A,B,W)
      DOUBLE PRECISION  PB          ,PZ         ,P1     ,
     1                  A           ,B          ,W      ,
     2                  FNMM        ,FNPN       ,FNPM   ,
     3                  PZI         ,COST       ,THETA  ,
     4                  P1I         ,R          ,PHOLD
      DIMENSION         PB(1)       ,PZ(1)      ,P1(1),
     1                  A(1)        ,B(1)       ,W(1)
      SAVE LM1, LM2, LMM
      DO  10 N=1,L
         PB(N) = 0.
   10 CONTINUE
      IF( INIT.NE.0) GO TO 70
      LM1 = L-1
      LM2 = L-2
      LMM = L-M
      PI = 4.d0*DATAN(1.D0)
      CALL DALFK(LM2,M,P1)
      CALL DALFK(LM1,M,PZ)
      K = 0
      N = L-1
      FNMM = N-M+1
      FNPN = N+N+1
      FNPM = N+M+1
   50 N = N-1
      IF(N.LT.M) RETURN
   60 FNMM = FNMM-1.d0
      FNPN = FNPN-2.d0
      FNPM = FNPM-1.d0
      K = K+1
      B(K) = DSQRT(FNMM*FNPM/(FNPN*(FNPN+2.d0)))
      GO TO  50
   70 CALL DLFPT(LM1,M,THETA,PZ,PZI)
      PB(L) = PZI
      IF(L.EQ.1) RETURN
      CALL DLFPT(LM2,M,THETA,P1,P1I)
      IF(LMM-2) 140,80,90
   80 PB(L-1) = P1I
      GO TO 140
   90 COST = DCOS(THETA)
      DO 100 K=1,LMM
         PB(K) = -COST
         A(K) = B(K)
  100 CONTINUE
      IF (ABS(PZI) .LT. ABS(P1I)) GO TO 110
      PB(1) = PZI
      R = -A(1)*PB(1)
      CALL DTRIH (LMM-1,A,PB(2),A(2),R)
      GO TO 120
  110 PB(1) = PZI
      PB(2) = P1I
      IF(LMM.EQ.1) GO TO 120
      R = -A(2)*PB(2)
      CALL DTRIH (LMM-2,A(2),PB(3),A(3),R)
  120 NDO = (L+1)/2
      DO 130 N=1,NDO
         N1 = L-N
         PHOLD = PB(N1+1)
         PB(N1+1) = PB(N)
         PB(N) = PHOLD
  130 CONTINUE
  140 RETURN
      END
