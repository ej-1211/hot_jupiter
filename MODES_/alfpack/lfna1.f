      SUBROUTINE LFNA1 (INIT,M,L,I,PB,PZ,P1,A,B,W)
      DIMENSION       PB(1)       ,PZ(1)      ,P1(1)      ,A(1)       ,
     1                B(1)        ,W(1)
      SAVE LM1, LMM, PI
      IF(INIT.NE.0) GO TO 70
      LM1 = L-1
      LM2 = L-2
      LMM = L-M
      PI = 3.14159265358979
      CALL ALFK (LM2,M,B)
      CALL LFP (0,LM2,M,L,B,PB,W)
      CALL LFP (1,LM2,M,L,B,PB,W)
      CALL ALFK (LM1,M,B)
      CALL LFP (0,LM1,M,L,B,PZ,W)
      CALL LFP (1,LM1,M,L,B,PZ,W)
      DO 40 N=1,L
         P1(N) = PB(N)
   40 CONTINUE
      K = 0
      N = L-1
      FNMM = N-M+1
      FNPN = N+N+1
      FNPM = N+M+1
   50 N = N-1
      IF(N.LT.M) RETURN
   60 FNMM = FNMM-1.
      FNPN = FNPN-2.
      FNPM = FNPM-1.
      K = K+1
      B(K) = SQRT(FNMM*FNPM/(FNPN*(FNPN+2.)))
      GO TO  50
   70 PB(L) = PZ(I)
      IF (LMM-2) 140, 80, 90
   80 PB(L-1) = P1(I)
      GO TO 140
   90 THETA = (I-1)*PI/LM1
      COST = COS(THETA)
      DO 100 K=1,LMM
         PB(K) = -COST
         A(K) = B(K)
  100 CONTINUE
      IF (ABS(PZ(I)) .LT. ABS(P1(I))) GO TO 110
      PB(1) = PZ(I)
      R = -A(1)*PB(1)
      CALL TRIH (LMM-1,A,PB(2),A(2),R)
      GO TO 120
  110 PB(1) = PZ(I)
      PB(2) = P1(I)
      R = -A(2)*PB(2)
      CALL TRIH (LMM-2,A(2),PB(3),A(3),R)
  120 NDO = (L+1)/2
      DO 130 N=1,NDO
         N1 = L-N
         PHOLD = PB(N1+1)
         PB(N1+1) = PB(N)
         PB(N) = PHOLD
  130 CONTINUE
  140 IF(M.EQ.0) RETURN
      DO 150 N=1,M
         PB(N) = 0.
  150 CONTINUE
      RETURN
      END
