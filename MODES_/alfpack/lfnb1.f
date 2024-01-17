      SUBROUTINE LFNB1 (INIT,M,L,I,PB,PZ,P1,A,B)
      DIMENSION       PB(1)       ,PZ(1)      ,P1(1)      ,A(1)       ,
     1                B(1)
      SAVE COST
      IF(INIT.NE.0) GO TO 50
      LM1 = L-1
      LM2 = L-2
      PI = 3.14159265358979
      THETA = (I-1)*PI/LM1
      COST = COS(THETA)
      LS = 0
      DO  40 N=LM2,LM1
         LS = LS+1
         CALL ALFK (N,0,A)
         CALL LFPT (N,0,THETA,A,PZ(LS))
         CALL ALFK (N,1,A)
         CALL LFPT (N,1,THETA,A,P1(LS))
   40 CONTINUE
      PHOLD1 = PZ(1)
      PHOLD2 = P1(2)
      PZ(1) = PZ(2)
      P1(2) = P1(1)
      P1(1) = PHOLD1
      PZ(2) = PHOLD2
      CALL RECM (L,LM2,I,P1,A,B)
      CALL RECM (L,LM1,I,PZ,A,B)
      RETURN
   50 K = 0
      LMM = L-M
      IF (LMM-1)  60, 60, 70
   60 PB(L) = PZ(M+1)
      GO TO 160
   70 N = L-1
      FNMM = N-M+1
      FNPN = N+N+1
      FNPM = N+M+1
   80 N = N-1
      IF (N-M) 100, 90, 90
   90 FNMM = FNMM-1.
      FNPN = FNPN-2.
      FNPM = FNPM-1.
      K = K+1
      A(K) = SQRT(FNMM*FNPM/(FNPN*(FNPN+2.)))
      B(K) = -COST
      GO TO  80
  100 IF (ABS(PZ(M+1)) .LT. ABS(P1(M+1))) GO TO 120
      PB(1) = PZ(M+1)
      PB(2) = -A(1)*PB(1)
      IF(LMM-1) 140,140,110
  110 CALL TRIVEC (LMM-1,A,B,A(2),PB(2))
      GO TO 140
  120 PB(1) = PZ(M+1)
      PB(2) = P1(M+1)
      PB(3) = -A(2)*PB(2)
      IF(LMM-2) 140,140,130
  130 CALL TRIVEC (LMM-2,A(2),B(2),A(3),PB(3))
  140 NDO=(L+1)/2
      DO 150 N=1,NDO
         N1 = L-N
         PHOLD = PB(N1+1)
         PB(N1+1) = PB(N)
         PB(N) = PHOLD
  150 CONTINUE
  160 IF(M.EQ.0) RETURN
      DO 170 N=1,M
         PB(N) = 0.
  170 CONTINUE
      RETURN
      END
