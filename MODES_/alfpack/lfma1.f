      SUBROUTINE LFMA1 (INIT,N,L,I,PB,PZ,P1,A,B,W)
      DIMENSION       PB(1)       ,PZ(1)      ,P1(1)      ,A(1)       ,
     1                B(1)        ,W(1)
      SAVE DT, NM1
      IF(INIT.NE.0) GO TO 60
      IF(N.LE.0) RETURN
      CALL ALFK (N,0,PB)
      CALL LFP (0,N,0,L,PB,PZ,W)
      CALL LFP (1,N,0,L,PB,PZ,W)
      CALL ALFK (N,1,PB)
      CALL LFP (0,N,1,L,PB,P1,W)
      CALL LFP (1,N,1,L,PB,P1,W)
      NM1 = N-1
      NP1 = N+1
      PI = 3.14159265358979
      DT = PI/(L-1)
      FNPM = N
      FNMM = N+1
      FMPM = 0.
      DO  50 M=1,N
         FNPM = FNPM+1.
         FNMM = FNMM-1.
         FMPM = FMPM-2.
         B(M) = SQRT(FNMM*FNPM)
   50 CONTINUE
      RETURN
   60 IF(N-1) 61,91,65
   61 PB(1) = 1./SQRT(2.)
      RETURN
   65 THETA = FLOAT(I-1)*DT
      COST = COS(THETA)
      SINT = SIN(THETA)
      FMPM = 0.
      DO  70 M=1,N
         FMPM = FMPM-2.
         A(M) = B(M)*SINT
         PB(M+1) = FMPM*COST
   70 CONTINUE
      IF (ABS(PZ(I))-ABS(P1(I)))  90, 90, 80
   80 R = -A(1)*PZ(I)
      CALL TRIH (N,A,PB(2),A(2),R)
      PB(1) = PZ(I)
      GO TO 100
   90 R = -A(2)*P1(I)
      CALL TRIH (NM1,A(2),PB(3),A(3),R)
   91 PB(1) = PZ(I)
      PB(2) = P1(I)
  100 RETURN
      END
