      SUBROUTINE RECM (L,N,I,PB,A,B)
      DIMENSION       PB(1)       ,A(1)       ,B(1)
      PI = 3.14159265358979
      DT = PI/(L-1)
      FNPM = N
      FNMM = N+1
      FMPM = 0.
      THETA = (I-1)*DT
      COST = COS(THETA)
      SINT = SIN(THETA)
      DO  10 M=1,N
         FNPM = FNPM+1.
         FNMM = FNMM-1.
         FMPM = FMPM-2.
         A(M) = SINT*SQRT(FNMM*FNPM)
         B(M) = FMPM*COST
   10 CONTINUE
      IF (ABS(PB(1))-ABS(PB(2)))  30, 30, 20
   20 PB(2) = -A(1)*PB(1)
      CALL TRIVEC (N,A,B,A(2),PB(2))
      GO TO  40
   30 PB(3) = -A(2)*PB(2)
      NM1 = N-1
      CALL TRIVEC (NM1,A(2),B(2),A(3),PB(3))
   40 RETURN
      END
