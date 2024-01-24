      SUBROUTINE DRECN (L,M,THETA,B,A)
      DIMENSION       A(1)       ,B(1)
      DOUBLE PRECISION
     1                A          ,B          ,COST       ,R1         ,
     2                R2         ,FNMM       ,FNPN       ,FNPM       ,
     3                THETA      ,BHOLD
      COST = DCOS(THETA)
      LMM = L-M
      R1 = B(1)
      R2 = B(2)
      K = 0
      N = L-1
      FNMM = N-M+1
      FNPN = N+N+1
      FNPM = N+M+1
   10 N = N-1
      IF (N-M)  30, 20, 20
   20 FNMM = FNMM-1.
      FNPN = FNPN-2.
      FNPM = FNPM-1.
      K = K+1
      A(K) = DSQRT(FNMM*FNPM/(FNPN*(FNPN+2.)))
      B(K+1) = -COST
      GO TO  10
   30 IF (DABS(R1) .LT. DABS(R2)) GO TO  40
      B(1) = R1
      R1 = -A(1)*R1
      CALL DTRIH (LMM-1,A,B(2),A(2),R1)
      GO TO  50
   40 B(1) = R1
      B(2) = R2
      R2 = -A(2)*R2
      CALL DTRIH (LMM-2,A(2),B(3),A(3),R2)
   50 NDO = (L+1)/2
      DO  60 N=1,NDO
         N1 = L-N
         BHOLD = B(N1+1)
         B(N1+1) = B(N)
         B(N) = BHOLD
   60 CONTINUE
      RETURN
      END
