      SUBROUTINE DLFMA1 (INIT,N,L,I,PB,PZ,P1,A,B)
      DOUBLE PRECISION 
     1                PB         ,PZ         ,P1         ,A          ,
     2                B          ,PI         ,DT         ,FNPM       ,
     3                FNMM       ,FMPM       ,THETA      ,COST       ,
     4                SINT       ,R
      DIMENSION       PB(1)       ,PZ(1)      ,P1(1)      ,A(1)       ,
     1                B(1)
      SAVE DT, NM1
      IF(INIT.NE.0) GO TO 60
      IF(N.LE.0) RETURN
      PI = 4.*DATAN(1.D0)
      DT = PI/(L-1)
      CALL DALFK (N,0,PB)
      DO 10 II=1,L
      THETA = (II-1)*DT
      CALL DLFPT (N,0,THETA,PB,PZ(II))
   10 CONTINUE
      CALL DALFK (N,1,PB)
      DO 15 II=1,L
      THETA = (II-1)*DT
      CALL DLFPT (N,1,THETA,PB,P1(II))
   15 CONTINUE
      NM1 = N-1
      NP1 = N+1
      FNPM = N
      FNMM = N+1
      FMPM = 0.
      DO  50 M=1,N
         FNPM = FNPM+1.
         FNMM = FNMM-1.
         FMPM = FMPM-2.
         B(M) = DSQRT(FNMM*FNPM)
   50 CONTINUE
      RETURN
   60 IF(N-1) 61,91,65
   61 PB(1) = 1./DSQRT(2.D0)
      RETURN
   65 THETA = FLOAT(I-1)*DT
      COST = DCOS(THETA)
      SINT = DSIN(THETA)
      FMPM = 0.
      DO  70 M=1,N
         FMPM = FMPM-2.
         A(M) = B(M)*SINT
         PB(M+1) = FMPM*COST
   70 CONTINUE
      IF (DABS(PZ(I))-DABS(P1(I)))  90, 90, 80
   80 R = -A(1)*PZ(I)
      CALL DTRIH (N,A,PB(2),A(2),R)
      PB(1) = PZ(I)
      GO TO 100
   90 R = -A(2)*P1(I)
      CALL DTRIH (NM1,A(2),PB(3),A(3),R)
   91 PB(1) = PZ(I)
      PB(2) = P1(I)
  100 RETURN
      END
