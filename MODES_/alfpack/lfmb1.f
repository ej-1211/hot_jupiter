      SUBROUTINE LFMB1 (INIT,L,N,I,PB,PMZ,DPMZ,P1,PM1,DPM1,P2)
      DOUBLE PRECISION
     1                DPMZ       ,DPM1       ,DPMZ1      ,DPMZ2      ,
     2                DPM11      ,DPM12      ,DTH        ,DPI
      DIMENSION       P1(1)      ,P2(1)      ,PB(1)       ,PMZ(1)     ,
     1                PM1(1)     ,DPMZ(1)    ,DPM1(1)
      SAVE COST, SINT
      IF(INIT.NE.0) GO TO 50
      LM1 = L-1
      LM2 = L-2
      DPI = 4.*DATAN(1.D0)
      DTH = (I-1)*DPI/LM1
      THETA = DTH
      COST = COS(THETA)
      SINT = SIN(THETA)
      CALL DALFK (LM1,0,P1)
      CALL DLFPT (LM1,0,DTH,P1,DPMZ1)
      CALL DALFK (LM1,1,P1)
      CALL DLFPT (LM1,1,DTH,P1,DPMZ2)
      CALL DALFK (LM2,0,P1)
      CALL DLFPT (LM2,0,DTH,P1,DPM11)
      CALL DALFK (LM2,1,P1)
      CALL DLFPT (LM2,1,DTH,P1,DPM12)
      DPMZ(1) = DPMZ1
      DPMZ(2) = DPM11
      CALL DRECN (L,0,DTH,DPMZ,PM1)
      DO  30 NH=1,L
         PB(NH) = DPMZ(NH)
   30 CONTINUE
      DPM1(1) = DPMZ2
      DPM1(2) = DPM12
      CALL DRECN (L,1,DTH,DPM1,PMZ)
      DO  40 NH=1,L
         PMZ(NH) = PB(NH)
         PM1(NH) = DPM1(NH)
   40 CONTINUE
      RETURN
   50 IF (N) 130, 60, 70
   60 PB(1) = PMZ(1)
      GO TO 130
   70 IF (N-1) 130, 80, 90
   80 PB(1) = PMZ(2)
      PB(2) = PM1(2)
      GO TO 130
   90 FNPM = N
      FNMM = N+1
      FMPM = 0.
      DO 100 M=1,N
         FNPM = FNPM+1.
         FNMM = FNMM-1.
         FMPM = FMPM-2.
         P1(M) = SINT*SQRT(FNMM*FNPM)
         P2(M) = FMPM*COST
  100 CONTINUE
      IF (ABS(PMZ(N+1))-ABS(PM1(N+1))) 120,120,110
  110 PB(2) = -P1(1)*PMZ(N+1)
      CALL TRIVEC (N,P1,P2,P1(2),PB(2))
      PB(1) = PMZ(N+1)
      GO TO 130
  120 PB(3) = -P1(2)*PM1(N+1)
      CALL TRIVEC (N-1,P1(2),P2(2),P1(3),PB(3))
      PB(1) = PMZ(N+1)
      PB(2) = PM1(N+1)
  130 RETURN
      END
