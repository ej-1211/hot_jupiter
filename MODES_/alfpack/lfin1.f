      SUBROUTINE LFIN1(INIT,THETA,L,M,NM,ID,P3,PHZ,PH1,P1,P2,CP)
      DIMENSION       P1(L,1)    ,P2(L,1)    ,P3(ID,1)   ,PHZ(L,1)   ,
     1                PH1(L,1)   ,CP(1)      ,THETA(1)
      NMP1 = NM+1
      IF(INIT .NE. 0) GO TO 5
      SSQRT2 = 1./SQRT(2.)
      DO 10 I=1,L
      PHZ(I,1) = SSQRT2
   10 CONTINUE
      DO 15 NP1=2,NMP1
      NH = NP1-1
      CALL ALFK(NH,0,CP)
      DO 16 I=1,L
      CALL LFPT(NH,0,THETA(I),CP,PHZ(I,NP1))
   16 CONTINUE
      CALL ALFK(NH,1,CP)
      DO 17 I=1,L
      CALL LFPT(NH,1,THETA(I),CP,PH1(I,NP1))
   17 CONTINUE
   15 CONTINUE
      RETURN
    5 MP1 = M+1
      FM = FLOAT(M)
      TM = FM+FM
      IF(M-1)25,30,35
   25 DO 45 NP1=1,NMP1
      DO 45 I=1,L
      P3(I,NP1) = PHZ(I,NP1)
      P1(I,NP1) = PHZ(I,NP1)
   45 CONTINUE
      RETURN
   30 DO 50 NP1=2,NMP1
      DO 50 I=1,L
      P3(I,NP1) = PH1(I,NP1)
      P2(I,NP1) = PH1(I,NP1)
   50 CONTINUE
      RETURN
   35 TEMP = TM*(TM-1.)
      CC = SQRT((TM+1.)*(TM-2.)/TEMP)
      EE = SQRT(2./TEMP)
      DO 85 I=1,L
      P3(I,M+1) = CC*P1(I,M-1)-EE*P1(I,M+1)
   85 CONTINUE
      IF(M .EQ. NM) RETURN
      TEMP = TM*(TM+1.)
      CC = SQRT((TM+3.)*(TM-2.)/TEMP)
      EE = SQRT(6./TEMP)
      DO 70 I=1,L
      P3(I,M+2) = CC*P1(I,M)-EE*P1(I,M+2)
   70 CONTINUE
      MP3 = M+3
      IF(NMP1 .LT. MP3) GO TO 80
      DO 75 NP1=MP3,NMP1
      N = NP1-1
      FN = FLOAT(N)
      TN = FN+FN
      CN = (TN+1.)/(TN-3.)
      FNPM = FN+FM
      FNMM = FN-FM
      TEMP = FNPM*(FNPM-1.)
      CC = SQRT(CN*(FNPM-3.)*(FNPM-2.)/TEMP)
      DD = SQRT(CN*FNMM*(FNMM-1.)/TEMP)
      EE = SQRT((FNMM+1.)*(FNMM+2.)/TEMP)
      DO 75 I=1,L
      P3(I,NP1) = CC*P1(I,NP1-2)+DD*P3(I,NP1-2)-EE*P1(I,NP1)
   75 CONTINUE
   80 DO 90 NP1=M,NMP1
      DO 90 I=1,L
      P1(I,NP1) = P2(I,NP1)
      P2(I,NP1) = P3(I,NP1)
   90 CONTINUE
      RETURN
      END
