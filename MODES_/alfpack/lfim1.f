      SUBROUTINE LFIM1(INIT,THETA,L,N,NM,ID,P3,PHZ,PH1,P1,P2,CP)
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
    5 IF(N .GT. 2) GO TO 60
      IF(N-1)25,30,35
   25 DO 45 I=1,L
      P3(I,1)=PHZ(I,1)
   45 CONTINUE
      RETURN
   30 DO 50 I=1,L
      P3(I,1) = PHZ(I,2)
      P3(I,2) = PH1(I,2)
   50 CONTINUE
      RETURN
   35 SQ5S6 = SQRT(5./6.)
      SQ1S6 = SQRT(1./6.)
      DO 55 I=1,L
      P3(I,1) = PHZ(I,3)
      P3(I,2) = PH1(I,3)
      P3(I,3) = SQ5S6*PHZ(I,1)-SQ1S6*P3(I,1)
      P1(I,1) = PHZ(I,2)
      P1(I,2) = PH1(I,2)
      P2(I,1) = PHZ(I,3)
      P2(I,2) = PH1(I,3)
      P2(I,3) = P3(I,3)
   55 CONTINUE
      RETURN
   60 NM1 = N-1
      NP1 = N+1
      FN = FLOAT(N)
      TN = FN+FN
      CN = (TN+1.)/(TN-3.)
      DO 65 I=1,L
      P3(I,1) = PHZ(I,NP1)
      P3(I,2) = PH1(I,NP1)
   65 CONTINUE
      IF(NM1 .LT. 3) GO TO 71
      DO 70 MP1=3,NM1
      M = MP1-1
      FM = FLOAT(M)
      FNPM = FN+FM
      FNMM = FN-FM
      TEMP = FNPM*(FNPM-1.)
      CC = SQRT(CN*(FNPM-3.)*(FNPM-2.)/TEMP)
      DD = SQRT(CN*FNMM*(FNMM-1.)/TEMP)
      EE = SQRT((FNMM+1.)*(FNMM+2.)/TEMP)
      DO 70 I=1,L
      P3(I,MP1) = CC*P1(I,MP1-2)+DD*P1(I,MP1)-EE*P3(I,MP1-2)
   70 CONTINUE
   71 FNPM = FN+FN-1.
      TEMP = FNPM*(FNPM-1.)
      CC = SQRT(CN*(FNPM-3.)*(FNPM-2.)/TEMP)
      EE = SQRT(6./TEMP)
      DO 75 I=1,L
      P3(I,N) = CC*P1(I,N-2)-EE*P3(I,N-2)
   75 CONTINUE
      FNPM = FN+FN
      TEMP = FNPM*(FNPM-1.)
      CC = SQRT(CN*(FNPM-3.)*(FNPM-2.)/TEMP)
      EE = SQRT(2./TEMP)
      DO 80 I=1,L
      P3(I,N+1) = CC*P1(I,N-1)-EE*P3(I,N-1)
   80 CONTINUE
      DO 90 MP1=1,NP1
      DO 90 I=1,L
      P1(I,MP1) = P2(I,MP1)
      P2(I,MP1) = P3(I,MP1)
   90 CONTINUE
      RETURN
      END
