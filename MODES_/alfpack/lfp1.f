      SUBROUTINE LFP1(INIT,N,M,L,CP,P,WSAVE1,WSAVE2,WSAVE3)
      DIMENSION CP(1),P(1),WSAVE1(1),WSAVE2(1),WSAVE3(1)
      SAVE LC, LQ, LS
      IF(INIT.NE.0) GO TO 41
      LC=(L+1)/2
      LS=LC-2
      LQ=LC-1
      CALL SINTI(LS,WSAVE1)
      CALL COSTI(LC,WSAVE2)
      CALL COSQI(LQ,WSAVE3)
      RETURN
   41 IF (N)  10, 10, 40
   10 IF (M)  20, 20, 40
   20 SSQRT2 = 1./SQRT(2.)
      DO  30 I=1,L
         P(I) = SSQRT2
   30 CONTINUE
      RETURN
   40 LS2 = (L+1)/2
      LM1 = L-1
      NP1 = N+1
      PI = 4.*ATAN(1.)
      DT = PI/LM1
      NMOD = MOD(N,2)
      MMOD = MOD(M,2)
      IF (NMOD)  50, 50,120
   50 IF (MMOD)  60, 60, 90
   60 KDP = N/2+1
      DO 70 I=1,KDP
      P(I)=.5*CP(I)
   70 CONTINUE
      P(LC)=P(LC)+P(LC)
      CALL COST(LC,P,WSAVE2)
      DO 80 I=1,LC
      LMI=L-I
      P(LMI+1)=P(I)
   80 CONTINUE
      GO TO 190
   90 KDP=N/2
      DO 100 I=1,KDP
      P(I+1)=.5*CP(I)
  100 CONTINUE
      P(LS+2)=0.
      CALL SINT(LS,P(2),WSAVE1)
      DO 110 I=1,LS
      LMI=L-I
      P(LMI)=-P(I+1)
  110 CONTINUE
      P(L)=0.
      GO TO 190
  120 KDP=(N+1)/2
      IF(MMOD)140,140,160
  140 DO 130 I=1,KDP
      P(I)=.25*CP(I)
  130 CONTINUE
      CALL COSQB(LQ,P,WSAVE3)
      DO 150 I=1,LQ
      LMI=L-I
      P(LMI+1)=-P(I)
  150 CONTINUE
      GO TO 190
  160 DO 180 I=1,KDP
      P(I+1)=.25*CP(I)
  180 CONTINUE
      CALL SINQB(LQ,P(2),WSAVE3)
      DO 170 I=1,LQ
      LMI=L-I
      P(LMI)=P(I+1)
  170 CONTINUE
      P(L)=0.
  190 RETURN
      END
