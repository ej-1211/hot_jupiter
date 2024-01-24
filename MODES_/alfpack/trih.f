      SUBROUTINE TRIH (N,A,B,C,R)
      DIMENSION       A(1)       ,B(1)       ,C(1)
      IF (N) 120,120, 10
   10 IF (N-2)  20, 30, 40
   20 B(1) = R/B(1)
      RETURN
   30 QIH = A(2)
      BIH = B(2)
      GO TO  70
   40 QIH = A(N)
      BIH = B(N)
      DO  60 IDO=3,N
         I = N-IDO+2
         IF (ABS(BIH) .LT. ABS(C(I))) GO TO  50
         RATIO = C(I)/BIH
         C(I) = 0.
         B(I+1) = QIH/BIH
         BIH = B(I)-RATIO*QIH
         QIH = A(I)
         GO TO  60
   50    B(I+1) = B(I)/C(I)
         C(I) = A(I)/C(I)
         BIH1 = QIH-BIH*B(I+1)
         QIH = -BIH*C(I)
         BIH = BIH1
   60 CONTINUE
   70 IF (ABS(BIH) .LT. ABS(C(1))) GO TO  80
      Q2 = QIH/BIH
      BIH = B(1)-C(1)/BIH*QIH
      B(1) = R/BIH
      B(2) = -Q2*B(1)
      GO TO  90
   80 RATIO = BIH/C(1)
      BIH = QIH-RATIO*B(1)
      RIH = -RATIO*R
      B1 = RIH/BIH
      B(2) = (R-B(1)*B1)/C(1)
      B(1) = B1
   90 IF (N-3) 120,100,100
  100 DO 110 I=3,N
         B(I) = -B(I)*B(I-1)-C(I-1)*B(I-2)
  110 CONTINUE
  120 RETURN
      END
