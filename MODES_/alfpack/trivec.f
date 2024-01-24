      SUBROUTINE TRIVEC (N,A,B,C,E)
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,E(1)
      IF (N) 120,120, 10
   10 R = E(1)
      IF (N-2)  20, 30, 40
   20 E(1) = R/B(1)
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
C
C     PB(I+1)=0.
C     Q(I+1)=QIH/BIH
C
         BIH = B(I)-RATIO*QIH
         QIH = A(I)
         GO TO  60
   50    B(I+1) = B(I)/C(I)
         C(I) = A(I)/C(I)
         BIH1 = QIH-BIH*B(I+1)
         QIH = -BIH*C(I)
         BIH = BIH1
C
C     PB(I+1)=A(I)/C(I)
C     Q(I+1)=B(I)/C(I)
C
   60 CONTINUE
   70 IF (ABS(BIH) .LT. ABS(C(1))) GO TO  80
      Q2 = QIH/BIH
C
C     Q(2)=QIH/BIH
C
      BIH = B(1)-C(1)/BIH*QIH
      E(1) = R/BIH
      E(2) = -Q2*E(1)
C
C     E(2)=-Q(2)*E(1)
C
      GO TO  90
   80 RATIO = BIH/C(1)
      BIH = QIH-RATIO*B(1)
      RIH = -RATIO*R
      E(1) = RIH/BIH
      E(2) = (R-B(1)*E(1))/C(1)
   90 IF (N-3) 120,100,100
  100 DO 110 I=3,N
         E(I) = -B(I)*E(I-1)-C(I-1)*E(I-2)
C
C     E(I)=-Q(I)*E(I-1)-PB(I)*E(I-2)
C
  110 CONTINUE
  120 RETURN
      END
