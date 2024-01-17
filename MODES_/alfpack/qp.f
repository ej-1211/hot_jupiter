      FUNCTION QP (L)
      QP = 1.
      FL = L
   10 IF (FL.LE.0.) GO TO 30
   20 QP = FL*QP/(FL+1.)
      FL = FL-2.
      GO TO  10
   30 RETURN
      END
