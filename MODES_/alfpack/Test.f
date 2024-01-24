      CALL TALFPK(IER)
      STOP
      END
      SUBROUTINE TALFPK(IER)
C Reference: "Tables of Normalized Associated Legendre Polynomials",
C            by S.L. Belousov, Pergamon Press 1962
C
      PARAMETER (N=42, M=15, DEGRS=50.0)
      REAL CP(N/2 + 1)
C
      IER = 0
      CALL ALFK (N,M,CP)
C
      PI = 4.*ATAN(1.0)                                                     
      THETA = DEGRS*PI/180.
      CALL LFPT (N,M,THETA,CP,PB)
C
      WRITE (6,'(A)')      ' Belousov''s tables show PBAR(42,15,50.)'//
     1 '= -0.766020'
      WRITE (6,'(A,F9.6)') ' Routine LFPT in ALFPACK returned '//
     1 'value =',PB
C
      RETURN
      END
