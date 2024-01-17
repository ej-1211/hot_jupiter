C SUBROUTINE LFP (INIT,N,M,L,CP,PB,W)
C
C DIMENSION OF           CP((N/2)+1), PB(L), W(5*L+41)
C ARGUMENTS
C
C PURPOSE                ROUTINE LFP USES COEFFICIENTS COMPUTED BY
C                        ROUTINE ALFK TO CALCULATE THE SINGLE PRECISION
C                        NORMALIZED ASSOCIATED LEGENDRE FUNCTION PBAR(N,
C                        M,THETA) AT COLATITUDES THETA=(I-1)*PI/(L-1),
C                        I=1,...,L. SUBROUTINE LFP EVALUATES PBAR
C                        USING ONE OF THE FOLLOWING TRIGONOMETRIC
C                        EXPANSIONS
C
C                        1) FOR N EVEN AND M EVEN, PBAR(M,N,THETA) =
C                           .5*CP(1) PLUS THE SUM FROM K=1 TO K=N/2
C                           OF CP(K)*COS(2*K*TH)
C
C                        2) FOR N EVEN AND M ODD, PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=N/2 OF
C                           CP(K)*SIN(2*K*TH)
C
C                        3) FOR N ODD AND M EVEN, PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=(N+1)/2 OF
C                           CP(K)*COS((2*K-1)*TH)
C
C                        4) FOR N ODD AND M ODD,  PBAR(M,N,THETA) =
C                           THE SUM FROM K=1 TO K=(N+1)/2 OF
C                           CP(K)*SIN((2*K-1)*TH)
C
C
C USAGE                  CALL LFP(INIT,N,M,L,CP,PB,W)
C
C ARGUMENTS
C
C ON INPUT               INIT
C                          = 0 INITIALIZATION ONLY
C                          = 1 COMPUTE PBAR(N,M,THETA)
C
C                          LFP CALL WITH INIT = 0 INITIALIZES ARRAY W;
C                          NO VALUES OF PBAR(N,M,THETA) ARE COMPUTED.
C                          INIT=0 SHOULD BE USED ON THE FIRST CALL, OR
C                          IF L OR W VALUES DIFFER FROM THOSE IN THE
C                          PREVIOUS CALL.
C
C                        N
C                          NONNEGATIVE INTEGER, LESS THAN L, SPECIFYING
C                          THE DEGREE OF PBAR(N,M,THETA)
C
C                        M
C                          IS THE ORDER OF PBAR(N,M,THETA). M CAN BE
C                          ANY INTEGER HOWEVER PBAR(N,M,THETA) = 0
C                          IF ABS(M) IS GREATER THAN N AND
C                          PBAR(N,M,THETA) = (-1)**M*PBAR(N,-M,THETA)
C                          FOR NEGATIVE M.
C
C                        L
C                          NUMBER OF COLATITUDES THETA=(I-1)*PI/(L-1)
C                          FOR I=1,...,L WHERE L IS GREATER THAN 1.
C                          L MUST BE AN ODD INTEGER.
C
C                        CP
C                          SINGLE PRECISION ARRAY OF LENGTH (N/2)+1
C                          CONTAINING COEFFICIENTS COMPUTED BY ROUTINE
C                          ALFK
C
C                        W
C                          A SINGLE PRECISION WORK ARRAY WITH AT
C                          LEAST 5*L+41 LOCATIONS
C
C ON OUTPUT              PB
C                          SINGLE PRECISION ARRAY OF LENGTH L CONTAINING
C                          PBAR(N,M,THETA), THETA=(I-1)*PI/(L-1) FOR I=1
C                          ,...,L.
C
C                        W
C                          A SINGLE PRECISION ARRAY CONTAINING VALUES
C                          WHICH MUST NOT BE DESTROYED IF THE NEXT CALL
C                          WILL HAVE THE SAME VALUE OF INPUT PARAMETER N
C
C SPECIAL CONDITIONS     CALLS TO ROUTINE LFP MUST BE PRECEDED BY AN
C                        APPROPRIATE CALL TO ROUTINE ALFK.
C
C PRECISION              SINGLE
C
C ALGORITHM              THE TRIGONOMETRIC SERIES FORMULA USED BY
C                        ROUTINE LFP TO CALCULATE PBAR(N,M,THETA) FOR
C                        THETA=(I-1)*PI/(L-1), I=1,...,N, DEPENDS ON
C                        M AND N AS FOLLOWS:
C
C                           1) FOR N EVEN AND M EVEN, THE FORMULA IS
C                              .5*CP(1) PLUS THE SUM FROM K=1 TO K=N/2
C                              OF CP(K)*COS(2*K*THETA)
C                           2) FOR N EVEN AND M ODD. THE FORMULA IS
C                              THE SUM FROM K=1 TO K=N/2 OF
C                              CP(K)*SIN(2*K*THETA)
C                           3) FOR N ODD AND M EVEN, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*COS((2*K-1)*THETA)
C                           4) FOR N ODD AND M ODD, THE FORMULA IS
C                              THE SUM FROM K=1 TO K=(N+1)/2 OF
C                              CP(K)*SIN((2*K-1)*THETA)
C
C ACCURACY               COMPARISON BETWEEN ROUTINES LFP AND DOUBLE
C                        PRECISION DLFP ON THE CRAY1 INDICATES GREATER
C                        ACCURACY FOR SMALLER VALUES OF INPUT PARAMETER
C                        N.  AGREEMENT TO 12 PLACES WAS OBTAINED FOR
C                        N=10 AND TO 11 PLACES FOR N=100.
C
C TIMING                 TIME PER CALL TO ROUTINE LFP IS DEPENDENT ON
C                        THE INPUT PARAMETERS L AND N.
C
      SUBROUTINE LFP (INIT,N,M,L,CP,PB,W)
      DIMENSION       CP(1)       ,PB(1)    ,W(1)
C
      DO 10 I=1,L
      PB(I) = 0.
   10 CONTINUE
      MA = IABS(M)
      IF(MA .GT. N) RETURN
      IW1 = L+L+12
      IW2 = IW1+3*(L+1)/2+15
      CALL LFP1(INIT,N,MA,L,CP,PB,W,W(IW1),W(IW2))
      RETURN
      END
