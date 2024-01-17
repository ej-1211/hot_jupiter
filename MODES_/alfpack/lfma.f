C SUBROUTINE LFMA (INIT,N,L,I,PB,W)
C
C DIMENSION OF
C ARGUMENTS              PB(N+1), W(9*L+41)
C
C PURPOSE                ROUTINES LFMA AND LFMB BOTH CALCULATE SINGLE
C                        PRECISION. NORMALIZED ASSOCIATED LEGENDRE
C                        FUNCTIONS PBAR(N,M,THETA) FOR M=0,...,N AT
C                        COLATITUDE THETA=(I-1)*PI/(L-1).  UNDER SPECIAL
C                        CONDITIONS (SEE BELOW) ROUTINE LFMA IS MORE
C                        EFFICIENT THAN LFMB.
C
C USAGE                  CALL LFMA(INIT,N,L,I,PB,W)
C
C ARGUMENTS
C
C ON INPUT               INIT
C                          = 0   INITIALIZATION ONLY
C                          = 1   COMPUTE PBAR(N,M,THETA)
C
C                          LFMA CALL WITH INIT = 0 INITIALIZES ARRAY W;
C                          NO VALUES OF PBAR(N,M,THETA) ARE COMPUTED.
C                          INIT=0 SHOULD BE USED ON THE FIRST CALL, OR
C                          IF N, L, OR W VALUES DIFFER FROM THOSE IN THE
C                          PREVIOUS CALL.
C
C                        N
C                          NONNEGATIVE INTEGER, LESS THAN THE INPUT
C                          PARAMETER L, SPECIFYING THE DEGREE AND
C                          MAXIMUM ORDER OF PBAR(N,M,THETA).
C
C                        L
C                          INTEGER GREATER THAN 1 USED IN SPECIFYING
C                          COLATITUDE THETA=(I-1)*PI/(L-1)
C                          L MUST BE AN ODD INTEGER
C
C                        I
C                          POSITIVE INTEGER, LESS THAN L+1, USED IN
C                          SPECIFYING COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        W
C                          SINGLE PRECISION ARRAY REQUIRING SPECIFIC
C                          VALUES ONLY WHEN THE VALUE OF THE INPUT
C                          PARAMETERS L AND N ARE THE SAME IN THE
C                          PREVIOUS CALL.  IN THIS CASE THE VALUES OF W
C                          MUST BE AS RETURNED IN THE PREVIOUS CALL.
C                          W MUST HAVE 9*L+41 LOCATIONS
C
C
C ON OUTPUT              PB
C                          SINGLE PRECISION STORING PBAR(N,M,THETA)
C                          AS PBAR(M+1), M=0,...,N
C
C                        W
C                          SINGLE PRECISION CONTAINING VALUES WHICH
C                          MUST NOT BE DESTROYED UNTIL THE NEXT CALL
C                          DIFFERS IN VALUE OF INPUT PARAMETERS L OR N.
C
C
C SPECIAL CONDITIONS     ROUTINE LFMA IS MOST EFFICIENT IF THE VALUES
C                        OF THE INPUT PARAMETERS L AND N REMAIN FIXED ON
C                        CONSECUTIVE CALLS
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFMA SOLVES A TRIDIAGONAL SYSTEM OF
C                        EQUATIONS TO OBTAIN VALUES OF PBAR(N,M,THETA).
C
C ACCURACY               COMPARISON BETWEEN ROUTINES LFMA AND DOUBLE
C                        PRECISION DLFMA ON THE CRAY1 INDICATES GREATER
C                        ACCURACY FOR SMALLER VALUES OF INPUT PARAMETER
C                        N.  AGREEMENT TO 12 PLACES WAS OBTAINED FOR
C                        N=10 AND TO 11 PLACES FOR N=100.
C
C TIMING                 TIME PER CALL TO ROUTINE LFMA IS DEPENDENT ON
C                        THE INPUT PARAMETER N.
C
      SUBROUTINE LFMA (INIT,N,L,I,PB,W)
      DIMENSION       PB(1)       ,W(1)
C
      IW1 = L+1
      IW2 = IW1+L
      IW3 = IW2+L
      IW4 = IW3+L
      CALL LFMA1 (INIT,N,L,I,PB,W,W(IW1),W(IW2),W(IW3),W(IW4))
      RETURN
      END
