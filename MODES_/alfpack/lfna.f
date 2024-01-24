C SUBROUTINE LFNA (INIT,M,L,I,PB,W)
C
C DIMENSION OF           W(9*L+41), PB(L)
C ARGUMENTS
C
C PURPOSE                ROUTINES LFNA AND LFNB BOTH CALCULATE SINGLE
C                        PRECISION, NORMALIZED ASSOCIATED LEGENDRE
C                        FUNCTIONS PBAR(N,M,THETA) FOR N=M,...,L-1 AT
C                        COLATITUDE THETA=(I-1)*PI/(L-1).  UNDER SPECIAL
C                        CONDITIONS (SEE BELOW) ROUTINE LFNA IS MORE
C                        EFFICIENT THAN LFNB.
C
C USAGE                  CALL LFNA(INIT,M,L,I,PB,W)
C
C ARGUMENTS
C
C ON INPUT               INIT
C                          = 0   INITIALIZATION ONLY
C                          = 1   COMPUTE PBAR(N,M,THETA)
C
C                          LFNA CALL WITH INIT = 0 INITIALIZES ARRAY W;
C                          NO VALUES OF PBAR(N,M,THETA) ARE COMPUTED.
C                          INIT=0 SHOULD BE USED ON THE FIRST CALL, OR
C                          IF L, M, OR W VALUES DIFFER FROM THOSE IN THE
C                          PREVIOUS CALL.
C
C                        M
C                          NONNEGATIVE INTEGER, LESS THAN L, SPECIFYING
C                          THE ORDER OF PBAR(N,M,THETA), N=M,...,L-1.
C
C                        L
C                          INTEGER GREATER THAN 1 USED IN SPECIFYING
C                          THE DEGREE OF PBAR(N,M,THETA), N=M,...,L-1,
C                          AND COLATITUDE THETA=(I-1)*PI/(L-1)
C                          L MUST BE AN ODD INTEGER
C
C                        I
C                          POSITIVE INTEGER, LESS THAN L+1, USED IN
C                          SPECIFYING COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        W
C                          SINGLE PRECISION ARRAY REQUIRING SPECIFIC
C                          VALUES ONLY WHEN THE VALUE OF THE INPUT
C                          PARAMETERS L AND M ARE THE SAME IN THE
C                          PREVIOUS CALL.  IN THIS CASE THE VALUES OF
C                          W MUST BE AS RETURNED IN THE PREVIOUS CALL.
C                          W MUST HAVE 9*L+41 LOCATIONS
C
C ON OUTPUT              PB
C                          SINGLE PRECISION ARRAY STORING PBAR(N,M,THETA
C                          IN PB(N+1) , N=M,...,L-1
C
C                        W
C                          SINGLE PRECISION ARRAY CONTAINING VALUES
C                          WHICH MUST NOT BE DESTROYED UNTIL THE NEXT
C                          CALL DIFFERS IN VALUE OF INPUT PARAMETERS
C                          L OR M
C
C SPECIAL CONDITIONS     ROUTINE LFNA IS MOST EFFICIENT IF THE VALUES
C                        OF THE INPUT PARAMETERS L AND M REMAIN FIXED ON
C                        CONSECUTIVE CALLS
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFNA SOLVES A TRIDIAGONAL SYSTEM OF
C                        EQUATIONS TO OBTAIN VALUES OF PBAR(N,M,THETA)
C
C ACCURACY               COMPARISON BETWEEN ROUTINES LFNA AND DOUBLE
C                        PRECISION DLFNA ON THE CRAY1 INDICATES GREATER
C                        ACCURACY FOR SMALLER VALUES OF INPUT PARAMETER
C                        M.  AGREEMENT TO 12 PLACES WAS OBTAINED FOR
C                        M=10 AND TO 11 PLACES FOR M=100.
C
C TIMING                 TIME PER CALL TO ROUTINE LFNA IS DEPENDENT ON
C                        THE DIFFERENCE L-M OF THE INPUT PARAMETERS L
C                        AND M.
C
      SUBROUTINE LFNA (INIT,M,L,I,PB,W)
      DIMENSION       PB(1)       ,W(1)
C
      IW1 = L+1
      IW2 = IW1+L
      IW3 = IW2+L
      IW4 = IW3+L
      CALL LFNA1 (INIT,M,L,I,PB,W,W(IW1),W(IW2),W(IW3),W(IW4))
      RETURN
      END
