C SUBROUTINE LFNB (INIT,M,L,I,PB,W)
C
C DIMENSION OF
C ARGUMENTS              W(4*L), PB(L)
C
C PURPOSE                ROUTINES LFNB AND LFNA BOTH CALCULATE SINGLE
C                        PRECISION, NORMALIZED ASSOCIATED LEGENDRE
C                        FUNCTIONS PBAR(N,M,THETA) FOR N=M,...,L-1 AT
C                        COLATITUDE THETA=(I-1)*PI/(L-1).  UNDER SPECIAL
C                        CONDITIONS (SEE BELOW) ROUTINE LFNB IS MORE
C                        EFFICIENT THAN LFNA.
C
C USAGE                  CALL LFNB (INIT,M,L,I,PB,W)
C
C ARGUMENTS
C ON INPUT               INIT
C                          = 0   INITIALIZATION ONLY
C                          = 1   COMPUTE PBAR(N,M,THETA)
C
C                          LFNB CALL WITH INIT = 0 INITIALIZES ARRAY W;
C                          NO VALUES OF PBAR(N,M,THETA) ARE COMPUTED.
C                          INIT=0 SHOULD BE USED ON THE FIRST CALL, OR
C                          IF I, L, OR W VALUES DIFFER FROM THOSE IN THE
C                          PREVIOUS CALL.
C
C                        M
C                          NONNEGATIVE INTEGER, LESS THAN L, SPECIFYING
C                          THE ORDER OF PBAR(N,M,THETA), N=M,...,L-1.
C                        L
C                          INTEGER GREATER THAN 1 USED IN SPECIFYING
C                          THE DEGREE OF PBAR(N,M,THETA), N=M,...,L-1,
C                          AND COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        I
C                          POSITIVE INTEGER, LESS THAN L+1, USED IN
C                          SPECIFYING COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        W
C                          SINGLE PRECISION ARRAY REQUIRING SPECIFIC
C                          VALUES ONLY WHEN THE VALUE OF THE INPUT
C                          PARAMETERS L AND I ARE THE SAME IN THE
C                          PREVIOUS CALL.  IN THIS CASE THE VALUES OF W
C                          MUST BE AS RETURNED IN THE PREVIOUS CALL.
C
C ON OUTPUT              PB
C                          SINGLE PRECISION ARRAY STORING PBAR(N,M,
C                          THETA) IN PB(N+1), N=M,...,L-1
C
C                        W
C                          SINGLE PRECISION ARRAY CONTAINING VALUES
C                          WHICH MUST NOT BE DESTROYED UNTIL THE NEXT
C                          CALL DIFFERS IN VALUE OF INPUT PARAMETERS L
C                          OR I
C
C SPECIAL CONDITIONS     ROUTINE LFNB IS MOST EFFICIENT IF M IS THE ONLY
C                        INPUT PARAMETER VALUE CHANGING IN CONSECUTIVE
C                        CALLS.
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFNB SOLVES A TRIDIAGONAL SYSTEM OF
C                        EQUATIONS TO OBTAIN VALUES OF PBAR(N,M,THETA)
C
C ACCURACY               COMPARISON BETWEEN ROUTINES LFNB AND DOUBLE
C                        PRECISION DLFNB ON THE CRAY1 INDICATE GREATER
C                        ACCURACY FOR SMALLER VALUES OF INPUT PARAMETER
C                        M.  AGREEMENT TO 11 PLACES WAS OBTAINED FOR
C                        M=10 AND TO 10 PLACES FOR M=100.
C
C TIMING                 TIME PER CALL TO ROUTINE LFNB IS DEPENDENT ON
C                        THE DIFFERENCE L-M OF THE INPUT PARAMETERS L
C                        AND M.
C
      SUBROUTINE LFNB (INIT,M,L,I,PB,W)
      DIMENSION       PB(1)       ,W(1)
C
      IW1 = L+1
      IW2 = IW1+L
      IW3 = IW2+L
      CALL LFNB1 (INIT,M,L,I,PB,W,W(IW1),W(IW2),W(IW3))
      RETURN
      END
