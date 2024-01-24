C SUBROUTINE LFMB (INIT,N,L,I,PB,W)
C
C DIMENSION OF
C ARGUMENTS              W(4*L), PB(L)
C
C PURPOSE                ROUTINES LFMB AND LFMA BOTH CALCULATE SINGLE
C                        PRECISION, NORMALIZED ASSOCIATED LEGENDRE
C                        FUNCTIONS PBAR(N,M,THETA) FOR M=0,...,N AT
C                        COLATITUDE THETA=(I-1)*PI/(L-1).  UNDER SPECIAL
C                        CONDITIONS (SEE BELOW) ROUTINE LFMB IS MORE
C                        EFFICIENT THAN LFMA.
C
C USAGE                  CALL LFMB(INIT,N,L,I,PB,W)
C
C ARGUMENTS
C ON INPUT               INIT
C                          = 0   INITIALIZATION ONLY
C                          = 1   COMPUTE PBAR(N,M,THETA)
C
C                          LFMA CALL WITH INIT = 0 INITIALIZES ARRAY W;
C                          NO VALUES OF PBAR(N,M,THETA) ARE COMPUTED.
C                          INIT=0 SHOULD BE USED ON THE FIRST CALL, OR
C                          IF L, I, OR W VALUES DIFFER FROM THOSE IN THE
C                          PREVIOUS CALL.
C
C                        N
C                          NONNEGATIVE INTEGER, LESS THAN L, SPECIFYING
C                          THE DEGREE AND MIXIMUM ORDER OF PBAR(N,M,
C                          THETA)
C
C                        L
C                          INTEGER GREATER THAN 1 USED IN SPECIFYING
C                          COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        I
C                          POSITIVE INTEGER, LESS THAN L+1, USED IN
C                          SPECIFYING COLATITUDE THETA=(I-1)*PI/(L-1)
C
C                        W
C                          SINGLE PRECISION ARRAY REQUIRING SPECIFIC
C                          VALUES ONLY WHEN THE VALUE OF THE INPUT
C                          PARAMETERS I AND L ARE THE SAME IN THE
C                          PREVIOUS CALL.  IN THIS CASE THE VALUES OF W
C                          MUST BE AS RETURNED IN THE PREVIOUS CALL.
C
C ON OUTPUT              PB
C                          SINGLE ARRAY STORING PBAR(N,M,THETA) IN
C                          PB(M+1), M=0,...,N
C
C                        W
C                          SINGLE PRECISION ARRAY CONTAINING VALUES
C                          WHICH MUST NOT BE DESTROYED UNTIL THE NEXT
C                          CALL DIFFERS IN VALUE OF INPUT PARAMETERS I
C                          OR L
C
C SPECIAL CONDITIONS     ROUTINE LFMB IS MOST EFFICIENT IF THE VALUES
C                        OF THE INPUT PARAMETERS I AND L REMAIN FIXED
C                        ON CONSECUTIVE CALLS.
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFMB SOLVES A TRIDIAGONAL SYSTEM OF
C                        EQUATIONS TO OBTAIN VALUES OF PBAR(N,M,THETA).
C
C TIMING                 TIME PER CALL TO ROUTINE LFMB IS DEPENDENT ON
C                        THE INPUT PARAMETER N.
C
      SUBROUTINE LFMB (INIT,N,L,I,PB,W)
      DIMENSION       PB(1)       ,W(1)
C
      IW1 = L+1
      IW2 = IW1+L
      IW3 = IW2+L
      CALL LFMB1 (INIT,L,N,I,PB,W,W,W(IW1),W(IW2),W(IW2),W(IW3))
      RETURN
      END
