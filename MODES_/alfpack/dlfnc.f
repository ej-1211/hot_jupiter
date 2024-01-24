C SUBROUTINE LFNC (INIT,M,L,THETA,PB,W)
C
C DIMENSION OF           W(5*L), PB(L+1)
C ARGUMENTS
C
C PURPOSE                ROUTINES LFNA AND LFNC BOTH CALCULATE SINGLE
C                        PRECISION, NORMALIZED ASSOCIATED LEGENDRE
C                        FUNCTIONS PBAR(N,M,THETA) FOR N=M,...,L-1
C                        SUBROUTINE LFNC DIFFERS FROM LFNA IN THAT THETA
C                        IS USER SPECIFIED RATHER THAN EQUALLY SPACED
C
C USAGE                  CALL LFNC(INIT,M,L,THETA,PB,W)
C
C ARGUMENTS
C
C ON INPUT               INIT
C                          = 0   INITIALIZATION ONLY
C                          = 1   COMPUTE PBAR(N,M,THETA)
C
C                          LFNC CALL WITH INIT = 0 INITIALIZES ARRAY W;
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
C                          THE MAXIMUM VALUE OF N+1, WHERE M IS DEFINED
C                          IN THE 'PURPOSE' SECTION ABOVE.
C
C                        THETA
C                          COLATITUDE, IN RADIANS, FOR PBAR(N,M,THETA).
C
C                        W
C                          SINGLE PRECISION ARRAY REQUIRING SPECIFIC
C                          VALUES ONLY WHEN THE VALUE OF THE INPUT
C                          PARAMETERS L AND M ARE THE SAME IN THE
C                          PREVIOUS CALL.  IN THIS CASE THE VALUES OF
C                          W MUST BE AS RETURNED IN THE PREVIOUS CALL.
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
C SPECIAL CONDITIONS     ROUTINE LFNC IS MOST EFFICIENT IF THE VALUES
C                        OF THE INPUT PARAMETERS L AND M REMAIN FIXED ON
C                        CONSECUTIVE CALLS
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFNC SOLVES A TRIDIAGONAL SYSTEM OF
C                        EQUATIONS TO OBTAIN VALUES OF PBAR(N,M,THETA)
C
C ACCURACY               THE ACCURACY IS GREATER FOR THE SMALLER VALUES
C                        OF INPUT PARAMETER M.  AGREEMENT TO 12 PLACES
C                        WAS OBTAINED FOR M=10 AND TO 11 PLACES FOR
C                        M=100
C
C TIMING                 TIME PER CALL TO ROUTINE LFNC IS DEPENDENT ON
C                        THE DIFFERENCE L-M OF THE INPUT PARAMETERS L
C                        AND M.
C
      SUBROUTINE DLFNC (INIT,M,L,THETA,PB,W)
      DOUBLE PRECISION  PB        ,W
      DIMENSION         PB(1)     ,W(1)
C
      IW1 = L+1
      IW2 = IW1+L
      IW3 = IW2+L
      IW4 = IW3+L
      CALL DLFNC1(INIT,M,L,THETA,PB,W,W(IW1),W(IW2),W(IW3),W(IW4))
      RETURN
      END
