C SUBROUTINE LFIN (INIT,THETA,L,M,NM,PB,ID,WLFIN)
C
C DIMENSION OF           THETA(L),  PB(ID,NM+1),  WLFIN(4*L*(NM+1))
C ARGUMENTS
C
C PURPOSE                GIVEN M AND L, ROUTINE LFIN CALCULATES
C                        THE NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
C                        PBAR(N,M,THETA) FOR N=M,...,NM AND THETA(I)
C                        FOR I=1,...,L WHERE
C
C                        PBAR(M,N,THETA) = SQRT((2*N+1)*FACTORIAL(N-M)
C                        /(2*FACTORIAL(N+M)))*SIN(THETA)**M/(2**N*
C                        FACTORIAL(N)) TIMES THE (N+M)TH DERIVATIVE OF
C                        (X**2-1)**N WITH RESPECT TO X=COS(THETA)
C
C USAGE                  CALL LFIN (INIT,THETA,L,M,NM,PB,ID,WLFIN)
C
C ARGUMENTS
C ON INPUT               INIT
C                        = 0
C                            INITIALIZATION ONLY - USING PARAMETERS
C                            L, NM AND THE ARRAY THETA, SUBROUTINE LFIN
C                            INITIALIZES THE ARRAY WLFIN FOR SUBSEQUENT
C                            USE IN THE COMPUTATION OF THE ASSOCIATED
C                            LEGENDRE FUNCTIONS PB. INITIALIZATION DOES
C                            NOT HAVE TO BE REPEATED UNLESS L, NM OR
C                            THE ARRAY THETA ARE CHANGED.
C                        = 1
C                            SUBROUTINE LFIN USES THE ARRAY WLFIN THAT
C                            WAS COMPUTED WITH INIT = 0 TO COMPUTE PB
C
C                        THETA
C                          AN ARRAY THAT CONTAINS THE COLATITUDES
C                          AT WHICH THE ASSOCIATED LEGENDRE FUNCTIONS
C                          WILL BE COMPUTED. THE COLATITUDES MUST BE
C                          SPECIFIED IN RADIANS.
C
C                        L
C                          THE LENGTH OF THE THETA ARRAY. LFIN IS
C                          VECTORIZED WITH VECTOR LENGTH L.
C
C                        M
C                          NONNEGATIVE INTEGER, LESS THAN NM, SPECIFYING
C                          DEGREE OF PBAR(N,M,THETA). SUBROUTINE LFIN
C                          MUST BE CALLED STARTING WITH N=0. N MUST BE
C                          INCREMENTED BY ONE IN SUBSEQUENT CALLS AND
C                          MUST NOT EXCEED NM.
C
C                        NM
C                          THE MAXIMUM VALUE OF N AND M
C
C                        ID
C                          THE FIRST DIMENSION OF THE TWO DIMENSIONAL
C                          ARRAY PB AS IT APPEARS IN THE PROGRAM THAT
C                          CALLS LFIN. (SEE OUTPUT PARAMETER PB)
C
C                        WLFIN
C                          AN ARRAY WITH LENGTH 4*L*(NM+1) WHICH
C                          MUST BE INITIALIZED BY CALLING LFIN
C                          WITH INIT=0 (SEE PARAMETER INIT)  IT
C                          MUST NOT BE ALTERED BETWEEN CALLS TO
C                          LFIN.
C
C
C ON OUTPUT              PB
C                          A TWO DIMENSIONAL ARRAY WITH FIRST
C                          DIMENSION ID IN THE PROGRAM THAT CALLS
C                          LFIN. THE SECOND DIMENSION OF PB MUST
C                          BE AT LEAST NM+1. STARTING WITH M=0
C                          LFIN IS CALLED REPEATEDLY WITH M BEING
C                          INCREASED BY ONE BETWEEN CALLS. ON EACH
C                          CALL, SUBROUTINE LFIN COMPUTES PB(I,N+1)
C                          = PBAR(M,N,THETA(I)) FOR N=M,...,NM AND
C                          I=1,...L.
C
C                        WLFIN
C                          ARRAY CONTAINING VALUES WHICH MUST NOT
C                          BE ALTERED UNLESS L, NM OR THE ARRAY THETA
C                          ARE CHANGED IN WHICH CASE LFIN MUST BE
C                          CALLED WITH INIT=0 TO REINITIALIZE THE
C                          WLFIN ARRAY.
C
C SPECIAL CONDITIONS     M MUST BE INCREASED BY ONE BETWEEN CALLS
C                        OF LFIN IN WHICH M IS NOT ZERO.
C
C PRECISION              SINGLE
C
C ALGORITHM              ROUTINE LFIN CALCULATES PBAR(N,M,THETA) USING
C                        A FOUR TERM RECURRENCE RELATION. (UNPUBLISHED
C                        NOTES BY PAUL N. SWARZTRAUBER)
C
      SUBROUTINE LFIN (INIT,THETA,L,M,NM,PB,ID,WLFIN)
      DIMENSION       PB(1)        ,WLFIN(1)
C
C     TOTAL LENGTH OF WLFIN IS 4*L*(NM+1)
C
      LNX = L*(NM+1)
      IW1 = LNX+1
      IW2 = IW1+LNX
      IW3 = IW2+LNX
      CALL LFIN1(INIT,THETA,L,M,NM,ID,PB,WLFIN,WLFIN(IW1),
     1                WLFIN(IW2),WLFIN(IW3),WLFIN(IW2))
      RETURN
      END
