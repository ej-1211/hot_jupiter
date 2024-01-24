C SUBROUTINE LFIM (INIT,THETA,L,N,NM,PB,ID,WLFIM)
C
C DIMENSION OF           THETA(L),  PB(ID,NM+1),  WLFIM(4*L*(NM+1))
C ARGUMENTS
C
C PURPOSE                GIVEN N AND L, ROUTINE LFIM CALCULATES
C                        THE NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS
C                        PBAR(N,M,THETA) FOR M=0,...,N AND THETA(I)
C                        FOR I=1,...,L WHERE
C
C                        PBAR(M,N,THETA) = SQRT((2*N+1)*FACTORIAL(N-M)
C                        /(2*FACTORIAL(N+M)))*SIN(THETA)**M/(2**N*
C                        FACTORIAL(N)) TIMES THE (N+M)TH DERIVATIVE OF
C                        (X**2-1)**N WITH RESPECT TO X=COS(THETA)
C
C USAGE                  CALL LFIM (INIT,THETA,L,N,NM,PB,ID,WLFIM)
C
C ARGUMENTS
C ON INPUT               INIT
C                        = 0
C                            INITIALIZATION ONLY - USING PARAMETERS
C                            L, NM AND ARRAY THETA, SUBROUTINE LFIM
C                            INITIALIZES ARRAY WLFIM FOR SUBSEQUENT
C                            USE IN THE COMPUTATION OF THE ASSOCIATED
C                            LEGENDRE FUNCTIONS PB. INITIALIZATION
C                            DOES NOT HAVE TO BE REPEATED UNLESS
C                            L, NM, OR ARRAY THETA ARE CHANGED.
C                        = 1
C                            SUBROUTINE LFIM USES THE ARRAY WLFIM THAT
C                            WAS COMPUTED WITH INIT = 0 TO COMPUTE PB.
C
C                        THETA
C                          AN ARRAY THAT CONTAINS THE COLATITUDES
C                          AT WHICH THE ASSOCIATED LEGENDRE FUNCTIONS
C                          WILL BE COMPUTED. THE COLATITUDES MUST BE
C                          SPECIFIED IN RADIANS.
C
C                        L
C                          THE LENGTH OF THE THETA ARRAY. LFIM IS
C                          VECTORIZED WITH VECTOR LENGTH L.
C
C                        N
C                          NONNEGATIVE INTEGER, LESS THAN NM, SPECIFYING
C                          DEGREE OF PBAR(N,M,THETA). SUBROUTINE LFIM
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
C                          CALLS LFIM. (SEE OUTPUT PARAMETER PB)
C
C                        WLFIM
C                          AN ARRAY WITH LENGTH 4*L*(NM+1) WHICH
C                          MUST BE INITIALIZED BY CALLING LFIM
C                          WITH INIT=0 (SEE PARAMETER INIT)  IT
C                          MUST NOT BE ALTERED BETWEEN CALLS TO
C                          LFIM.
C
C
C ON OUTPUT              PB
C                          A TWO DIMENSIONAL ARRAY WITH FIRST
C                          DIMENSION ID IN THE PROGRAM THAT CALLS
C                          LFIM. THE SECOND DIMENSION OF PB MUST
C                          BE AT LEAST NM+1. STARTING WITH N=0
C                          LFIM IS CALLED REPEATEDLY WITH N BEING
C                          INCREASED BY ONE BETWEEN CALLS. ON EACH
C                          CALL, SUBROUTINE LFIM COMPUTES
C                          = PBAR(M,N,THETA(I)) FOR M=0,...,N AND
C                          I=1,...L.
C
C                        WLFIM
C                          ARRAY CONTAINING VALUES WHICH MUST NOT
C                          BE ALTERED UNLESS L, NM OR THE ARRAY THETA
C                          ARE CHANGED IN WHICH CASE LFIM MUST BE
C                          CALLED WITH INIT=0 TO REINITIALIZE THE
C                          WLFIM ARRAY.
C
C SPECIAL CONDITIONS     N MUST BE INCREASED BY ONE BETWEEN CALLS
C                        OF LFIM IN WHICH N IS NOT ZERO.
C
C PRECISION              SINGLE
C
C
C ALGORITHM              ROUTINE LFIM CALCULATES PBAR(N,M,THETA) USING
C                        A FOUR TERM RECURRENCE RELATION. (UNPUBLISHED
C                        NOTES BY PAUL N. SWARZTRAUBER)
C
      SUBROUTINE LFIM (INIT,THETA,L,N,NM,PB,ID,WLFIM)
      DIMENSION       PB(1)        ,WLFIM(1)
C
C     TOTAL LENGTH OF WLFIM IS 4*L*(NM+1)
C
      LNX = L*(NM+1)
      IW1 = LNX+1
      IW2 = IW1+LNX
      IW3 = IW2+LNX
      CALL LFIM1(INIT,THETA,L,N,NM,ID,PB,WLFIM,WLFIM(IW1),
     1                WLFIM(IW2),WLFIM(IW3),WLFIM(IW2))
      RETURN
      END
