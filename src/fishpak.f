C
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION 3.1 , OCTOBER 1980)                  *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     cosgen.f

*DECK COSGEN
      SUBROUTINE COSGEN (N, IJUMP, FNUM, FDEN, A)
C***BEGIN PROLOGUE  COSGEN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (COSGEN-S, CMPCSG-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine computes required cosine values in ascending
C     order.  When IJUMP .GT. 1 the routine computes values
C
C        2*COS(J*PI/L) , J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1)
C
C     where L = IJUMP*(N/IJUMP+1).
C
C
C     when IJUMP = 1 it computes
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     where
C        FNUM = 0.5, FDEN = 0.0, for regular reduction values.
C        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C                                in POISN2 only.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  COSGEN
      DIMENSION       A(*)
C
C
C***FIRST EXECUTABLE STATEMENT  COSGEN
      PI = PIMACH(DUM)
      IF (N .EQ. 0) GO TO 105
      IF (IJUMP .EQ. 1) GO TO 103
      K3 = N/IJUMP+1
      K4 = K3-1
      PIBYN = PI/(N+IJUMP)
      DO 102 K=1,IJUMP
         K1 = (K-1)*K3
         K5 = (K-1)*K4
         DO 101 I=1,K4
            X = K1+I
            K2 = K5+I
            A(K2) = -2.*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
      GO TO 105
  103 CONTINUE
      NP1 = N+1
      Y = PI/(N+FDEN)
      DO 104 I=1,N
         X = NP1-I-FNUM
         A(I) = 2.*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
      RETURN
      END

!============================================================================

c     genbun.f

*DECK GENBUN
      SUBROUTINE GENBUN (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
     +   IERROR, W)
C***BEGIN PROLOGUE  GENBUN
C***PURPOSE  Solve by a cyclic reduction algorithm the linear system
C            of equations that results from a finite difference
C            approximation to certain 2-d elliptic PDE's on a centered
C            grid .
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C)
C***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine GENBUN solves the linear system of equations
C
C          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C
C          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C               for I = 1,2,...,M  and  J = 1,2,...,N.
C
C     The indices I+1 and I-1 are evaluated modulo M, i.e.,
C     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
C     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or
C     X(I,1) depending on an input parameter.
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     NPEROD
C       Indicates the values that X(I,0) and X(I,N+1) are assumed to
C       have.
C
C       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1).
C       = 1  If X(I,0) = X(I,N+1) = 0  .
C       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1).
C       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1).
C       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0.
C
C     N
C       The number of unknowns in the J-direction.  N must be greater
C       than 2.
C
C     MPEROD
C       = 0 if A(1) and C(M) are not zero.
C       = 1 if A(1) = C(M) = 0.
C
C     M
C       The number of unknowns in the I-direction.  M must be greater
C       than 2.
C
C     A,B,C
C       One-dimensional arrays of length M that specify the
C       coefficients in the linear equations given above.  If MPEROD = 0
C       the array elements must not depend upon the index I, but must be
C       constant.  Specifically, the subroutine checks the following
C       condition
C
C             A(I) = C(1)
C             C(I) = C(1)
C             B(I) = B(1)
C
C       for I=1,2,...,M.
C
C     IDIMY
C       The row (or first) dimension of the two-dimensional array Y as
C       it appears in the program calling GENBUN.  This parameter is
C       used to specify the variable dimension of Y.  IDIMY must be at
C       least M.
C
C     Y
C       A two-dimensional array that specifies the values of the right
C       side of the linear system of equations given above.  Y must be
C       dimensioned at least M*N.
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M
C       locations.  The actual number of locations used is computed by
C       GENBUN and is returned in location W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     Y
C       Contains the solution X.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for number zero, a solution is not attempted.
C
C       = 0  No error.
C       = 1  M .LE. 2
C       = 2  N .LE. 2
C       = 3  IDIMY .LT. M
C       = 4  NPEROD .LT. 0 or NPEROD .GT. 4
C       = 5  MPEROD .LT. 0 or MPEROD .GT. 1
C       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for
C            some I=1,2,...,M.
C       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list)
C     Arguments
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3,
C     Required       PIMACH
C
C     Special        NONE
C     Conditions
C
C     Common         NONE
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Standardized April 1, 1973
C                    Revised August 20,1973
C                    Revised January 1, 1976
C
C     Algorithm      The linear system is solved by a cyclic reduction
C                    algorithm described in the reference.
C
C     Space          4944(decimal) = 11520(octal) locations on the NCAR
C     Required       Control Data 7600.
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine GENBUN is roughly proportional
C                    to M*N*log2(N), but also depends on the input
C                    parameter NPEROD.  Some typical values are listed
C                    in the table below.  More comprehensive timing
C                    charts may be found in the reference.
C                       To measure the accuracy of the algorithm a
C                    uniform random number generator was used to create
C                    a solution array X for the system given in the
C                    'PURPOSE' with
C
C                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
C
C                    and, when MPEROD = 1
C
C                       A(1) = C(M) = 0
C                       A(M) = C(1) = 2.
C
C                    The solution X was substituted into the given sys-
C                    tem and, using double precision, a right side Y was
C                    computed.  Using this array Y subroutine GENBUN was
C                    called to produce an approximate solution Z.  Then
C                    the relative error, defined as
C
C                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C
C                    where the two maxima are taken over all I=1,2,...,M
C                    and J=1,2,...,N, was computed.  The value of E is
C                    given in the table below for some typical values of
C                    M and N.
C
C
C                       M (=N)    MPEROD    NPEROD    T(MSECS)    E
C                       ------    ------    ------    --------  ------
C
C                         31        0         0          36     6.E-14
C                         31        1         1          21     4.E-13
C                         31        1         3          41     3.E-13
C                         32        0         0          29     9.E-14
C                         32        1         1          32     3.E-13
C                         32        1         3          48     1.E-13
C                         33        0         0          36     9.E-14
C                         33        1         1          30     4.E-13
C                         33        1         3          34     1.E-13
C                         63        0         0         150     1.E-13
C                         63        1         1          91     1.E-12
C                         63        1         3         173     2.E-13
C                         64        0         0         122     1.E-13
C                         64        1         1         128     1.E-12
C                         64        1         3         199     6.E-13
C                         65        0         0         143     2.E-13
C                         65        1         1         120     1.E-12
C                         65        1         3         138     4.E-13
C
C     Portability    American National Standards Institute Fortran.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For
C                    Solving Block Tridiagonal Systems Of Arbitrary
C                    Dimensions,' SIAM J. on Numer. Anal.,
C                    14(Sept., 1977), PP. 706-720.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
C                 block tridiagonal systems of arbitrary dimensions,
C                 SIAM Journal on Numerical Analysis 14, (September
C                 1977), pp. 706-720.
C***ROUTINES CALLED  POISD2, POISN2, POISP2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  GENBUN
C
C
      DIMENSION       Y(IDIMY,*)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C***FIRST EXECUTABLE STATEMENT  GENBUN
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.0 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 102
      DO 101 I=2,M
         IF (A(I) .NE. C(1)) GO TO 103
         IF (C(I) .NE. C(1)) GO TO 103
         IF (B(I) .NE. B(1)) GO TO 103
  101 CONTINUE
      GO TO 104
  102 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
      GO TO 104
  103 IERROR = 6
  104 IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      IWBA = MP1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      MP = MPEROD+1
      NP = NPEROD+1
      GO TO (114,107),MP
  107 GO TO (108,109,110,111,123),NP
  108 CALL POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  109 CALL POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),
     1             W(IWD),W(IWTCOS),W(IWP))
      GO TO 112
  110 CALL POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  111 CALL POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
  112 IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 124
  113 GO TO (127,133),MP
  114 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 119 J=1,N
         DO 115 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116),MODD
  116    W(M) = 2.*Y(M,J)
  117    CONTINUE
         DO 118 I=1,M
            Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (120,121),MODD
  120 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 IREV = 1
      NBY2 = N/2
  124 DO 126 J=1,NBY2
         MSKIP = N+1-J
         DO 125 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
      GO TO (110,113),IREV
  127 CONTINUE
      DO 132 J=1,N
         DO 128 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (130,129),MODD
  129    W(M) = .5*Y(M,J)
  130    CONTINUE
         DO 131 I=1,M
            Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END

!============================================================================

c     hwscrt.f

*DECK HWSCRT
      SUBROUTINE HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSCRT
C***PURPOSE  Solves the standard five-point finite difference
C            approximation to the Helmholtz equation in Cartesian
C            coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSCRT-S)
C***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSCRT solves the standard five-point finite
C     difference approximation to the Helmholtz equation in Cartesian
C     coordinates:
C
C          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y).
C
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     A,B
C       The range of X, i.e., A .LE. X .LE. B.  A must be less than B.
C
C     M
C       The number of panels into which the interval (A,B) is
C       subdivided.  Hence, there will be M+1 grid points in the
C       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1,
C       where DX = (B-A)/M is the panel width. M must be greater than 3.
C
C     MBDCND
C       Indicates the type of boundary conditions at X = A and X = B.
C
C       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J).
C       = 1  If the solution is specified at X = A and X = B.
C       = 2  If the solution is specified at X = A and the derivative of
C            the solution with respect to X is specified at X = B.
C       = 3  If the derivative of the solution with respect to X is
C            specified at X = A and X = B.
C       = 4  If the derivative of the solution with respect to X is
C            specified at X = A and the solution is specified at X = B.
C
C     BDA
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to X at X = A.
C       When MBDCND = 3 or 4,
C
C            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDA is a dummy variable.
C
C     BDB
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to X at X = B.
C       When MBDCND = 2 or 3,
C
C            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value BDB is a dummy variable.
C
C     C,D
C       The range of Y, i.e., C .LE. Y .LE. D.  C must be less than D.
C
C     N
C       The number of panels into which the interval (C,D) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where
C       DY = (D-C)/N is the panel width.  N must be greater than 3.
C
C     NBDCND
C       Indicates the type of boundary conditions at Y = C and Y = D.
C
C       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J).
C       = 1  If the solution is specified at Y = C and Y = D.
C       = 2  If the solution is specified at Y = C and the derivative of
C            the solution with respect to Y is specified at Y = D.
C       = 3  If the derivative of the solution with respect to Y is
C            specified at Y = C and Y = D.
C       = 4  If the derivative of the solution with respect to Y is
C            specified at Y = C and the solution is specified at Y = D.
C
C     BDC
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Y at Y = C.
C       When NBDCND = 3 or 4,
C
C            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDC is a dummy variable.
C
C     BDD
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Y at Y = D.
C       When NBDCND = 2 or 3,
C
C            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDD is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .GT. 0, a solution may not exist.  However, HWSCRT will
C       attempt to find a solution.
C
C     F
C       A two-dimensional array which specifies the values of the right
C       side of the Helmholtz equation and boundary values (if any).
C       For I = 2,3,...,M and J = 2,3,...,N
C
C            F(I,J) = F(X(I),Y(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND     F(1,J)        F(M+1,J)
C            ------     ---------     --------
C
C              0        F(A,Y(J))     F(A,Y(J))
C              1        U(A,Y(J))     U(B,Y(J))
C              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1
C              3        F(A,Y(J))     F(B,Y(J))
C              4        F(A,Y(J))     U(B,Y(J))
C
C
C            NBDCND     F(I,1)        F(I,N+1)
C            ------     ---------     --------
C
C              0        F(X(I),C)     F(X(I),C)
C              1        U(X(I),C)     U(X(I),D)
C              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1
C              3        F(X(I),C)     F(X(I),D)
C              4        F(X(I),C)     U(X(I),D)
C
C       F must be dimensioned at least (M+1)*(N+1).
C
C       NOTE:
C
C       If the table calls for both the solution U and the right side F
C       at a corner then the solution must be specified.
C
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSCRT.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*(N+1) +
C       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
C       locations used is computed by HWSCRT and is returned in location
C       W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1,
C       J = 1,2,...,N+1  .
C
C     PERTRB
C       If a combination of periodic or derivative boundary conditions
C       is specified for a Poisson equation (LAMBDA = 0), a solution may
C       not exist.  PERTRB is a constant, calculated and subtracted from
C       F, which ensures that a solution exists.  HWSCRT then computes
C       this solution, which is a least squares solution to the original
C       approximation.  This solution plus any constant is also a
C       solution.  Hence, the solution is not unique.  The value of
C       PERTRB should be small compared to the right side F.  Otherwise,
C       a solution is obtained to an essentially different problem.
C       This comparison should always be made to insure that a
C       meaningful solution has been obtained.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for numbers 0 and 6, a solution is not attempted.
C
C       = 0  No error.
C       = 1  A .GE. B.
C       = 2  MBDCND .LT. 0 or MBDCND .GT. 4  .
C       = 3  C .GE. D.
C       = 4  N .LE. 3
C       = 5  NBDCND .LT. 0 or NBDCND .GT. 4  .
C       = 6  LAMBDA .GT. 0  .
C       = 7  IDIMF .LT. M+1  .
C       = 8  M .LE. 3
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSCRT, the user should test IERROR after the call.
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C
C     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
C     Arguments      W(see argument list)
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
C     Required       TRIX,TRI3,PIMACH
C
C     Special        NONE
C     Conditions
C
C     Common         NONE
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Standardized September 1, 1973
C                    Revised April 1, 1976
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    GENBUN to solve the system.
C
C     Space          13110(octal) = 5704(decimal) locations on the NCAR
C     Required       Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HWSCRT is roughly proportional
C                    to M*N*log2(N), but also depends on the input
C                    parameters NBDCND and MBDCND.  Some typical values
C                    are listed in the table below.
C                       The solution process employed results in a loss
C                    of no more than three significant digits for N and
C                    M as large as 64.  More detailed information about
C                    accuracy can be found in the documentation for
C                    subroutine GENBUN which is the routine that
C                    solves the finite difference equations.
C
C
C                       M(=N)    MBDCND    NBDCND    T(MSECS)
C                       -----    ------    ------    --------
C
C                        32        0         0          31
C                        32        1         1          23
C                        32        3         3          36
C                        64        0         0         128
C                        64        1         1          96
C                        64        3         3         142
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
C                    Subprograms for The Solution Of Elliptic Equations'
C                    NCAR TN/IA-109, July, 1975, 138 pp.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C***ROUTINES CALLED  GENBUN
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HWSCRT
C
C
      DIMENSION       F(IDIMF,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C***FIRST EXECUTABLE STATEMENT  HWSCRT
      IERROR = 0
      IF (A .GE. B) IERROR = 1
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 2
      IF (C .GE. D) IERROR = 3
      IF (N .LE. 3) IERROR = 4
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 5
      IF (IDIMF .LT. M+1) IERROR = 7
      IF (M .LE. 3) IERROR = 8
      IF (IERROR .NE. 0) RETURN
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND .GT. 0) MPEROD = 1
      DELTAX = (B-A)/M
      TWDELX = 2./DELTAX
      DELXSQ = 1./DELTAX**2
      DELTAY = (D-C)/N
      TWDELY = 2./DELTAY
      DELYSQ = 1./DELTAY**2
      NP = NBDCND+1
      NP1 = N+1
      MP = MBDCND+1
      MP1 = M+1
      NSTART = 1
      NSTOP = N
      NSKIP = 1
      GO TO (104,101,102,103,104),NP
  101 NSTART = 2
      GO TO 104
  102 NSTART = 2
  103 NSTOP = NP1
      NSKIP = 2
  104 NUNK = NSTOP-NSTART+1
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      MSTART = 1
      MSTOP = M
      MSKIP = 1
      GO TO (117,105,106,109,110),MP
  105 MSTART = 2
      GO TO 107
  106 MSTART = 2
      MSTOP = MP1
      MSKIP = 2
  107 DO 108 J=NSTART,NSTOP
         F(2,J) = F(2,J)-F(1,J)*DELXSQ
  108 CONTINUE
      GO TO 112
  109 MSTOP = MP1
      MSKIP = 2
  110 DO 111 J=NSTART,NSTOP
         F(1,J) = F(1,J)+BDA(J)*TWDELX
  111 CONTINUE
  112 GO TO (113,115),MSKIP
  113 DO 114 J=NSTART,NSTOP
         F(M,J) = F(M,J)-F(MP1,J)*DELXSQ
  114 CONTINUE
      GO TO 117
  115 DO 116 J=NSTART,NSTOP
         F(MP1,J) = F(MP1,J)-BDB(J)*TWDELX
  116 CONTINUE
  117 MUNK = MSTOP-MSTART+1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (127,118,118,120,120),NP
  118 DO 119 I=MSTART,MSTOP
         F(I,2) = F(I,2)-F(I,1)*DELYSQ
  119 CONTINUE
      GO TO 122
  120 DO 121 I=MSTART,MSTOP
         F(I,1) = F(I,1)+BDC(I)*TWDELY
  121 CONTINUE
  122 GO TO (123,125),NSKIP
  123 DO 124 I=MSTART,MSTOP
         F(I,N) = F(I,N)-F(I,NP1)*DELYSQ
  124 CONTINUE
      GO TO 127
  125 DO 126 I=MSTART,MSTOP
         F(I,NP1) = F(I,NP1)-BDD(I)*TWDELY
  126 CONTINUE
C
C    MULTIPLY RIGHT SIDE BY DELTAY**2.
C
  127 DELYSQ = DELTAY*DELTAY
      DO 129 I=MSTART,MSTOP
         DO 128 J=NSTART,NSTOP
            F(I,J) = F(I,J)*DELYSQ
  128    CONTINUE
  129 CONTINUE
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2+MUNK
      ID4 = ID3+MUNK
      S = DELYSQ*DELXSQ
      ST2 = 2.*S
      DO 130 I=1,MUNK
         W(I) = S
         J = ID2+I
         W(J) = -ST2+ELMBDA*DELYSQ
         J = ID3+I
         W(J) = S
  130 CONTINUE
      IF (MP .EQ. 1) GO TO 131
      W(1) = 0.
      W(ID4) = 0.
  131 CONTINUE
      GO TO (135,135,132,133,134),MP
  132 W(ID2) = ST2
      GO TO 135
  133 W(ID2) = ST2
  134 W(ID3+1) = ST2
  135 CONTINUE
      PERTRB = 0.
      IF (ELMBDA) 144,137,136
  136 IERROR = 6
      GO TO 144
  137 IF ((NBDCND.EQ.0 .OR. NBDCND.EQ.3) .AND.
     1    (MBDCND.EQ.0 .OR. MBDCND.EQ.3)) GO TO 138
      GO TO 144
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  138 A1 = 1.
      A2 = 1.
      IF (NBDCND .EQ. 3) A2 = 2.
      IF (MBDCND .EQ. 3) A1 = 2.
      S1 = 0.
      MSP1 = MSTART+1
      MSTM1 = MSTOP-1
      NSP1 = NSTART+1
      NSTM1 = NSTOP-1
      DO 140 J=NSP1,NSTM1
         S = 0.
         DO 139 I=MSP1,MSTM1
            S = S+F(I,J)
  139    CONTINUE
         S1 = S1+S*A1+F(MSTART,J)+F(MSTOP,J)
  140 CONTINUE
      S1 = A2*S1
      S = 0.
      DO 141 I=MSP1,MSTM1
         S = S+F(I,NSTART)+F(I,NSTOP)
  141 CONTINUE
      S1 = S1+S*A1+F(MSTART,NSTART)+F(MSTART,NSTOP)+F(MSTOP,NSTART)+
     1     F(MSTOP,NSTOP)
      S = (2.+(NUNK-2)*A2)*(2.+(MUNK-2)*A1)
      PERTRB = S1/S
      DO 143 J=NSTART,NSTOP
         DO 142 I=MSTART,MSTOP
            F(I,J) = F(I,J)-PERTRB
  142    CONTINUE
  143 CONTINUE
      PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
  144 CALL GENBUN (NPEROD,NUNK,MPEROD,MUNK,W(1),W(ID2+1),W(ID3+1),
     1             IDIMF,F(MSTART,NSTART),IERR1,W(ID4+1))
      W(1) = W(ID4+1)+3*MUNK
C
C     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
C
      IF (NBDCND .NE. 0) GO TO 146
      DO 145 I=MSTART,MSTOP
         F(I,NP1) = F(I,1)
  145 CONTINUE
  146 IF (MBDCND .NE. 0) GO TO 148
      DO 147 J=NSTART,NSTOP
         F(MP1,J) = F(1,J)
  147 CONTINUE
      IF (NBDCND .EQ. 0) F(MP1,NP1) = F(1,NP1)
  148 CONTINUE
      RETURN
      END

!============================================================================

c     merge.f

      SUBROUTINE MERGE (TCOS,I1,M1,I2,M2,I3)
      DIMENSION       TCOS(1)
C
C
C     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
C     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
C     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
C     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
C
C
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 .EQ. 0) GO TO 107
      IF (M2 .EQ. 0) GO TO 104
  101 J = J+1
      L = J1+I1
      X = TCOS(L)
      L = J2+I2
      Y = TCOS(L)
      IF (X-Y) 102,102,103
  102 TCOS(J) = X
      J1 = J1+1
      IF (J1 .GT. M1) GO TO 106
      GO TO 101
  103 TCOS(J) = Y
      J2 = J2+1
      IF (J2 .LE. M2) GO TO 101
      IF (J1 .GT. M1) GO TO 109
  104 K = J-J1+1
      DO 105 J=J1,M1
         M = K+J
         L = J+I1
         TCOS(M) = TCOS(L)
  105 CONTINUE
      GO TO 109
  106 CONTINUE
      IF (J2 .GT. M2) GO TO 109
  107 K = J-J2+1
      DO 108 J=J2,M2
         M = K+J
         L = J+I2
         TCOS(M) = TCOS(L)
  108 CONTINUE
  109 CONTINUE
      RETURN
      END

!============================================================================

c     pimach.f

*DECK PIMACH
      FUNCTION PIMACH (DUM)
C***BEGIN PROLOGUE  PIMACH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PIMACH-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subprogram supplies the value of the constant PI correct to
C     machine precision where
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
C***SEE ALSO  HSTCSP, HSTSSP, HWSCSP
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PIMACH
C
C***FIRST EXECUTABLE STATEMENT  PIMACH
      PIMACH = 3.14159265358979
      RETURN
      END

!============================================================================

c     poisd2.f

*DECK POISD2
      SUBROUTINE POISD2 (MR, NR, ISTAG, BA, BB, BC, Q, IDIMQ, B, W, D,
     +   TCOS, P)
C***BEGIN PROLOGUE  POISD2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISD2-S, CMPOSD-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation for Dirichlet boundary
C     conditions.
C
C     ISTAG = 1 if the last diagonal block is the matrix A.
C     ISTAG = 2 if the last diagonal block is the matrix A+I.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  COSGEN, S1MERG, TRIX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine S1MERG rather than deleted
C           routine MERGE.  (WRB)
C***END PROLOGUE  POISD2
C
      DIMENSION       Q(IDIMQ,*) ,BA(*)      ,BB(*)      ,BC(*)      ,
     1                TCOS(*)    ,B(*)       ,D(*)       ,W(*)       ,
     2                P(*)
C***FIRST EXECUTABLE STATEMENT  POISD2
      M = MR
      N = NR
      JSH = 0
      FI = 1./ISTAG
      IP = -M
      IPSTOR = 0
      GO TO (101,102),ISTAG
  101 KR = 0
      IRREG = 1
      IF (N .GT. 1) GO TO 106
      TCOS(1) = 0.
      GO TO 103
  102 KR = 1
      JSTSAV = 1
      IRREG = 2
      IF (N .GT. 1) GO TO 106
      TCOS(1) = -1.
  103 DO 104 I=1,M
         B(I) = Q(I,1)
  104 CONTINUE
      CALL TRIX (1,0,M,BA,BB,BC,B,TCOS,D,W)
      DO 105 I=1,M
         Q(I,1) = B(I)
  105 CONTINUE
      GO TO 183
  106 LR = 0
      DO 107 I=1,M
         P(I) = 0.
  107 CONTINUE
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 L = 2*JST
      NODD = 2-2*((NUN+1)/2)+NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      GO TO (110,109),NODD
  109 JSP = JSP-L
      GO TO 111
  110 JSP = JSP-JST
      IF (IRREG .NE. 1) JSP = JSP-L
  111 CONTINUE
C
C     REGULAR REDUCTION
C
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IF (L .GT. JSP) GO TO 118
      DO 117 J=L,JSP,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         JM3 = JM2-JSH
         JP3 = JP2+JSH
         IF (JST .NE. 1) GO TO 113
         DO 112 I=1,M
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
         GO TO 115
  113    DO 114 I=1,M
            T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
            Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
         CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
         DO 116 I=1,M
            Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
C
C     REDUCTION FOR LAST UNKNOWN
C
  118 GO TO (119,136),NODD
  119 GO TO (152,120),IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (123,121),ISTAG
  121 CONTINUE
      IF (JST .NE. 1) GO TO 123
      DO 122 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = 0.
  122 CONTINUE
      GO TO 130
  123 GO TO (124,126),NODDPR
  124 DO 125 I=1,M
         IP1 = IP+I
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
      GO TO 128
  126 DO 127 I=1,M
         B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
      IP = IP+M
      IPSTOR = MAX(IPSTOR,IP+M)
      DO 131 I=1,M
         IP1 = IP+I
         P(IP1) = Q(I,J)+B(I)
         B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
      IF (LR .NE. 0) GO TO 133
      DO 132 I=1,JST
         KRPI = KR+I
         TCOS(KRPI) = TCOS(I)
  132 CONTINUE
      GO TO 134
  133 CONTINUE
      CALL COSGEN (LR,JSTSAV,0.,FI,TCOS(JST+1))
      CALL S1MERG (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL TRIX (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
      DO 135 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
      LR = KR
      KR = KR+L
      GO TO 152
C
C     EVEN NUMBER OF UNKNOWNS
C
  136 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (137,138),IRREG
  137 CONTINUE
      JSTSAV = JST
      IDEG = JST
      KR = L
      GO TO 139
  138 CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
      KR = KR+JST
  139 IF (JST .NE. 1) GO TO 141
      IRREG = 2
      DO 140 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  140 CONTINUE
      GO TO 150
  141 DO 142 I=1,M
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
      GO TO (143,145),IRREG
  143 DO 144 I=1,M
         Q(I,J) = Q(I,JM2)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
      IRREG = 2
      GO TO 150
  145 CONTINUE
      GO TO (146,148),NODDPR
  146 DO 147 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
      IP = IP-M
      GO TO 150
  148 DO 149 I=1,M
         Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      DO 151 I=1,M
         Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN .GE. 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      DO 153 I=1,M
         B(I) = Q(I,J)
  153 CONTINUE
      GO TO (154,155),IRREG
  154 CONTINUE
      CALL COSGEN (JST,1,0.5,0.0,TCOS)
      IDEG = JST
      GO TO 156
  155 KR = LR+JST
      CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
      IDEG = KR
  156 CONTINUE
      CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      JM1 = J-JSH
      JP1 = J+JSH
      GO TO (157,159),IRREG
  157 DO 158 I=1,M
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
      GO TO 164
  159 GO TO (160,162),NODDPR
  160 DO 161 I=1,M
         IP1 = IP+I
         Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
      IP = IP-M
      GO TO 164
  162 DO 163 I=1,M
         Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN .GT. N) GO TO 183
      DO 182 J=JST,N,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         IF (J .GT. JST) GO TO 166
         DO 165 I=1,M
            B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
         GO TO 170
  166    IF (JP2 .LE. N) GO TO 168
         DO 167 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
         IF (JST .LT. JSTSAV) IRREG = 1
         GO TO (170,171),IRREG
  168    DO 169 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
         CALL COSGEN (JST,1,0.5,0.0,TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    IF (J+L .GT. N) LR = LR-JST
         KR = JST+LR
         CALL COSGEN (KR,JSTSAV,0.0,FI,TCOS)
         CALL COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
         IF (JST .GT. 1) GO TO 174
         DO 173 I=1,M
            Q(I,J) = B(I)
  173    CONTINUE
         GO TO 182
  174    IF (JP2 .GT. N) GO TO 177
  175    DO 176 I=1,M
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
         GO TO 182
  177    GO TO (175,178),IRREG
  178    IF (J+JSH .GT. N) GO TO 180
         DO 179 I=1,M
            IP1 = IP+I
            Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
         IP = IP-M
         GO TO 182
  180    DO 181 I=1,M
            Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
      L = L/2
      GO TO 164
  183 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END

!============================================================================

c     poisn2.f

*DECK POISN2
      SUBROUTINE POISN2 (M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2,
     +   B3, W, W2, W3, D, TCOS, P)
C***BEGIN PROLOGUE  POISN2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISN2-S, CMPOSN-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation with Neumann boundary
C     conditions.
C
C     ISTAG = 1 if the last diagonal block is A.
C     ISTAG = 2 if the last diagonal block is A-I.
C     MIXBND = 1 if have Neumann boundary conditions at both boundaries.
C     MIXBND = 2 if have Neumann boundary conditions at bottom and
C     Dirichlet condition at top.  (for this case, must have ISTAG = 1.)
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine S1MERG rather than deleted
C           routine MERGE.  (WRB)
C***END PROLOGUE  POISN2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
C***FIRST EXECUTABLE STATEMENT  POISN2
      FISTAG = 3-ISTAG
      FNUM = 1./ISTAG
      FDEN = 0.5*(ISTAG-1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103),ISTAG
  101 CONTINUE
      DO 102 I=1,MR
         Q(I,N) = .5*Q(I,N)
  102 CONTINUE
      GO TO (103,104),MIXBND
  103 IF (N .LE. 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      GO TO (105,106),MIXBND
  105 JSTART = 1
      GO TO 107
  106 JSTART = JR
      NROD = 1-NROD
  107 CONTINUE
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      CALL COSGEN (I2R,1,0.5,0.0,TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 108
      J = JR
      GO TO 116
  108 CONTINUE
C
C     REGULAR REDUCTION.
C
      DO 115 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 109
         JM1 = JP1
         JM2 = JP2
         JM3 = JP3
  109    CONTINUE
         IF (I2R .NE. 1) GO TO 111
         IF (J .EQ. 1) JM2 = JP2
         DO 110 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
         GO TO 113
  111    CONTINUE
         DO 112 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 114 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  115 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  116 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 128
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 118
      DO 117 I=1,MR
         B(I) = FISTAG*Q(I,J)
         Q(I,J) = Q(I,JM2)
  117 CONTINUE
      GO TO 126
  118 DO 119 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
      IF (NRODPR .NE. 0) GO TO 121
      DO 120 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
      IP = IP-MR
      GO TO 123
  121 CONTINUE
      DO 122 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 IF (LR .EQ. 0) GO TO 124
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      GO TO 126
  124 CONTINUE
      DO 125 I=1,MR
         B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 127 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
      KR = KR+I2R
      GO TO 151
  128 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 135
      DO 129 I=1,MR
         B(I) = Q(I,J)
  129 CONTINUE
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      GO TO (133,130),ISTAG
  130 DO 131 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  131 CONTINUE
      TCOS(1) = 1.
      TCOS(2) = 0.
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 132 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
      GO TO 150
  133 CONTINUE
      DO 134 I=1,MR
         P(I) = B(I)
         Q(I,J) = Q(I,JM2)+2.*Q(I,JP2)+3.*B(I)
  134 CONTINUE
      GO TO 150
  135 CONTINUE
      DO 136 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
      IF (NRODPR .NE. 0) GO TO 138
      DO 137 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  137 CONTINUE
      GO TO 140
  138 CONTINUE
      DO 139 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX(IPSTOR,IP+MR)
      DO 141 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
      IF (LR .EQ. 0) GO TO 142
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(I2R+1))
      CALL S1MERG (TCOS,0,I2R,I2R,LR,KR)
      GO TO 144
  142 DO 143 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  143 CONTINUE
  144 CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      IF (LR .NE. 0) GO TO 145
      GO TO (146,145),ISTAG
  145 CONTINUE
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      GO TO 148
  146 CONTINUE
      DO 147 I=1,MR
         B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
      DO 149 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
      LR = KR
      KR = KR+JR
  151 CONTINUE
      GO TO (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 155
      GO TO 154
  153 NR = NLAST/JR
      IF (NR .LE. 1) GO TO 192
  154 I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 184
      IF (LR .NE. 0) GO TO 170
      IF (N .NE. 3) GO TO 161
C
C     CASE N = 3.
C
      GO TO (156,168),ISTAG
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,2)
  157 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 158 I=1,MR
         Q(I,2) = B(I)
         B(I) = 4.*B(I)+Q(I,1)+2.*Q(I,3)
  158 CONTINUE
      TCOS(1) = -2.
      TCOS(2) = 2.
      I1 = 2
      I2 = 0
      CALL TRIX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
      DO 159 I=1,MR
         Q(I,2) = Q(I,2)+B(I)
         B(I) = Q(I,1)+2.*Q(I,2)
  159 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 160 I=1,MR
         Q(I,1) = B(I)
  160 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
C
C     CASE N = 2**P+1
C
  161 CONTINUE
      GO TO (162,170),ISTAG
  162 CONTINUE
      DO 163 I=1,MR
         B(I) = Q(I,J)+.5*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 164 I=1,MR
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
         B(I) = Q(I,1)+2.*Q(I,NLAST)+4.*Q(I,J)
  164 CONTINUE
      JR2 = 2*JR
      CALL COSGEN (JR,1,0.0,0.0,TCOS)
      DO 165 I=1,JR
         I1 = JR+I
         I2 = JR+1-I
         TCOS(I1) = -TCOS(I2)
  165 CONTINUE
      CALL TRIX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  166 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 167 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
      GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168 DO 169 I=1,MR
         B(I) = Q(I,2)
         Q(I,2) = 0.
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  169 CONTINUE
      JR = 1
      I2R = 0
      J = 2
      GO TO 177
  170 CONTINUE
      DO 171 I=1,MR
         B(I) = .5*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
      IF (NROD .NE. 0) GO TO 173
      DO 172 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  172 CONTINUE
      GO TO 175
  173 DO 174 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
      DO 176 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+2.*T
  176 CONTINUE
  177 CONTINUE
      K1 = KR+2*JR-1
      K2 = KR+JR
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (K2+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+K2+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL S1MERG (TCOS,K1,K2,K1+K2,JR-1,0)
      K3 = K1+K2+LR
      CALL COSGEN (JR,1,0.5,0.0,TCOS(K3+1))
      K4 = K3+JR+1
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
      CALL S1MERG (TCOS,K3,JR,K3+JR,KR,K1)
      IF (LR .EQ. 0) GO TO 178
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      CALL S1MERG (TCOS,K3,JR,K3+JR,LR,K3-LR)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K4))
  178 K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 179 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 180 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  180 CONTINUE
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 182
      DO 181 I=1,MR
         Q(I,1) = B(I)
  181 CONTINUE
      GO TO 194
  182 CONTINUE
      DO 183 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
      GO TO 194
  184 CONTINUE
      IF (N .NE. 2) GO TO 188
C
C     CASE  N = 2
C
      DO 185 I=1,MR
         B(I) = Q(I,1)
  185 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 186 I=1,MR
         Q(I,1) = B(I)
         B(I) = 2.*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
      TCOS(1) = -FISTAG
      TCOS(2) = 2.
      CALL TRIX (2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
  188 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 189 I=1,MR
         II = IP+I
         B3(I) = 0.
         B(I) = Q(I,1)+2.*P(II)
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)
         B2(I) = 2.*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
      K1 = KR+JR-1
      TCOS(K1+1) = -2.
      K4 = K1+3-ISTAG
      CALL COSGEN (KR+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+KR+1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
      CALL S1MERG (TCOS,K1,KR,K1+KR,JR-1,0)
      CALL COSGEN (KR,1,0.5,FDEN,TCOS(K1+1))
      K2 = KR
      K4 = K1+K2+1
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 190 I=1,MR
         B(I) = B(I)+B2(I)
  190 CONTINUE
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 191 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
      GO TO 194
  192 DO 193 I=1,MR
         B(I) = Q(I,NLAST)
  193 CONTINUE
      GO TO 196
  194 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 195 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 198
      DO 197 I=1,MR
         Q(I,NLAST) = 0.
  197 CONTINUE
      GO TO 202
  198 CONTINUE
      IF (NROD .NE. 0) GO TO 200
      DO 199 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  199 CONTINUE
      IP = IP-MR
      GO TO 202
  200 DO 201 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
      CALL COSGEN (KR,1,0.5,FDEN,TCOS)
      CALL COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
      IF (LR .NE. 0) GO TO 204
      DO 203 I=1,MR
         B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 205 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 222
      GO TO (207,208),MIXBND
  207 JSTART = 1+JR
      GO TO 209
  208 JSTART = JR
  209 CONTINUE
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 210
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 211
  210 CONTINUE
      JSTOP = NLAST-JR
  211 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      DO 221 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 213
         DO 212 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
         GO TO 215
  213    CONTINUE
         DO 214 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
         IF (JR .NE. 1) GO TO 217
         DO 216 I=1,MR
            Q(I,J) = 0.
  216    CONTINUE
         GO TO 219
  217    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 218 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 220 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END

!============================================================================

c     poisp2.f

*DECK POISP2
      SUBROUTINE POISP2 (M, N, A, BB, C, Q, IDIMQ, B, B2, B3, W, W2, W3,
     +   D, TCOS, P)
C***BEGIN PROLOGUE  POISP2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISP2-S, CMPOSP-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson equation with periodic boundary
C     conditions.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  POISD2, POISN2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  POISP2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                P(*)
C***FIRST EXECUTABLE STATEMENT  POISP2
      MR = M
      NR = (N+1)/2
      NRM1 = NR-1
      IF (2*NR .NE. N) GO TO 107
C
C     EVEN NUMBER OF UNKNOWNS
C
      DO 102 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 101 I=1,MR
            S = Q(I,NRMJ)-Q(I,NRPJ)
            T = Q(I,NRMJ)+Q(I,NRPJ)
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  101    CONTINUE
  102 CONTINUE
      DO 103 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
         Q(I,N) = 2.*Q(I,N)
  103 CONTINUE
      CALL POISD2 (MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 105 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 104 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,NRMJ))
            T = .5*(Q(I,NRPJ)-Q(I,NRMJ))
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
         Q(I,N) = .5*Q(I,N)
  106 CONTINUE
      GO TO 118
  107 CONTINUE
C
C     ODD  NUMBER OF UNKNOWNS
C
      DO 109 J=1,NRM1
         NRPJ = N+1-J
         DO 108 I=1,MR
            S = Q(I,J)-Q(I,NRPJ)
            T = Q(I,J)+Q(I,NRPJ)
            Q(I,J) = S
            Q(I,NRPJ) = T
  108    CONTINUE
  109 CONTINUE
      DO 110 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
  110 CONTINUE
      LH = NRM1/2
      DO 112 J=1,LH
         NRMJ = NR-J
         DO 111 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  111    CONTINUE
  112 CONTINUE
      CALL POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 114 J=1,NRM1
         NRPJ = NR+J
         DO 113 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,J))
            T = .5*(Q(I,NRPJ)-Q(I,J))
            Q(I,NRPJ) = T
            Q(I,J) = S
  113    CONTINUE
  114 CONTINUE
      DO 115 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
  115 CONTINUE
      DO 117 J=1,LH
         NRMJ = NR-J
         DO 116 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  116    CONTINUE
  117 CONTINUE
  118 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END

!============================================================================

c     s1merg.f

*DECK S1MERG
      SUBROUTINE S1MERG (TCOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  S1MERG
C***SUBSIDIARY
C***PURPOSE  Merge two strings of ascending real numbers.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   This subroutine merges two ascending strings of numbers in the
C   array TCOS.  The first string is of length M1 and starts at
C   TCOS(I1+1).  The second string is of length M2 and starts at
C   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  SCOPY
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did
C           not compile correctly with optimization on the IBM RS6000.
C           (RWC)
C   920130  Code name changed from MERGE to S1MERG.  (WRB)
C***END PROLOGUE  S1MERG
      INTEGER I1, I2, I3, M1, M2
      REAL TCOS(*)
C
      INTEGER J1, J2, J3
C
C***FIRST EXECUTABLE STATEMENT  S1MERG
      IF (M1.EQ.0 .AND. M2.EQ.0) RETURN
C
      IF (M1.EQ.0 .AND. M2.NE.0) THEN
         CALL SCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      IF (M1.NE.0 .AND. M2.EQ.0) THEN
         CALL SCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      J1 = 1
      J2 = 1
      J3 = 1
C
   10 IF (TCOS(I1+J1) .LE. TCOS(I2+J2)) THEN
         TCOS(I3+J3) = TCOS(I1+J1)
         J1 = J1+1
         IF (J1 .GT. M1) THEN
            CALL SCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ELSE
         TCOS(I3+J3) = TCOS(I2+J2)
         J2 = J2+1
         IF (J2 .GT. M2) THEN
            CALL SCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ENDIF
      J3 = J3+1
      GO TO 10
      END

!============================================================================

c     scopy.f

*DECK SCOPY
      SUBROUTINE SCOPY (N, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      SINGLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  copy of vector SX (unchanged if N .LE. 0)
C
C     Copy single precision SX to single precision SY.
C     For I = 0 to N-1, copy  SX(LX+I*INCX) to SY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SCOPY
      REAL SX(*), SY(*)
C***FIRST EXECUTABLE STATEMENT  SCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I+1) = SX(I+1)
        SY(I+2) = SX(I+2)
        SY(I+3) = SX(I+3)
        SY(I+4) = SX(I+4)
        SY(I+5) = SX(I+5)
        SY(I+6) = SX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SY(I) = SX(I)
   70 CONTINUE
      RETURN
      END

!============================================================================

c     tri3.f

*DECK TRI3
      SUBROUTINE TRI3 (M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
C***BEGIN PROLOGUE  TRI3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRI3-S, CMPTR3-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve three linear systems whose common coefficient
C     matrix is a rational function in the matrix given by
C
C                  TRIDIAGONAL (...,A(I),B(I),C(I),...)
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRI3
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,K(4)       ,
     1                TCOS(*)    ,Y1(*)      ,Y2(*)      ,Y3(*)      ,
     2                D(*)       ,W1(*)      ,W2(*)      ,W3(*)
      INTEGER K1P1, K2P1, K3P1, K4P1
C
C***FIRST EXECUTABLE STATEMENT  TRI3
      MM1 = M-1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      K1P1 = K1+1
      K2P1 = K2+1
      K3P1 = K3+1
      K4P1 = K4+1
      K2K3K4 = K2+K3+K4
      IF (K2K3K4 .EQ. 0) GO TO 101
      L1 = (K1+1)/(K2+1)
      L2 = (K1+1)/(K3+1)
      L3 = (K1+1)/(K4+1)
      LINT1 = 1
      LINT2 = 1
      LINT3 = 1
      KINT1 = K1
      KINT2 = KINT1+K2
      KINT3 = KINT2+K3
  101 CONTINUE
      DO 115 N=1,K1
         X = TCOS(N)
         IF (K2K3K4 .EQ. 0) GO TO 107
         IF (N .NE. L1) GO TO 103
         DO 102 I=1,M
            W1(I) = Y1(I)
  102    CONTINUE
  103    IF (N .NE. L2) GO TO 105
         DO 104 I=1,M
            W2(I) = Y2(I)
  104    CONTINUE
  105    IF (N .NE. L3) GO TO 107
         DO 106 I=1,M
            W3(I) = Y3(I)
  106    CONTINUE
  107    CONTINUE
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO 108 I=2,M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
  108    CONTINUE
         DO 109 IP=1,MM1
            I = M-IP
            Y1(I) = Y1(I)-D(I)*Y1(I+1)
            Y2(I) = Y2(I)-D(I)*Y2(I+1)
            Y3(I) = Y3(I)-D(I)*Y3(I+1)
  109    CONTINUE
         IF (K2K3K4 .EQ. 0) GO TO 115
         IF (N .NE. L1) GO TO 111
         I = LINT1+KINT1
         XX = X-TCOS(I)
         DO 110 I=1,M
            Y1(I) = XX*Y1(I)+W1(I)
  110    CONTINUE
         LINT1 = LINT1+1
         L1 = (LINT1*K1P1)/K2P1
  111    IF (N .NE. L2) GO TO 113
         I = LINT2+KINT2
         XX = X-TCOS(I)
         DO 112 I=1,M
            Y2(I) = XX*Y2(I)+W2(I)
  112    CONTINUE
         LINT2 = LINT2+1
         L2 = (LINT2*K1P1)/K3P1
  113    IF (N .NE. L3) GO TO 115
         I = LINT3+KINT3
         XX = X-TCOS(I)
         DO 114 I=1,M
            Y3(I) = XX*Y3(I)+W3(I)
  114    CONTINUE
         LINT3 = LINT3+1
         L3 = (LINT3*K1P1)/K4P1
  115 CONTINUE
      RETURN
      END

!============================================================================

c     trix.f

*DECK TRIX
      SUBROUTINE TRIX (IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
C***BEGIN PROLOGUE  TRIX
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve a system of linear equations where the
C     coefficient matrix is a rational function in the matrix given by
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRIX
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,
     1                TCOS(*)    ,D(*)       ,W(*)
      INTEGER KB, KC
C***FIRST EXECUTABLE STATEMENT  TRIX
      MM1 = M-1
      KB = IDEGBR+1
      KC = IDEGCR+1
      L = (IDEGBR+1)/(IDEGCR+1)
      LINT = 1
      DO 108 K=1,IDEGBR
         X = TCOS(K)
         IF (K .NE. L) GO TO 102
         I = IDEGBR+LINT
         XX = X-TCOS(I)
         DO 101 I=1,M
            W(I) = Y(I)
            Y(I) = XX*Y(I)
  101    CONTINUE
  102    CONTINUE
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO 103 I=2,MM1
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  103    CONTINUE
         Z = B(M)-X-A(M)*D(MM1)
         IF (Z .NE. 0.) GO TO 104
         Y(M) = 0.
         GO TO 105
  104    Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  105    CONTINUE
         DO 106 IP=1,MM1
            I = M-IP
            Y(I) = Y(I)-D(I)*Y(I+1)
  106    CONTINUE
         IF (K .NE. L) GO TO 108
         DO 107 I=1,M
            Y(I) = Y(I)+W(I)
  107    CONTINUE
         LINT = LINT+1
         L = (LINT*KB)/KC
  108 CONTINUE
      RETURN
      END

!============================================================================
