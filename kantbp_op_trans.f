      SUBROUTINE KANTBP(TITLE,IPTYPE,ISC,ISCAT,NROOT,MDIM,IDIM,NPOL,
     1                  RTOL,NITEM,SHIFT,IPRINT,IPRSTP,NMESH,RMESH,
     2                  NDIR,NDIL,NMDIL,THRSHL,THRSHR,IBOUND,FNOUT,
     3                  IOUT,POTEN,IOUP,FMATR,IOUM,EVWFN,IOUF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                                                                      C
C             #     #    #    #     # ####### ######  ######           C
C             #    #    # #   ##    #    #    #     # #     #          C
C             #   #    #   #  # #   #    #    #     # #     #          C
C             ####    ####### #  #  #    #    ######  ######           C
C             #   #   #     # #   # #    #    #     # #                C
C             #    #  #     # #    ##    #    #     # #                C
C             #     # #     # #     #    #    ######  #                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                                                                      C
C                    version  3.1  from  28.12.2021                    C
C                                                                      C
C                                                                      C
C                              written  by                             C
C                                                                      C
C                                                                      C
C                O. CHULUUNBAATAR, A.A. GUSEV, S.I. VINITSKY           C
C                                                                      C
C                   JOINT INSTITUTE FOR NUCLEAR RESEARCH,              C
C                    DUBNA, 141980 MOSCOW REGION, RUSSIA               C
C                                                                      C
C                           A.G. ABRASHKEVICH                          C
C                                                                      C
C                   IBM TORONTO LAB, 8200 WARDEN AVENUE,               C
C                       MARKHAM, ON L6G 1C7, CANADA                    C
C                                                                      C
C                          P.W. Wen, C.J. Lin                          C
C                                                                      C
C                  China Institute of Atomic Energy,                   C
C                         102413 Beijing, China                        C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
C  Program KANTBP calculates energy values, reflection and             C
C  transmission matrices, and corresponding wave functions in a        C
C  coupled-channel approximation of the adiabatic approach. In this    C
C  approach, a multi-dimensional Schrodinger equation is reduced to    C
C  a system of the coupled second-order ordinary differential          C
C  equations on the finite interval with homogeneous boundary          C
C  conditions of the third type. The resulting system of equations     C
C  which contains the potential matrix elements and first-derivative   C
C  coupling terms is solved using high-order accuracy approximations   C
C  of the finite element method.                                       C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
C                        I N P U T  D A T A                            C
C                        ------------------                            C
C                                                                      C
C TITLE   CHARACTER  title of the run to be printed on the output      C
C                    listing. The title should be no longer than       C
C                    70 characters.                                    C
C IPTYPE  INTEGER    IPTYPE contains information about type of the     C
C                    problem solved. If IPTYPE = 0 the program         C
C                    calculates the eigenvalue problem; otherwise,     C
C                    it calculates the scattering problem.             C
C ISC     INTEGER    flag for performing the calculation of the        C
C                    reflection and transmission matrices:             C
C                    = 1 -- if NOPENR > 0, the calculation of the      C
C                    reflection matrix is carried out only with the    C
C                    direction <--. If NOPENL > 0, the calculation of  C
C                    the transmission matrix is carried out            C
C                    additionally. The properties of the reflection    C
C                    and transmission matrices are also verified.      C
C                    If NOPENR = 0, the message NO RIGHT OPEN CHANNELS C
C                    is printed;                                       C
C                    = 2 -- if NOPENL > 0, the calculation of the      C
C                    reflection matrix is carried out only with the    C
C                    direction -->. If NOPENR > 0, the calculation of  C
C                    the transmission matrix is carried out            C
C                    additionally. The properties of the reflection    C
C                    and transmission matrices are also verified.      C
C                    If NOPENL = 0, the message NO LEFT OPEN CHANNELS  C
C                    is printed;                                       C
C                    = 3 -- if NOPENR > 0 and NOPENL > 0, the          C
C                    calculation of the reflection and transmission    C
C                    matrices is carried out with both directions <--  C
C                    and -->. The properties of the reflection and     C
C                    transmission matrices are also verified.          C
C                    Here NOPENL and NOPENR are the number of the left C
C                    and right open channels.                          C
C ISCAT   INTEGER    flag for performing the calculation of the left   C
C                    THRSHL and right THRSHR threshold values:         C
C                    = 1 -- V(z) and Q(z) matrices have the special    C
C                    asymptotic behaviour at Z <= Z_MIN and            C
C                    Z >= Z_MAX. The user must set both the left       C
C                    THRSHL and right THRSHR threshold values;         C
C                    = 2 -- V(z) and Q(z) matrices have the special    C
C                    asymptotic behaviour at Z <= Z_MIN, and Q(z) = 0  C
C                    and V(z) is the constant matrix at Z >= Z_MAX.    C
C                    The user must set the left THRSHL threshold       C
C                    values, and the right THRSHR values are           C
C                    calculated by the program;                        C
C                    = 3 --  Q(z) = 0 and V(z) is the constant matrix  C
C                    at Z <= Z_MIN, the V(z) and Q(z) matrices have    C
C                    the special asymptotic behaviour at Z >= Z_MAX.   C
C                    The user must set the right THRSHR threshold      C
C                    values, and the left THRSHL threshold values are  C
C                    calculated by the program;                        C
C                    = 4 -- Q(z)=0 and V(z) is the constant matrix at  C
C                    Z <= Z_MIN and Z >= Z_MAX. The user does not set  C
C                    the left THRSHL and right THRSHR threshold        C
C                    values, they are calculated by the program.       C
C NROOT   INTEGER    number of eigenvalues (energy levels) and         C
C                    eigenvectors (radial wave functions) required.    C
C                    NROOT should be equals to 1 in case of IBOUND > 4 C
C                    and not used for the scattering problem.          C
C MDIM    INTEGER    maximum number of coupled differential equations. C
C IDIM    INTEGER    dimension of the envelope space and always should C
C                    be equal to 1 for the scattering problem.         C
C NPOL    INTEGER    order of finite-element shape functions           C
C                    (interpolating Lagrange polynomials).             C
C                    Usually set to 6.                                 C
C RTOL    REAL(8)    convergence tolerance on eigenvalues (1.D-06 or   C
C                    smaller). This value is not used for the          C
C                    scattering problem.                               C
C NITEM   INTEGER    maximum number of iterations permitted (usually   C
C                    set to 16). This value is not used for the        C
C                    scattering problem.                               C
C SHIFT   REAL(8)    For the eigenvalue problem, SHIFT contains the    C
C                    energy spectrum. If SHIFT = 0 the value of the    C
C                    energy shift is determined automatically by the   C
C                    program; otherwise, the NROOT eigenvalues and     C
C                    eigenvectors closest to the shift given are       C
C                    calculated (the nonzero value of SHIFT is         C
C                    recommended since it significantly speeds up the  C
C                    calculation). For the scattering problem, SHIFT   C
C                    contains the given double energy spectrum.        C
C IPRINT  INTEGER    level of print:                                   C
C                    = 0 -- minimal level of print. The initial data,  C
C                    short information about the numerical scheme      C
C                    parameters, main flags and keys, and energy       C
C                    values calculated are printed out;                C
C                    = 1 -- radial functions calculated are printed    C
C                    out with step IPRSTP additionally;                C
C                    = 2 -- potential matrix is printed out with step  C
C                    IPRSTP;                                           C
C                    = 3 -- information about nodal point distribution C
C                    is printed out;                                   C
C                    = 4 -- global matrices A and B are printed out    C
C                    additionally;                                     C
C                    = 5 -- the highest level of print. The local      C
C                    stiffness and mass matrices together with all     C
C                    current information about the course of the       C
C                    subspace iteration method solution of the         C
C                    generalized eigenvalue problem are printed out.   C
C IPRSTP  INTEGER    step with which potential matrix and radial wave  C
C                    functions are printed out.                        C
C NMESH   INTEGER    dimension of array RMESH. NMESH always should be  C
C                    odd and > = 3.                                    C
C RMESH   REAL(8)     array RMESH contains information about           C
C                    subdivision of interval [Z_MIN,Z_MAX] of          C
C                    hyperradius RHO on subintervals. The whole        C
C                    interval [Z_MIN,Z_MAX] is divided as follows:     C
C                    RMESH(1) = Z_MIN, RMESH(NMESH) = Z_MAX, and the   C
C                    values of RMESH(I) set the number of elements     C
C                    for each subinterval [RMESH(I-1), RMESH(I+1)],    C
C                    where I=2, 4, ... , NMESH-1.                      C
C NDIR    INTEGER    dimension of array NDIL. If NDIR > MDIM the       C
C                    message about the error is printed and the        C
C                    execution of the program is stopped.              C
C NDIL    INTEGER    array NDIL containing information about the set   C
C                    of numbers of coupled differential equations and  C
C                    always should be  NDIL(NDIR) = < MDIM.            C
C NMDIL   INTEGER    key parameter. If NMDIL = 0 the potential matrix  C
C                    elements of radial coupling are calculated and    C
C                    written to file POTEN; otherwise, they are read   C
C                    from file POTEN.                                  C
C THRSHL  REAL(8)    array THRSHL of dimension MDIM containing values  C
C                    of the left thresholds. This array is not used    C
C                    for the eigenvalue problem.                       C
C THRSHR  REAL(8)    array THRSHL of dimension MDIM containing values  C
C                    of the right thresholds. This array is not used   C
C                    for the eigenvalue problem.                       C
C IBOUND  INTEGER    parameter defining the type of boundary           C
C                    conditions set in the boundary points Z = Z_MIN   C
C                    conditions set in the boundary points Z = Z_MIN   C
C                    and Z = Z_MAX:                                    C
C                    = 1 -- the Dirichlet -- Dirichlet boundary        C
C                    conditions;                                       C
C                    = 2 -- the Dirichlet -- Neumann boundary          C
C                    conditions;                                       C
C                    = 3 -- the Neumann -- Dirichlet boundary          C
C                    conditions;                                       C
C                    =  4 -- the Neumann -- Neumann boundary           C
C                    conditions;                                       C
C                    = 6 -- the Dirichlet -- third type boundary       C
C                    conditions (only for NROOT = 1):                  C
C                    at Z = Z_MIN the Dirichlet boundary condition     C
C                    (see the case 2) is used and at Z_MAX the         C
C                    user-supplied subroutine ASYMEV for the           C
C                    calculation of initial value LAMBDA^(0)(Z_MAX);   C
C                    = 8 -- the Neumann -- third type boundary         C
C                    conditions (only for NROOT = 1):                  C
C                    at Z = Z_MIN the Neumann boundary condition (see  C
C                    the case 4) is used and at Z_MAX the same         C
C                    boundary condition is used as in case 6.          C
C                    IBOUND always should be equal to 8 for the        C
C                    scattering problem. In this case user-supplied    C
C                    subroutines ASYMSL AND ASYMSR for the calculation C
C                    of the regular, irregular asymptotic              C
C                    matrix-solutions and their derivatives            C
C                    at Z = Z_MIN and Z = Z_MAX are used.              C
C FNOUT   CHARACTER  name of the output file (up to 55 characters)     C
C                    for printing out the results of the calculation.  C
C                    It is system specific and may include a complete  C
C                    path to the file location.                        C
C IOUT    INTEGER    number of the output logical device for printing  C
C                    out the results of the calculation (usually set   C
C                    to 7).                                            C
C POTEN   CHARACTER  name of the input/output file (up to 55           C
C                    characters) containing potential matrix elements  C
C                    of radial coupling.                               C
C IOUP    INTEGER    number of the logical device for reading/storing  C
C                    data from/into file POTEN.                        C
C FMATR   CHARACTER  name of the scratch file (up to 55 characters)    C
C                    for storing stiffness matrix.                     C
C IOUM    INTEGER    number of the logical device for storing          C
C                    stiffness matrix.                                 C
C EVWFN   CHARACTER  name of the output file (up to 55 characters) for C
C                    storing the results of the calculation, namely,   C
C                    the energy values or reaction matrix,             C
C                    finite-element grid points, and radial wave       C
C                    functions. It is used only if IOUF > 0.           C
C IOUF    INTEGER    number of the logical device for storing data     C
C                    into file EVWFN.                                  C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
C                        O U T P U T  D A T A                          C
C                        --------------------                          C
C                                                                      C
C The results of the calculation of energy values or reaction matrix   C
C and radial wave functions are written using unformatted segmented    C
C records into file EVWFN, according to the following operator:        C
C                                                                      C
C     WRITE(IOUF) NDIM,NN,NROOT,NGRID,(EIGV(I),I=1,NROOT),             C
C    1            (XGRID(I),I=1,NGRID),((R(I,J),I=1,NN),J=1,NROOT)     C
C                                                                      C
C or                                                                   C
C                                                                      C
C      WRITE(IOUF) NDIM,NN,NL,NR,NGRID,((RR(I,J),I=1,NR),J=1,NR),      C
C     1                                ((TT(I,J),I=1,NL),J=1,NR),      C
C     1            (XGRID(I),I=1,NGRID),((R(I,J),I=1,NN),J=1,NR)       C
C                                                                      C
C In the above, parameters presented in the WRITE statement have the   C
C following meaning:                                                   C
C NDIM is the number of equations.                                     C
C NGRID is the number of finite-element grid points.                   C
C NN = NDIM * NGRID, if IBOIND = 4 or 8;                               C
C NN = NDIM * (NGRID - 1), if IBOIND = 3 or 2 or 6;                    C
C NN = NDIM * (NGRID - 2), if IBOIND = 1.                              C
C NROOT is the number of roots (energy levels).                        C
C Array EIGV contains the energy values calculated.                    C
C NL and NR are the numbers of open channels: NL = NOPENL, NR = NOPENR C
C       for <--, NL = NOPENR, NR = NOPENL for -->.                     C
C Arrays RR and TT contain the reflection and transmission matrices    C
C       values calculated, respectively.                               C
C Array XGRID contains the values of the finite-element grid points.   C
C Array R contains NROOT or NR eigenfunctions each per NN elements     C
C       in length stored in the following way: for each of the NGRID   C
C       mesh points per NDIM elements of eigenfunction                 C
C       (see scheme below):                                            C
C                                                                      C
C        1-st                  2-nd        . . .       last            C
C                                                                      C
C                1                     1                      1        C
C                2                     2                      2        C
C    1-st point  .         1-st point  .   . . .   1-st point .        C
C                .                     .                      .        C
C                .                     .                      .        C
C             NDIM                  NDIM                   NDIM        C
C                                                                      C
C                1                     1                      1        C
C                2                     2                      2        C
C    2-nd point  .         2-nd point  .   . . .   2-nd point .        C
C                .                     .                      .        C
C                .                     .                      .        C
C             NDIM                  NDIM                   NDIM        C
C                                                                      C
C                .                     .                      .        C
C                .                     .                      .        C
C                .                     .                      .        C
C                                                                      C
C                1                     1                      1        C
C                2                     2                      2        C
C    last point  .         last point  .   . . .   last point .        C
C                .                     .                      .        C
C                .                     .                      .        C
C             NDIM                  NDIM                   NDIM        C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C                                                                      C
C                         PUBLISHED PROGRAMS USED                      C
C                         -----------------------                      C
C                                                                      C
C          Program SSPACE to solve the generalized eigenvalue          C
C          problem by the subspace iteration method. Published in:     C
C          K.- J. Bathe, Finite ELement Procedures in Engineering      C
C          Analysis (Prentice-Hall, Englewood Cliffs, NY, 1982)        C
C                                                                      C
C          Program GAULEG to calculate nodes and weights of the        C
C          Gauss-Legendre quadrature. Published in: W.H. Press,        C
C          B.F. Flanery, S.A. Teukolsky, and W.T. Vetterley,           C
C          Numerical Recipes: The Art of Scientific Computing          C
C          (Cambridge University Press, Cambridge, 1986)               C
C                                                                      C
C          Program GAUSSJ to calculate solution of the linear          C
C          equation by Gauss-Jordan elimination. Published in:         C
C          W.H. Press, B.F. Flanery, S.A. Teukolsky, and               C
C          W.T. Vetterley, Numerical Recipes: The Art of Scientific    C
C          Computing (Cambridge University Press, Cambridge, 1986)     C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      IMPLICIT REAL(8) (A-H,O-Y)
      IMPLICIT COMPLEX(8) (Z)
      DIMENSION RMESH(NMESH),NDIL(NDIR),THRSHL(MDIM),THRSHR(MDIM)
      REAL(8),     DIMENSION(:), ALLOCATABLE :: XGT
      REAL(8),     DIMENSION(:), ALLOCATABLE :: WGT
      REAL(8),     DIMENSION(:), ALLOCATABLE :: GRT
      REAL(8),     DIMENSION(:), ALLOCATABLE :: TOT
      COMPLEX(8),  DIMENSION(:), ALLOCATABLE :: ZTOT,ZTOTM
      INTEGER,     DIMENSION(:), ALLOCATABLE :: ITOT
      INTEGER ISTAT
      CHARACTER*70 TITLE,PROBTYPE
      CHARACTER*55 FNOUT,POTEN,FMATR,EVWFN
      DATA DZERO / 0.D0 /, ZERO / 0.D0 /, ONE / 1.D0 /
C
C      *******************************************
C      ***                                     ***
C      ***    I N P U T  D A T A  P H A S E    ***
C      ***                                     ***
C      *******************************************
C
C      THE FOLLOWING FILES ARE USED
C
C          IOUT   =  OUTPUT LOGICAL DEVICE NUMBER FOR
C                    PRINTING OUT RESULTS OF THE CALCULATION
C          IOUP   =  INPUT DEVICE WITH POTENTIAL MATRIX ELEMENTS
C          IOUM   =  OUTPUT STIFFNESS MATRIX
C          IOUF   =  OUTPUT LOGICAL DEVICE FOR STORING
C                    RADIAL WAVE FUNCTIONS
C
      OPEN (UNIT=IOUT,STATUS='UNKNOWN',FILE=FNOUT)
      OPEN (UNIT=IOUP,STATUS='UNKNOWN',FILE=POTEN,FORM='UNFORMATTED')
      OPEN (UNIT=IOUM,STATUS='UNKNOWN',FILE=FMATR,FORM='UNFORMATTED')
      IF (IOUF .GT. 0)
     1OPEN (UNIT=IOUF,STATUS='UNKNOWN',FILE=EVWFN,FORM='UNFORMATTED')
C
      IF (IPTYPE .EQ. 0) THEN
        PROBTYPE = 'EIGEN'
      ELSE
        PROBTYPE = 'SCATTER'
      END IF
C
C      THE RELATIVE ERROR OF THE CALCULATION
C
      A = 4.D0/3.0D0
   11 B = A - ONE
      C = B + B + B
      EPSY = DABS(C - ONE)
cx      IF (EPSY .EQ. DZERO) GO TO 11
      EPSY = 2.4D-16
C
C      FORM AND CHECK INPUT DATA
C
      CALL CHECKN(PROBTYPE,SHIFT,THRSHL,RTOL,EPSY,MDIM,NDIR,NDIL,
     1            IDIM,NROOT,IBOUND,IOUT,IOUP,IOUM)
C
      IF (PROBTYPE .EQ. 'EIGEN') THEN
        IF (NROOT .LT. 25) THEN
          NC = MIN0(2 * NROOT,NROOT + 8)
        ELSE
          NC = NROOT + 5
        END IF
        NNC  = NC * (NC + 1) / 2
      END IF
      NNODE  = NPOL + 1
      NGQ    = NNODE
C
C      SET GAUSSIAN QUADRATURE NODES AND WEIGHTS
C
      IXG = 1
      IWG = IXG + NGQ
C
C          Allocate space
C
      ALLOCATE (XGT(NGQ), STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2010) 'XGT', 'GAULEG', ISTAT
          STOP
      ENDIF
      ALLOCATE (WGT(NGQ), STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2010) 'WGT', 'GAULEG', ISTAT
          DEALLOCATE (XGT)
          STOP
      ENDIF
C
      CALL GAULEG(-ONE,ONE,XGT,WGT,NGQ,EPSY)
C
C      *****************************************************
C      ***                                               ***
C      ***   N O D A L  G E N E R A T I N G  P H A S E   ***
C      ***                                               ***
C      *****************************************************
C
      IGRID = IWG + NGQ
      NELEM = 0
      NG    = 0
      NL    = 0
      DO 25 I = 3 , NMESH, 2
          XM  = RMESH(I) - RMESH(I-2)
          XS  = XM / IDINT(RMESH(I-1)) / DFLOAT(NPOL)
          NM  = IDINT(XM / XS + 0.000000000001)
          NL  = NL + IDINT(RMESH(I-1))
          DO 15 J = 1 , NM
              NG = NG + 1
   15     CONTINUE
   25 CONTINUE
      NELEM = NL
      NGRID = NG + 1
C
C          Allocate space
C
      ALLOCATE (GRT(NGRID), STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2010) 'GRT', 'FEGRID', ISTAT
          DEALLOCATE (XGT)
          DEALLOCATE (WGT)
          STOP
      ENDIF
C
      CALL FEGRID(RMESH,NMESH,GRT,NGRID,NELEM,NPOL)
C
      IF (IPRINT .GE. 0) THEN
        WRITE(IOUT,1000) TITLE
        IF (PROBTYPE .EQ. 'EIGEN') THEN
          WRITE(IOUT,1010) MDIM,NROOT,NELEM,NGRID,NPOL,NGQ,NC,
     1                     IDIM,IBOUND,SHIFT,RTOL
        ELSE
          WRITE(IOUT,1020) MDIM,NELEM,NGRID,NPOL,NGQ,IDIM,
     1                     IBOUND,SHIFT
        END IF
        WRITE(IOUT,1030)
        NO    = 0
        DO 20 I = 3 , NMESH , 2
           NO  = NO + 1
           NM  = IDINT(RMESH(I-1))
           XM  = (RMESH(I) - RMESH(I-2)) / NM
           XS  = XM / NPOL
           WRITE(IOUT,1040) NO,NM,RMESH(I-2),XM,XS,RMESH(I)
   20   CONTINUE
      END IF
C
      DO 300 NDF = 1 , NDIR
C
        NDIM = NDIL(NDF)
        NNODES = NNODE * NDIM
        NSTF   = NNODES * (NNODES + 1) / 2
C
        IF (IPRINT .GE. 0) THEN
          WRITE(IOUT,1060)
          WRITE(IOUT,1090) NDIM, MDIM
          WRITE(*   ,1060)
          WRITE(   *,1090) NDIM, MDIM
        END IF
C
        NN    = NDIM * NGRID
        NNE   = NN
        NNM   = NN + 1
C
        INODE  = 1
        INODES = INODE  +  NNODE  * NELEM
        IMAXA  = INODES +  NNODES * NELEM
        IMHT   = IMAXA  +  NNM
        INHT   = IMHT   +  NNM
        ILAST  = INHT   +  NNM
C
C          Allocate space
C
        ALLOCATE (ITOT(ILAST), STAT=ISTAT)
        IF (ISTAT .NE. 0) THEN
            WRITE (*,2010) 'ITOT', 'NODGEN', ISTAT
            DEALLOCATE (XGT)
            DEALLOCATE (WGT)
            DEALLOCATE (GRT)
            STOP
        ENDIF
        DO I = 1 , ILAST
          ITOT(I) = 0
        END DO
C
C      GENERATE AND STORE NODAL DATA
C
        CALL NODGEN(ITOT(INODE),NNODE,ITOT(INODES),NNODES,NDIM,NELEM,
     1              IPRINT,IOUT)
C
        IF (IPRINT .GE. 2) WRITE(IOUT,1070) ILAST
C
C
C      SET THE DIRICHLET AND/OR NEUMANN BOUNDARY CONDITIONS
C
        IBOUN = IBOUND
        IF (IBOUND .GE. 6) IBOUN = IBOUN - 4
        CALL BOUNDC(ITOT(INODE),ITOT(INODES),ITOT(IMAXA),ITOT(IMHT),
     1              ITOT(INHT),IBOUN,NNODE,NNODES,NELEM,NDIM,NN,NNE,
     2              IPRINT,IOUT)
C
C      CALCULATE ADDRESSES OF DIAGONAL ELEMENTS
C      IN BANDED STIFFNESS AND MASS MATRICES
C
        CALL MAXHT(ITOT(IMAXA),ITOT(IMHT),ITOT(INODES),NNODES,NELEM,
     1             NN,NNM,NWK,NWM,MK,IPRINT,IOUT)
C
        MMK  = NWK / NN
        IF (IPRINT .GE. 0) WRITE(IOUT,1050) NN,NWK,MK,MMK
C
        NWT = NWK
      
        IF (PROBTYPE .EQ. 'SCATTER') THEN
C
C
C      **********************************************************
C      ***                                                    ***
C      ***    E L E M E N T  G E N E R A T I N G  P H A S E   ***
C      ***                                                    ***
C      **********************************************************
C
C      CALCULATE AND ASSEMBLE STIFFNESS AND MASS
C      ELEMENT MATRICES INTO THE GLOBAL ONES
C
          IPSI   = IGRID  +  NGRID
          IDPSI  = IPSI   +  NGQ    * NNODE
          IEMS1  = IDPSI  +  NGQ    * NNODE
          IXNODE = IEMS1  +  NNODE  * NNODE
          IUNOD1 = IXNODE +  NNODE
          IQNOD1 = IUNOD1 +  NNODE
          IAA    = IQNOD1 +  NNODE
          IEIGL  = IAA    +  NDIM   * NDIM
          IEIGR  = IEIGL  +  NDIM   * NDIM
          ILAST  = IEIGR  +  NDIM   * NDIM
C
          IF (IPRINT .GE. 2) WRITE(IOUT,1080) ILAST
C
C        Allocate space
C
          ALLOCATE (TOT(ILAST), STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2010) 'TOT', 'ASSMBC', ISTAT
              DEALLOCATE (XGT)
              DEALLOCATE (WGT)
              DEALLOCATE (GRT)
              DEALLOCATE (ITOT)
              STOP
          ENDIF
          DO I = 1 , ILAST
            TOT(I) = DZERO
          END DO
C
C          Copy data from earlier allocated arrays to the TOT array
C
          DO 35 I = 1 , NGQ
             TOT(IXG + I - 1) = XGT(I)
             TOT(IWG + I - 1) = WGT(I)
   35     CONTINUE
          DO 45 I = 1 , NGRID
             TOT(IGRID + I - 1) = GRT(I)
   45     CONTINUE
C
C                COMPLEX ARRAY
C
C
C         Allocate space
C
          ALLOCATE (ZTOTM(NWK), STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2010) 'ZTOTM', 'ASSMBC', ISTAT
              DEALLOCATE (XGT)
              DEALLOCATE (WGT)
              DEALLOCATE (GRT)
              DEALLOCATE (TOT)
              DEALLOCATE (ITOT)
              STOP
          ENDIF

          NOPENL = 0
C
          NOPENR = 0
          DO I = 1 , NDIM
            IF (SHIFT .GT. THRSHR(I)) NOPENR = NOPENR + 1
CX            WRITE(*,*) I,THRSHR(I)
          END DO
C
          NOPEN = MAX(NOPENL,NOPENR)

          NL = MAX(1,NOPENL)
          NR = MAX(1,NOPENR)
          IF (IPRINT .GE. 0) THEN
            WRITE(*   ,*)
            WRITE(IOUT,*)
            WRITE(   *,2105) NOPENR
            WRITE(IOUT,2105) NOPENR
            DO I = 1 , NOPENR
              WRITE(   *,2110) I,DSQRT(SHIFT - THRSHR(I))
              WRITE(IOUT,2110) I,DSQRT(SHIFT - THRSHR(I))
            END DO
          END IF
C
          IUPOT  = 1          
          IQPOT  = IUPOT   +  MDIM   * MDIM   * NGQ
          IVV    = IQPOT   +  MDIM   * MDIM   * NGQ
          IQQ    = IVV     +  MDIM   * MDIM
          IESTFA = IQQ     +  MDIM   * MDIM
          IESTF  = IESTFA  +  NSTF
          IEST1  = IESTF   +  NNODES * NNODES
          IUNODE = IEST1   +  NNODE  * NNODE
          IQNODE = IUNODE  +  NNODE
          IR     = IQNODE  +  NNODE
          IPIRR  = IR      +   NN     * NOPEN
          IPREG  = IPIRR   +   NDIM   * NOPEN
          IDPIRR = IPREG   +   NDIM   * NOPEN
          IDPREG = IDPIRR  +   NDIM   * NOPEN
          IRR_R  = IDPREG  +   NDIM   * NOPEN
          ITT_R  = IRR_R   +   NL     * NL
          IRR_L  = ITT_R   +   NR     * NL
          ITT_L  = IRR_L   +   NR     * NR
          IWW    = ITT_L   +   NL     * NR
          IWS    = IWW     +   NDIM   * NOPEN
          IRR0   = IWS     +   NDIM   * NOPEN
          IFAC   = IRR0    +   NDIM   * NDIM
          IWL    = IFAC    +   NN     * NDIM
          IWR    = IWL     +   NDIM   * NDIM
          ITEMP  = IWR     +   NDIM   * NDIM
          ILAST  = ITEMP   +   NOPEN  * NOPEN
C
          IF (IPRINT .GE. 2) WRITE(IOUT,1085) ILAST+NWK
C
C         Allocate space
C
          ALLOCATE (ZTOT(ILAST), STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2010) 'ZTOT', 'ASSMBC', ISTAT
              DEALLOCATE (XGT)
              DEALLOCATE (WGT)
              DEALLOCATE (GRT)
              DEALLOCATE (TOT)
              DEALLOCATE (ITOT)
              DEALLOCATE (ZTOTM)
              STOP
          ENDIF
          DO I = 1 , ILAST
            ZTOT(I) = ZERO
          END DO
C
C      CALCULATE AND ASSEMBLE STIFFNESS MATRIX
C
          CALL ASSMBC(ITOT(IMAXA),ZTOTM,ZTOT(IESTF),ZTOT(IESTFA),
     1                ZTOT(IEST1),TOT(IEMS1),TOT(IGRID),ZTOT(IUPOT),
     2                ZTOT(IQPOT),TOT(IXNODE),ZTOT(IUNODE),ZTOT(IQNODE),
     3                TOT(IUNOD1),TOT(IQNOD1),      
     3                TOT(IPSI),TOT(IDPSI),TOT(IXG),TOT(IWG),ZTOT(IVV),
     4                ZTOT(IQQ),ITOT(INODE),ITOT(INODES),NMDIL,NELEM,
     5                NNODE,NNODES,NGQ,NN,NNM,NWK,NWM,NWT,NDIM,MDIM,
     6                IDIM,NGRID,NSTF,IPRINT,IOUP,SHIFT,IPRSTP,IOUM,
     7                IOUT,ISC)
C
C           *****************************************
C           ***                                   ***
C           ***       SCATTERING PROBLEM          ***
C           ***                                   ***
C           *****************************************
C
          IF (ISC .EQ. 1 .OR. ISC .EQ. 3) THEN
C
            IF (IPRINT .GE. 0) THEN
              WRITE(IOUT,1060)
              WRITE(IOUT,1110)
              WRITE(*   ,1060)
              WRITE(*   ,1110)
            END IF
C
            IF (NOPENR .NE. 0) THEN
C
              IZER = 1
              IF (NOPENL .EQ. 0) IZER = 0
              CALL SCSOLC(ITOT(IMAXA),ZTOT(IR),ZTOT(IQQ),ZTOT(IFAC),
     1                    TOT(IGRID),TOT(IEIGL),TOT(IEIGR),NGRID,NN,NNM,
     2                    NWK,SHIFT,THRSHL,THRSHR,NDIM,MDIM,IDIM,IPRINT,
     3                    IPRSTP,IOUT,IOUF,NOPEN,NOPENL,NOPENR,NOPENL,
     4                    NOPENR,ZTOT(IPIRR),ZTOT(IPREG),ZTOT(IDPIRR),
     5                    ZTOT(IDPREG),ZTOT(IRR_L),ZTOT(ITT_L),
     6                    ZTOT(IWL),ZTOT(IWW),ZTOT(IWS),ZTOT(IRR0),
     7                    ZTOTM,ZTOT(IWR),ZTOT(ITEMP),IOUM,ISC,1,IZER)
C
            ELSE
C
              IF (IPRINT .GE. 0) THEN
                  WRITE(IOUT,2130)
                  WRITE(   *,2130)
              END IF
C
            END IF
C
          END IF
C
          DEALLOCATE (ZTOT, STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2020) 'ZTOT', ISTAT
          ENDIF
  321     DEALLOCATE (TOT, STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2020) 'TOT', ISTAT
          ENDIF
          DEALLOCATE (ITOT, STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2020) 'ITOT', ISTAT
          ENDIF
          DEALLOCATE (ZTOTM, STAT=ISTAT)
          IF (ISTAT .NE. 0) THEN
              WRITE (*,2020) 'ZTOTM', ISTAT
          ENDIF
C
        END IF
C
  300 CONTINUE
C
C       Deallocate used arrays
C
      DEALLOCATE (XGT, STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2020) 'XGT', ISTAT
      ENDIF
      DEALLOCATE (WGT, STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2020) 'WGT', ISTAT
      ENDIF
      DEALLOCATE (GRT, STAT=ISTAT)
      IF (ISTAT .NE. 0) THEN
          WRITE (*,2020) 'GRT', ISTAT
      ENDIF
C
      CLOSE (UNIT=IOUP)
      CLOSE (UNIT=IOUM)
      RETURN
 1000 FORMAT(/1X,'PROBLEM: ',A70/1X,8('*')/)
 1010 FORMAT(/15X,'C O N T R O L  I N F O R M A T I O N'/
     *        15X,'------------------------------------'//5X,
     * 'NUMBER OF DIFFERENTIAL EQUATIONS. . . . . (MDIM  ) = ',I6/5X,
     * 'NUMBER OF ENERGY LEVELS REQUIRED. . . . . (NROOT ) = ',I6/5X,
     * 'NUMBER OF FINITE ELEMENTS . . . . . . . . (NELEM ) = ',I6/5X,
     * 'NUMBER OF GRID POINTS . . . . . . . . . . (NGRID ) = ',I6/5X,
     * 'ORDER  OF SHAPE FUNCTIONS . . . . . . . . (NPOL  ) = ',I6/5X,
     * 'ORDER  OF GAUSS-LEGENDRE QUADRATURE . . . (NGQ   ) = ',I6/5X,
     * 'NUMBER OF SUBSPACE ITERATION VECTORS. . . (NC    ) = ',I6/5X,
     * 'DIMENSION OF ENVELOPE SPACE . . . . . . . (IDIM  ) = ',I6/5X,
     * 'BOUNDARY CONDITION CODE . . . . . . . . . (IBOUND) = ',I6/5X,
     * 'SHIFT OF DOUBLE ENERGY SPECTRUM . . . . . (SHIFT ) = ',G14.6/5X,
     * 'CONVERGENCE TOLERANCE . . . . . . . . . . (RTOL  ) = ',G14.6/)
 1020 FORMAT(/15X,'C O N T R O L  I N F O R M A T I O N'/
     *        15X,'------------------------------------'//5X,
     * 'NUMBER OF DIFFERENTIAL EQUATIONS. . . . . (MDIM  ) = ',I6/5X,
     * 'NUMBER OF FINITE ELEMENTS . . . . . . . . (NELEM ) = ',I6/5X,
     * 'NUMBER OF GRID POINTS . . . . . . . . . . (NGRID ) = ',I6/5X,
     * 'ORDER OF SHAPE FUNCTIONS. . . . . . . . . (NPOL  ) = ',I6/5X,
     * 'ORDER OF GAUSS-LEGENDRE QUADRATURE. . . . (NGQ   ) = ',I6/5X,
     * 'DIMENSION OF ENVELOPE SPACE . . . . . . . (IDIM  ) = ',I6/5X,
     * 'BOUNDARY CONDITION CODE . . . . . . . . . (IBOUND) = ',I6/5X,
     * 'DOUBLE ENERGY SPECTRUM. . . . . . . . . . (SHIFT ) = ',G14.6/)
 1030 FORMAT(/4X,
     * 'SUBDIVISION OF RHO-REGION ON THE FINITE-ELEMENT GROUPS:'/4X,
     * '****************************************************** ',//
     *2X,' NO OF  NUMBER OF  BEGIN OF  LENGTH OF    GRID     END OF  '/
     *2X,' GROUP  ELEMENTS   INTERVAL   ELEMENT     STEP    INTERVAL'/
     *2X,' -----  ---------  --------  ---------  --------  --------')
 1040 FORMAT(2X,I4,4X,I6,3X,F10.3,1X,F9.5,1X,F9.5,1X,F10.3)
 1050 FORMAT(//15X,'T O T A L  S Y S T E M  D A T A  '/
     *         15X,'-------------------------------  '//5X,
     * 'TOTAL NUMBER OF ALGEBRAIC EQUATIONS. . . . (NN ) =',I9/5X,
     * 'TOTAL NUMBER OF MATRIX ELEMENTS. . . . . . (NWK) =',I9/5X,
     * 'MAXIMUM HALF BANDWIDTH . . . . . . . . . . (MK ) =',I9/5X,
     * 'MEAN    HALF BANDWIDTH . . . . . . . . . . (MMK) =',I9/)
 1060 FORMAT(/80('*')/)
 1070 FORMAT(/5X,'LAST ADDRESS OF ARRAY ITOT USED = ',I10/)
 1080 FORMAT(/5X,'LAST ADDRESS OF ARRAY TOT USED = ',I10/)
 1085 FORMAT(/5X,'LAST ADDRESS OF ARRAY ZTOT USED = ',I10/)
 1090 FORMAT( 5X,'NDIM, MDIM=',2I10)
 1110 FORMAT(5X,'CALCULATION OF WAVE FUNCTION WITH DIRECTION <--'/)
 1120 FORMAT(5X,'CALCULATION OF WAVE FUNCTION WITH DIRECTION -->'/)
 2010 FORMAT(/"Couldn't allocate storage for array ",A4," before the ",
     *        A6, " call: STAT = ", I3/)
 2020 FORMAT(/"Couldn't deallocate storage for array ", A4,
     * ": STAT = ",I3/)
 2100 FORMAT(5X,
     /   'NUMBER OF LEFT OPEN CHANNELS. . . . . .  (NOPENL) = ',I6)
 2105 FORMAT(5X,
     /   'NUMBER OF RIGHT OPEN CHANNELS.  . . . .  (NOPENR) = ',I6)
 2110 FORMAT(5X,
     /   'VALUE OF I-TH MOMENTUM . . . . . . . . . (I,QR ) = ',I6,E12.4)
 2120 FORMAT(/5X,'NO LEFT AND RIGHT OPEN CHANNELS')
 2130 FORMAT(/5X,'NO RIGHT OPEN CHANNELS')
 2140 FORMAT(/5X,'NO LEFT OPEN CHANNELS')
      END
C
      SUBROUTINE CHECRT(NOPEN,NOPENL,NOPENR,ZR_L,ZT_L,ZR_R,ZT_R,ZR,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CHECK THE PROPERTIES OF REFLECTION AND        .
C .                  TRANSMISSION MATRICES                            .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Y), COMPLEX(8) (Z)
      DIMENSION ZR_L(NOPENR,NOPENR),ZT_L(NOPENL,NOPENR),ZR(NOPEN,NOPEN)
      DIMENSION ZR_R(NOPENL,NOPENL),ZT_R(NOPENR,NOPENL)
      DATA DZERO / 0.D0 /
C
      WRITE(IOUT,1010)
      WRITE(IOUT,1000)
      WRITE(*   ,1010)
      WRITE(*   ,1000)
C
      ERROR = DZERO
      DO I = 1 , NOPENL
        DO J = 1 , NOPENR
          ZR(I,J) = DZERO
          DO L = 1 , NOPENR
            ZR(I,J) = ZR(I,J) + DCONJG(ZT_R(L,I)) * ZR_L(L,J)
          END DO
          DO L = 1 , NOPENL
            ZR(I,J) = ZR(I,J) + DCONJG(ZR_R(L,I)) * ZT_L(L,J)
          END DO
          ERROR = DMAX1(ERROR,CDABS(ZR(I,J))    )
        END DO
      END DO
C
      WRITE(IOUT,1050)
      WRITE(*   ,1050)
      DO I = 1 , NOPENL
        WRITE(IOUT,1020) (DREAL(ZR(I,J)),J=1,NOPENR)
        WRITE(*   ,1020) (DREAL(ZR(I,J)),J=1,NOPENR)
      END DO
      WRITE(IOUT,1060)
      WRITE(*   ,1060)
      DO I = 1 , NOPENL
        WRITE(IOUT,1020) (DIMAG(ZR(I,J)),J=1,NOPENR)
        WRITE(*   ,1020) (DIMAG(ZR(I,J)),J=1,NOPENR)
      END DO
      WRITE(IOUT,1130) ERROR
      WRITE(*   ,1130) ERROR
C
      ERROR = DZERO
      DO I = 1 , NOPENL
        DO J = 1 , NOPENR
          ZR(I,J) = ZT_R(J,I)-ZT_L(I,J)
          ERROR = DMAX1(ERROR,CDABS(ZR(I,J))    )
        END DO
      END DO
      WRITE(IOUT,1110)
      WRITE(*   ,1110)
      DO I = 1 , NOPENL
        WRITE(IOUT,1020) (DREAL(ZR(I,J)),J=1,NOPENR)
        WRITE(*   ,1020) (DREAL(ZR(I,J)),J=1,NOPENR)
      END DO
      WRITE(IOUT,1120)
      WRITE(*   ,1120)
      DO I = 1 , NOPENL
        WRITE(IOUT,1020) (DIMAG(ZR(I,J)),J=1,NOPENR)
        WRITE(*   ,1020) (DIMAG(ZR(I,J)),J=1,NOPENR)
      END DO
      WRITE(IOUT,1130) ERROR
      WRITE(*   ,1130) ERROR
      WRITE(IOUT,1010)
      WRITE(*   ,1010)
      RETURN
 1000 FORMAT(/10X,'     C H E C K  P R O P E R T I E S     '/
     *        10X,'----------------------------------------')
 1010 FORMAT(/80('*')/)
 1020 FORMAT(1X,20(1X,G14.6))
 1050 FORMAT(/8X,'  R E  P A R T:  TT_->^1 * RR_<- + RR_->^1 * TT_<-  '/
     *        8X,'----------------------------------------------------')
 1060 FORMAT(/8X,'  I M  P A R T:  TT_->^1 * RR_<- + RR_->^1 * TT_<-  '/
     *        8X,'----------------------------------------------------')
 1110 FORMAT(/8X,'  R E  P A R T:  TT_->^T - TT_<-  '/
     *        8X,'----------------------------------')
 1120 FORMAT(/8X,'  I M  P A R T:  TT_->^T - TT_<-  '/
     *        8X,'----------------------------------')
 1130 FORMAT(/8X,'M A X I M A L  A B S O L U T E  E R R O R =',G14.6)
      END
C
      SUBROUTINE CHECKN(PROBTYPE,SHIFT,THRSHL,RTOL,EPSY,MDIM,NDIR,NDIL,
     1                  IDIM,NROOT,IBOUND,IOUT,IOUP,IOUM)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CHECK ALL INPUT DATA                          .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION NDIL(NDIR),THRSHL(MDIM)
      CHARACTER*70 PROBTYPE
      IF (NDIR .GT. MDIM .OR. NDIL(NDIR) .GT. MDIM) THEN
           WRITE(IOUT,1120) NDIR,NDIL(NDIR)
           STOP
      END IF
      IF (IBOUND .LT. 1 .OR. IBOUND .GT. 8 .OR.
     1    IBOUND .EQ. 5 .OR. IBOUND .EQ. 7) THEN
           WRITE(IOUT,1130) IBOUND
           STOP
      END IF
      IF (IBOUND .GE. 6 .AND. NROOT .GT. 1 .AND.
     1    PROBTYPE .EQ. 'EIGEN') THEN
           WRITE(IOUT,1140) IBOUND,NROOT,PROBTYPE
           STOP
      END IF
      IF (IBOUND .LT. 6 .AND. PROBTYPE .EQ. 'SCATTER') THEN
           WRITE(IOUT,1150) IBOUND,PROBTYPE
           STOP
      END IF
      IF (PROBTYPE .EQ. 'SCATTER') THEN
         DO I = 1 , MDIM
           IF (DABS(SHIFT - THRSHL(I)) .LT. EPSY) THEN
             WRITE(IOUT,1160) SHIFT,THRSHL(I),PROBTYPE
             STOP
           END IF
         END DO
      END IF
      IF (PROBTYPE .EQ. 'EIGEN' .AND. RTOL .LT. EPSY) THEN
           WRITE(IOUT,1170) RTOL,EPSY
           RTOL = EPSY
           WRITE(IOUT,1180) RTOL
      END IF
      IF (IDIM .NE. 1 .AND. PROBTYPE .EQ. 'SCATTER') THEN
           WRITE(IOUT,1190) IDIM,PROBTYPE
           STOP
      END IF
      IF (IOUT .EQ. 0) THEN
           IOUT = 7
      END IF
      IF (IOUP .EQ. 0) THEN
           IOUP = 10
      END IF
      IF (IOUM .EQ. 0) THEN
           IOUM = 11
      END IF
      RETURN
 1120 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF NDIR AND NDIL : ',2I4/)
 1130 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF IBOUND : ',I8/)
 1140 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF IBOUND, NROOT, PROBTYPE : ',
     *                                               2I4,2X,A70/)
 1150 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF IBOUND, PROBTYPE : ',I4,2X,A70/)
 1160 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF SHIFT, THRSH, PROBTYPE : ',
     *                                            2G12.4,2X,A70/)
 1170 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'WARNING: RTOL IS SMALLER THAN EPSY :',2G16.4)
 1180 FORMAT(/10X,'WARNING: HAS CHANGED RTOL : ',G12.4/)
 1190 FORMAT(/3X,'*****  KANTBP  ***** ',3X,
     *       'FORBIDDEN VALUE OF IDIM, PROBTYPE : ',I4,2X,A70/)
      END
C
      SUBROUTINE SCSOLC(MAXA,R,QQ,FAC,XGRID,EIGL,EIGR,NGRID,NN,NNM,NWK,
     1                  SHIFT,THRSHL,THRSHR,NDIM,MDIM,IDIM,IPRINT,
     2                  IPRSTP,IOUT,IOUF,NOPEN,NOPENL,NOPENR,NL,NR,PIRR,
     3                  PREG,DPIRR,DPREG,RR,TT,WL,WW,WS,RR0,ZZ,WR,TEMP,
     4                  IOUM,ISC,IDIR,IZER)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE REACTION MATRIX AND THE         .
C .                  EIGENSOLUTION AND PRINT OUT THE EIGENSOLUTION    .
C .                  CALCULATED AND STORE THE SOLUTIONS               .
C .                  OBTAINED INTO THE IOUF FILE, IF NECESSARY        .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION XGRID(NGRID),MAXA(NNM),THRSHL(MDIM),THRSHR(MDIM),
     1          EIGL(NDIM,NDIM),EIGR(NDIM,NDIM)
      COMPLEX(8) R(NN,NR),RR0(NDIM,NDIM),ZZ(NWK),QQ(MDIM,MDIM),
     1           PIRR(NDIM,NR),PREG(NDIM,NR),
     2           DPIRR(NDIM,NR),DPREG(NDIM,NR),
     3           TT(NL,NR),RR(NR,NR),WW(NDIM,NOPEN),WS(NDIM,NOPEN),
     4           FAC(NN,NDIM),WL(NDIM,NDIM),WR(NDIM,NDIM),
     5           TEMP(NOPEN,NOPEN)
      DATA DZERO / 0.D0 /, ONE / 1.D0 /
      COMMON/HION/AP,ZP,RP,AT,ZT,RT
      COMMON/CONST/HBAR,PI
      COMMON/TRANS/FTR,QTRANS,NTRANS
C
      ANMASS=938.D0
C
      IF (ISC .EQ. 3) THEN
        REWIND (IOUM)
        READ(IOUM)  (ZZ(I),I=1,NWK)
      END IF
C
      PI   = DACOS(-1.D0)
      IFPR = - 1
      ISH  =  1
C
      IF (IDIR .EQ. 1) THEN
C
        CALL ASYMSL(XGRID(1),XGRID(NGRID),NDIM,NOPENL,NOPENR,SHIFT,
     1              THRSHL,THRSHR,EIGL,EIGR,WL,WR,PREG,PIRR,DPREG,DPIRR,
     2              IOUT)
        XGI   = XGRID(NGRID) ** (IDIM - 1)
        XGI1  = XGRID(1) ** (IDIM - 1)
        KK = NN - NDIM
        KKM = 0
        ISGN = + 1
C
      ELSE IF (IDIR .EQ. 2) THEN
C
        CALL ASYMSR(XGRID(1),XGRID(NGRID),NDIM,NOPENL,NOPENR,SHIFT,
     1              THRSHL,THRSHR,EIGL,EIGR,WL,WR,PREG,PIRR,DPREG,DPIRR,
     2              IOUT)
        XGI  = XGRID(1) ** (IDIM - 1)
        XGI1  = XGRID(NGRID) ** (IDIM - 1)
        KK = 0
        KKM = NN - NDIM
        ISGN = - 1
C
      END IF
C
      DO I = 1 , NDIM
        DO J = 1 , NL
          WS(I,J) = WL(I,J)
        END DO
      END DO
C
      CALL DECOMC(ZZ,MAXA,NN,ISH,IFPR,IOUT)
      DO I = 1 , NN
        DO J = 1 , NDIM
          FAC(I,J) = DZERO
        END DO
      END DO
C
      DO I = 1 , NDIM
        FAC(KK + I,I) = ISGN * XGI
        CALL REDBAC(ZZ,FAC(1,I),MAXA,NN)
      END DO
C
      L = 0
      DO J = 1 , NDIM
        DO I = J , 1 , -1
          L = L + 1
          ZZ(L) = FAC(KK + I,J)
          WL(I,J) = DZERO
          WL(J,I) = DZERO
          IF (I .EQ. J) MAXA(J) = L
        END DO
        WL(J,J) = ONE
      END DO
      MAXA(NDIM + 1) = L + 1
      CALL DECOMC(ZZ,MAXA,NDIM,ISH,IFPR,IOUT)
C
      DO I = 1 , NDIM
        CALL REDBAC(ZZ,WL(1,I),MAXA,NDIM)
      END DO
C
      DO I = 1 , NDIM
        DO J = 1 , NDIM
          WL(I,J) = WL(I,J) + QQ(I,J)
        END DO
      END DO
C
C      CALCULATE THE WRONSKIAN
C
      DO I = 1 , NDIM
        DO J = 1 , NR
          R (I,J) =   DPIRR(I,J)
          WW(I,J) =   DPREG(I,J)
          DO L = 1 , NDIM
            R (I,J) = R (I,J) - QQ(I,L) * PIRR(L,J)
            WW(I,J) = WW(I,J) - QQ(I,L) * PREG(L,J)
          END DO
        END DO
      END DO
      DO I = 1 , NR
        DO J = 1 , NR
          RR(I,J) = DZERO
          DO K = 1 , NDIM
            RR(I,J) = RR(I,J) + PIRR(K,I) * WW(K,J) - R(K,I) * PREG(K,J)
          END DO
          RR(I,J) = RR(I,J) * XGI
        END DO
      END DO
C
C      PRINT THE WRONSKIAN
C
      IF (IPRINT .GE. 0) THEN
        WRITE(IOUT,1060)
        WRITE(*   ,1060)
        DO I = 1 , NR
          WRITE(IOUT,1020) (DIMAG(RR(I,J)),J=1,NR)
          WRITE(*   ,1020) (DIMAG(RR(I,J)),J=1,NR)
        END DO
      END IF
C
C      CALCULATE THE REFLECTION MATRIX
C
      DO I = 1 , NDIM
C      MODIFIED BY CHUKA 20.08.2024
        IF (I .LE. NDIM-NTRANS) THEN 
          RMASS1=AP*AT/(AP+AT)*ANMASS
        ELSE 
          RMASS1=(AP+2)*(AT-2)/(AP+AT)*ANMASS
        END IF  
        FACT = HBAR**2 / 2.D0 / RMASS1
        DO J = 1 , NR
          R (I,J) = - DPIRR(I,J) * FACT
          WW(I,J) =   DPREG(I,J) * FACT
          DO L = 1 , NDIM
            R (I,J) = R (I,J) + WL(I,L) * PIRR(L,J)
            WW(I,J) = WW(I,J) - WL(I,L) * PREG(L,J)
          END DO
        END DO
      END DO
C
      DO I = 1 , NR
        DO J = 1 , NR
          TEMP(I,J) = DZERO
          RR(I,J) = DZERO
          DO L = 1 , NDIM
            TEMP(I,J) = TEMP(I,J) + R (L,I) * R(L,J)
            RR(I,J) = RR(I,J) + R (L,I) * WW(L,J)
          END DO
        END DO
      END DO
C
      L = 0
      DO J = 1 , NR
        DO I = J , 1 , -1
          L = L + 1
          ZZ(L) = TEMP(I,J)
          IF (I .EQ. J) MAXA(J) = L
        END DO
      END DO
      MAXA(NR + 1) = L + 1
      CALL DECOMC(ZZ,MAXA,NR,ISH,IFPR,IOUT)
C
      DO I = 1 , NR
        CALL REDBAC(ZZ,RR(1,I),MAXA,NR)
      END DO
C
C      PRINT THE CALCULATED REFLECTION MATRIX
C
      IF (IPRINT .GE. 0) THEN
        WRITE(IOUT,1070)
        WRITE(*   ,1070)
        DO I = 1 , NR
          WRITE(IOUT,1020) (DREAL(RR(I,J)),J=1,NR)
          WRITE(*   ,1020) (DREAL(RR(I,J)),J=1,NR)
        END DO
        WRITE(IOUT,1075)
        WRITE(*   ,1075)
        DO I = 1 , NR
          WRITE(IOUT,1020) (DIMAG(RR(I,J)),J=1,NR)
          WRITE(*   ,1020) (DIMAG(RR(I,J)),J=1,NR)
        END DO
      END IF
C
C      CALCULATE THE RADIAL WAVEFUNCTIONS
C
      DO I = 1 , NDIM
        IJ = KK + I
        DO J = 1 , NR
          R(IJ,J) =  PREG(I,J)
          DO L = 1 , NR
            R(IJ,J) = R(IJ,J) + PIRR(I,L) * RR(L,J)
          END DO
        END DO
      END DO
C
      DO I = 1 , NDIM
        DO J = 1 , NR
          WW(I,J) = DZERO
          DO L = 1 , NDIM
            WW(I,J) = WW(I,J) + (WL(I,L) - QQ(I,L)) * R(KK + L,J)
          END DO
        END DO
      END DO
C
      DO II = 1 , NN - NDIM
        I = II
        IF (IDIR .EQ. 2) I = II + NDIM
        DO J = 1 , NR
          R(I,J) = DZERO
          DO L = 1 , NDIM
            R(I,J) = R(I,J) + FAC(I,L) * WW(L,J)
          END DO
        END DO
      END DO
C
      IF(IOUF .NE. 0) THEN
        NLR = NL
        IF (IZER .EQ. 0) NLR = 0
        WRITE(IOUF) NDIM,NN,NLR,NR,NGRID,((RR(I,J),I=1,NR ),J=1,NR)
      END IF
C
      RETURN
 1010 FORMAT(/80('*')/)
 1020 FORMAT(1X,9(1X,G20.12))
 1030 FORMAT(1X,F11.4,1X,9(1X,D11.4)/(11X,1X,9(1X,D11.4)))
 1040 FORMAT(1X,F11.4,1X,6(1X,D11.4)/(11X,1X,6(1X,D11.4)))
 1050 FORMAT(/8X,'Z',8X,' R E  P A R T  O F  F U N C T I O N S '/
     *        8X,'-',8X,'--------------------------------------')
 1055 FORMAT(/8X,'Z',8X,' I M  P A R T  O F  F U N C T I O N S '/
     *        8X,'-',8X,'--------------------------------------')
 1060 FORMAT(/8X,' I M  P A R T  O F  W R O N S K I A N '/
     *        8X,'--------------------------------------')
 1070 FORMAT(/8X,' R E  P A R T  O F  RR  M A T R I X '/
     *        8X,'------------------------------------')
 1075 FORMAT(/8X,' I M  P A R T  O F  RR  M A T R I X '/
     *        8X,'------------------------------------')
 1080 FORMAT(/8X,' R E  P A R T  O F  TT  M A T R I X '/
     *        8X,'------------------------------------')
 1085 FORMAT(/8X,' I M  P A R T  O F  TT  M A T R I X '/
     *        8X,'------------------------------------')
 2000 FORMAT(/10X,'     C H E C K  P R O P E R T I E S     '/
     *        10X,'----------------------------------------')
 2010 FORMAT(/80('*')/)
 2020 FORMAT(1X,10(1X,G14.6))
 2030 FORMAT(/8X,'  |RR|^2 + |TT|^2  '/
     *       8X,'-------------------------')
 2070 FORMAT(/8X,'  R E  P A R T:  RR^T - RR  '/
     *        8X,'----------------------------------')
 2080 FORMAT(/8X,'  I M  P A R T:  RR^T - RR  '/
     *        8X,'----------------------------------')
 2130 FORMAT(/8X,'M A X I M A L  A B S O L U T E  E R R O R =',G14.6)
      END
C
      SUBROUTINE ASSMBC(MAXA,A,ESTF,ESTFA,EST1,EMS1,XGRID,UPOT,QPOT,
     1                  XNODE,UNODE,QNODE,UNOD1,QNOD1,
     2                  PSI,DPSI,XG,WG,H,Q,NODE,
     2                  NODES,NMDIL,NELEM,NNODE,NNODES,NGQ,NN,NNM,
     3                  NWK,NWM,NWT,NDIM,MDIM,IDIM,NGRID,NSTF,IPRINT,
     4                  IOUP,SHIFT,IPRSTP,IOUM,IOUT,ISC)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE AND ASSEMBLE THE STIFFNESS AND      .
C .                  MASS ELEMENT MATRICES INTO THE GLOBAL            .
C .                  STIFFNESS AND MASS MATRICES                      .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      COMPLEX(8) A(NWT),UPOT(MDIM,MDIM,NGQ),QPOT(MDIM,MDIM,NGQ),
     1           H(MDIM,MDIM),Q(MDIM,MDIM),ESTF(NNODES,NNODES),
     2           ESTFA(NSTF),EST1(NNODE,NNODE),UNODE(NNODE),
     3           QNODE(NNODE)
      DIMENSION MAXA(NNM),EMS1(NNODE,NNODE),XGRID(NGRID),
     2          XNODE(NNODE),PSI(NNODE,NGQ),DPSI(NNODE,NGQ),
     3          UNOD1(NNODE),QNOD1(NNODE),
     4          XG(NGQ),WG(NGQ),NODE(NNODE,NELEM),
     5          NODES(NNODES,NELEM)
      DATA ZERO / 0.D0 /
C
      IF (NMDIL .EQ. 0) THEN
C
        CALL HQPOTN(XGRID,XG,NGQ,NGRID,UPOT,QPOT,H,Q,NELEM,MDIM,
     1              IOUT,IOUP,IPRINT,IPRSTP,NNODE,NODE,XNODE)
C
        NMDIL = 1
C
      END IF
C
      DO 15 IG = 1 , NGQ
C
          CALL SHAPEF(XG(IG),UNOD1,QNOD1,XNODE,NNODE)
C
          DO 14 IN = 1 , NNODE
             PSI (IN,IG) = UNOD1(IN)
             DPSI(IN,IG) = QNOD1(IN)
   14     CONTINUE
C
   15 CONTINUE
C
      IF (IPRINT .GE. 5) THEN
          WRITE(IOUT,1130)
          DO I = 1 , NGQ
             WRITE(IOUT,1140) (PSI(J,I),J=1,NNODE)
          END DO
          WRITE(IOUT,1150)
      END IF
C
C      BEGIN BIG LOOP ON ELEMENTS
C
      REWIND IOUP
C
C      CALCULATE MASS AND STIFFNESS MATRICES
C
C      INITIALIZE ARRAYS A(NWK), B(NWK)
C
      DO 10 I = 1 , NWT
         A(I) = ZERO
   10 CONTINUE
C
      DO 220 NEL = 1 , NELEM
C
        READ(IOUP) L,(((UPOT(I,J,IG),J=I  ,MDIM),I=1,MDIM),IG=1,NGQ),
     1               (((QPOT(I,J,IG),J=I+1,MDIM),I=1,MDIM),IG=1,NGQ)
C
         DO 40 II = 1 , NNODES
         DO 40 JJ = 1 , NNODES
            ESTF(II,JJ) = ZERO
   40    CONTINUE
C
         DO 50 LL = 1 , NNODE
            NE           = NODE(LL,NEL)
            XNODE(LL)    = XGRID(NE)
   50    CONTINUE
C
         IF (IPRINT .GE. 5) THEN
              WRITE(IOUT,1110) NEL
              WRITE(IOUT,1120) (I,XNODE(I),I=1,NNODE)
              WRITE(IOUT,1150)
         END IF
C
C      CALCULATE ELEMENT MATRICES FOR DIAGONAL
C      ELEMENTS OF MASS AND STIFFNESS MATRICES
C
         DO 100 IA = 1 , NDIM
C
            DO 60 NA = 1 , NNODE
               UNODE(NA) = UPOT(IA,IA,NA)
   60       CONTINUE
C
            CALL EMASSD(XNODE,XG,WG,NGQ,NNODE,EMS1,PSI,IDIM)
            CALL ESTIFD(XNODE,XG,WG,NGQ,NNODE,EST1,UNODE,PSI,DPSI,IDIM,
     1                                                         IA,NDIM)
C
            DO 70 NI = 1  , NNODE
            DO 70 NJ = NI , NNODE
C
               NEI   = (NI - 1) * NDIM + IA
               NEJ   = (NJ - 1) * NDIM + IA
C
               ESTF(NEI,NEJ) = ESTF(NEI,NEJ) + EST1(NI,NJ)
     1                               - SHIFT * EMS1(NI,NJ)
C
   70        CONTINUE
C
            IF (IPRINT .GE. 5) THEN
               WRITE(IOUT,1030) NEL
               DO 80 IP = 1 , NNODE
   80          WRITE(IOUT,1040) (EST1(IP,JP),JP=1,NNODE),
     1                          (EMS1(IP,JP),JP=1,NNODE)
            END IF
C
  100    CONTINUE
C
         IF (IPRINT .GE. 5) THEN
            WRITE(IOUT,1050) NEL
            DO 110 IP = 1 , NNODES
  110       WRITE(IOUT,1020) (ESTF(IP,JP),JP=1,NNODES)
         END IF
C
C      CALCULATE ELEMENT MATRIX FOR NONDIAGONAL
C      ELEMENTS OF STIFFNESS MATRIX
C
         DO 180 IA = 1 , NDIM - 1
         DO 180 IB = IA + 1 , NDIM
C
            DO 140 NAB = 1 , NNODE
               UNODE(NAB)   =   UPOT(IA,IB,NAB)
               QNODE(NAB)   =   QPOT(IA,IB,NAB)
  140       CONTINUE
C
            CALL ESTIFN(XNODE,XG,WG,NGQ,NNODE,EST1,UNODE,QNODE,
     1                  PSI,DPSI,IDIM)
C
            DO 150 NI = 1  , NNODE
            DO 150 NJ = NI , NNODE
C
               NEI  = (NI - 1) * NDIM + IA
               NEJ  = (NJ - 1) * NDIM + IB
C
               ESTF(NEI,NEJ) = ESTF(NEI,NEJ) + EST1(NI,NJ)
C
               IF (NJ .NE. NNODE) THEN
C
                 NEII = NEI + IB - IA
                 NEJJ = NEJ - IB + IA + NDIM
                 ESTF(NEII,NEJJ) = ESTF(NEII,NEJJ) + EST1(NJ+1,NI)
C
               END IF
C
  150       CONTINUE
C
            IF (IPRINT .GE. 5) THEN
               WRITE(IOUT,1070) NEL
               DO 160 IP = 1 , NNODE
  160          WRITE(IOUT,1040) (EST1(IP,JP),JP=1,NNODE)
            END IF
C
  180    CONTINUE
C
         IF (IPRINT .GE. 5) THEN
            WRITE(IOUT,1050) NEL
            DO 190 IP = 1 , NNODES
  190       WRITE(IOUT,1020) (ESTF(IP,JP),JP=1,NNODES)
         END IF
C
         NL  = 0
         DO 210 IL = 1  , NNODES
         DO 210 JL = IL , NNODES
C
            NL        = NL + 1
            ESTFA(NL) = ESTF(IL,JL)
C
  210    CONTINUE
C
C      ASSEMBLE ELEMENT MATRICES INTO GLOBAL
C      STIFFNESS AND MASS MATRICES
C
         CALL ADDVEK(A,NWT,MAXA,NNM,ESTFA,NSTF,NODES(1,NEL),NNODES,1)
C
  220 CONTINUE
C
      IF (ISC .EQ. 3) THEN
        REWIND(IOUM)
        WRITE(IOUM) (A(I),I=1,NWT)
      END IF
C
      READ(IOUP) L,((UPOT(I,J,1),J=I  ,MDIM),I=1,MDIM),
     1             ((Q   (I,J  ),J=I+1,MDIM),I=1,MDIM)
      DO I = 1 , MDIM
        DO J = I + 1 , MDIM
          Q(J,I) = - Q(I,J)
        END DO
        Q(I,I) = ZERO
      END DO
      READ(IOUP) L,((QPOT(I,J,1),J=I  ,MDIM),I=1,MDIM),
     1             ((H   (I,J  ),J=I+1,MDIM),I=1,MDIM)
      DO I = 1 , MDIM
        DO J = I + 1 , MDIM
          H(J,I) = - H(I,J)
        END DO
        H(I,I) = ZERO
      END DO
C
      IF (IPRINT .GE. 4) THEN
         WRITE(IOUT,1080) NN,NNM,NWK,NNODE,NNODES,NELEM,NDIM,NGRID
         WRITE(IOUT,1090) NWK
         WRITE(IOUT,1020) (A(I),I=1,NWK)
      END IF
C
      RETURN
 1020 FORMAT(6(2X,D11.4))
 1030 FORMAT(/5X,'ELEMENT STIFFNESS AND MASS MATRICES EST1 & ',
     1            'EMS1 (ELEMENT NO ',I3,' ) :'/)
 1040 FORMAT(5(3X,D13.6))
 1050 FORMAT(/5X,
     * 'ELEMENT STIFFNESS MATRIX ESTF (ELEMENT NO ',I3,' ) :'/)
 1070 FORMAT(/5X,
     * 'ELEMENT STIFFNESS MATRIX EST1 (ELEMENT NO ',I3,' ) :'/)
 1080 FORMAT(//15X,'DIMENSIONS SET IN ASSMBC ARE : ',//5X,
     1 'NN     = ',I5,5X,'NNM    = ',I5,5X,
     2 'NWK    = ',I5,5X,'NNODE  = ',I5,5X/5X,
     3 'NNODES = ',I5,5X,'NELEM  = ',I5,5X,
     4 'NDIM   = ',I5,5X,'NGRID  = ',I5,/)
 1090 FORMAT(//20X,'S T I F F N E S S  M A T R I X  A(',I5,') :'/)
 1110 FORMAT(/10X,'ARRAY XNODE (NEL = ',I3,' )'/)
 1120 FORMAT(5(1X,I2,2X,D11.4))
 1130 FORMAT(/15X,'SHAPE FUNCTIONS PSI(NNODE,NGQ) :'/)
 1140 FORMAT(6(1X,D13.6))
 1150 FORMAT(/)
      END
C
      SUBROUTINE EMASSD(XNODE,XG,WG,NGQ,NNODE,EMSS,PSI,IDIM)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE AN ELEMENT MASS MATRIX              .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION XNODE(NNODE),EMSS(NNODE,NNODE),XG(NGQ),WG(NGQ),
     1          PSI(NNODE,NGQ)
      DATA ZERO / 0.D0 /, TWO / 2.D0 /
C
C      INITIALIZE ELEMENT ARRAY
C
      DO 10 I = 1 , NNODE
      DO 10 J = 1 , NNODE
         EMSS(I,J) = ZERO
   10 CONTINUE
C
C      FORM JACOBI DETERMINANT
C
      DETJ = (XNODE(NNODE) - XNODE(1)) / TWO
      DETI = (XNODE(NNODE) + XNODE(1)) / TWO
C
C      BEGIN INTEGRATION POINT LOOP
C
      DO 50 IG = 1 , NGQ
C
        RG = DETJ * XG(IG) + DETI
        RR = RG ** (IDIM-1)
C
        DO 40 I = 1 , NNODE
        DO 40 J = I , NNODE
C
          EMSS(I,J) = EMSS(I,J)
     1              + PSI(I,IG) * PSI(J,IG) * RR * WG(IG) * DETJ
          EMSS(J,I) = EMSS(I,J)
C
   40  CONTINUE
C
   50 CONTINUE
      RETURN
      END
C
      SUBROUTINE ESTIFD(XNODE,XG,WG,NGQ,NNODE,EST,POT,PSI,DPSI,IDIM,
     1                                                         IA,NDIM)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE AN ELEMENT STIFFNESS MATRIX         .
C .                  FOR DIAGONAL PART                                .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION XNODE(NNODE),XG(NGQ),WG(NGQ),
     1          PSI(NNODE,NGQ),DPSI(NNODE,NGQ)
      COMPLEX(8) POT(NNODE),EST(NNODE,NNODE)
      DATA ZERO / 0.D0 /, TWO / 2.D0 /
C      ADDED BY CHUKA 20.08.2024
	COMMON/ANGULAR/L
      COMMON/HION/AP,ZP,RP,AT,ZT,RT
      COMMON/CONST/HBAR,PI
      COMMON/TRANS/FTR,QTRANS,NTRANS
C
      ANMASS=938.D0
C
      IF (IA .LE. NDIM-NTRANS) THEN 
        RMASS1=AP*AT/(AP+AT)*ANMASS
      ELSE 
        RMASS1=(AP+2)*(AT-2)/(AP+AT)*ANMASS
      END IF  
      FAC = HBAR**2 / 2.D0 / RMASS1
C
C      INITIALIZE ELEMENT ARRAY
C
      DO 10 I = 1 , NNODE
      DO 10 J = I , NNODE
C
         EST(I,J) = ZERO
C
   10 CONTINUE
C
C      FORM JACOBI DETERMINANT
C
      DETJ = (XNODE(NNODE) - XNODE(1)) / TWO
      DETI = (XNODE(NNODE) + XNODE(1)) / TWO
C
C      BEGIN INTEGRATION POINT LOOP
C
      DO 50 IG = 1 , NGQ
C
        RG = DETJ * XG(IG) + DETI
        RR = RG ** (IDIM - 1)
        VCENT = L * (L + 1D0) / RG ** 2
C
        DO 40 I = 1 , NNODE
        DO 40 J = I , NNODE
C
C      MODIFIED BY CHUKA 20.08.2024
         EST(I,J) = EST(I,J) + 
     1             (PSI(I,IG) * PSI(J,IG) * (POT(IG) + FAC * VCENT)
     1            + FAC * DPSI(I,IG) * DPSI(J,IG) / DETJ / DETJ)
     2            * RR * WG(IG) * DETJ
         EST(J,I) = EST(I,J)
C
   40   CONTINUE
C
   50 CONTINUE
      RETURN
      END
C
      SUBROUTINE ESTIFN(XNODE,XG,WG,NGQ,NNODE,EST,UPOT,QPOT,PSI,DPSI,
     1                  IDIM)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE AN ELEMENT STIFFNESS MATRIX         .
C .                  FOR NONDIAGONAL PART                             .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION XNODE(NNODE),XG(NGQ),WG(NGQ),
     1          PSI(NNODE,NGQ),DPSI(NNODE,NGQ)
      COMPLEX(8) EST(NNODE,NNODE),UPOT(NNODE),QPOT(NNODE)
      DATA ZERO / 0.D0 /, TWO / 2.D0 /
C
C      INITIALIZE ELEMENT ARRAY
C
      DO 10 I = 1 , NNODE
      DO 10 J = 1 , NNODE
C
         EST(I,J) = ZERO
C
   10 CONTINUE
C
C      FORM JACOBI DETERMINANT
C
      DETJ = (XNODE(NNODE) - XNODE(1)) / TWO
      DETI = (XNODE(NNODE) + XNODE(1)) / TWO
C
C      BEGIN INTEGRATION POINT LOOP
C
      DO 50 IG = 1 , NGQ
C
        RG = DETJ * XG(IG) + DETI
        RR = RG ** (IDIM-1)
C
        DO 40 I = 1 , NNODE
        DO 40 J = 1 , NNODE
C
          EST(I,J) = EST(I,J) +  RR * QPOT(IG) * (PSI(I,IG) * DPSI(J,IG)
     1             - DPSI(I,IG) * PSI (J,IG)) / DETJ * WG(IG) * DETJ
          EST(I,J) = EST(I,J) + PSI(I,IG) * UPOT(IG) * PSI(J,IG)
     1             * WG(IG) * DETJ * RR
   40   CONTINUE
C
   50 CONTINUE
      RETURN
      END
C
      SUBROUTINE HQPOTN(XGRID,XG,NGQ,NGRID,UPOT,QPOT,H,Q,NELEM,MDIM,
     1                  IOUT,IOUP,IPRINT,IPRSTP,NNODE,NODE,XNODE)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE POTENTIAL CURVES AND MATRIX         .
C .                  ELEMENTS AND WRITE TO FILE IOUP                  .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION XGRID(NGRID),NODE(NNODE,NELEM),XNODE(NNODE),XG(NGQ)
      COMPLEX(8) QPOT(MDIM,MDIM,NGQ),UPOT(MDIM,MDIM,NGQ),H(MDIM,MDIM),
     2          Q(MDIM,MDIM)
      DATA ZERO / 0.D0 /, TWO / 2.D0 /
C
      REWIND IOUP
      LC = 0
C
      DO 110 L = 1 , NELEM
C
        DO  LL = 1 , NNODE
          NE           = NODE(LL,L)
          XNODE(LL)    = XGRID(NE)
        END DO
C
        DETJ = (XNODE(NNODE) - XNODE(1)) / TWO
        DETI = (XNODE(NNODE) + XNODE(1)) / TWO
C
C      BEGIN INTEGRATION POINT LOOP
C
        DO 50 IG = 1 , NGQ
C
          RG = DETJ * XG(IG) + DETI
C
          CALL POTCAL(RG,H,Q,MDIM,IOUT)
C
          DO 170 I = 1 , MDIM
C
            DO 140 J = I , MDIM
C
              UPOT(I,J,IG) = H(I,J)
              UPOT(J,I,IG) = H(I,J)
C
              IF (I .EQ. J) THEN
C
                QPOT(I,I,IG) = ZERO
C
              ELSE
C
                QPOT(I,J,IG) =  Q(I,J)
                QPOT(J,I,IG) = -Q(I,J)
C
              END IF
C
  140       CONTINUE
C
  170     CONTINUE
C
C     PRINT POTENTIAL MATRIX ELEMENTS IN THE NGRID POINTS
C
          IF (IPRINT .GE. 2) THEN
            IF (LC / IPRSTP * IPRSTP .EQ. LC) THEN
              IF (MDIM .GT. 1) THEN
                WRITE(IOUT,1000)
                WRITE(IOUT,1010)  LC,RG
                DO 363 I = 1 , MDIM
 363            WRITE(IOUT,1020) (UPOT(I,J,IG),J=1,MDIM)
                WRITE(IOUT,1030)  LC,RG
                DO 364 I = 1 , MDIM
 364            WRITE(IOUT,1020) (QPOT(I,J,IG),J=1,MDIM)
              ELSE
                WRITE(IOUT,1040)
                WRITE(IOUT,1050) RG,UPOT(1,1,IG)
              END IF
            END IF
          END IF
C
          LC = LC + 1
C
   50   CONTINUE
C
        WRITE(IOUP) L,(((UPOT(I,J,IG),J=I  ,MDIM),I=1,MDIM),IG=1,NGQ),
     1                (((QPOT(I,J,IG),J=I+1,MDIM),I=1,MDIM),IG=1,NGQ)
C
  110 CONTINUE
C
         CALL POTCAL(XGRID(NGRID),H,Q,MDIM,IOUT)
         WRITE(IOUP) NGRID,((H(I,J),J=I  ,MDIM),I=1,MDIM),
     1                     ((Q(I,J),J=I+1,MDIM),I=1,MDIM)
         CALL POTCAL(XGRID(1    ),H,Q,MDIM,IOUT)
         WRITE(IOUP)      1,((H(I,J),J=I  ,MDIM),I=1,MDIM),
     1                      ((Q(I,J),J=I+1,MDIM),I=1,MDIM)
C
      ENDFILE IOUP
CX      STOP
C
      RETURN
 1000 FORMAT(//10X,'P O T E N T I A L  M A T R I C E S  V(I,J)',
     * ' A N D  Q(I,J):')
 1010 FORMAT(/5X,'V-MATRIX AT THE POINT NO = ',I5,
     * '  AND HYPERRADIUS RHO = ',F11.5)
 1020 FORMAT(6(2X,D11.4))
 1030 FORMAT(/5X,'Q-MATRIX AT THE POINT NO = ',I5,
     * '  AND HYPERRADIUS RHO = ',F11.5)
 1040 FORMAT(//10X,'P O T E N T I A L  M A T R I X  V(1,1):'//
     *        4(4X,'RHO',5X,'V(1,1)',2X))
 1050 FORMAT(4(1X,F7.3,2X,D11.4))
      END
C
      SUBROUTINE ADDVEK(A,NWT,MAXA,NNM,STIFF,NSTF,NODE,NNODE,KOD)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO ASSEMBLE ELEMENT VECTOR STIFF INTO GLOBAL     .
C .                  VECTOR A USING A COMPACT STORAGE FORM            .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION MAXA(NNM),NODE(NNODE)
      COMPLEX(8) A(NWT),STIFF(NSTF)
      NS = 0
      IF (KOD .NE. 1) NS = NWT / 2
      NDF = 0
      DO 20 I = 1 , NNODE
         NI = NODE(I)
         IF (NI .GT. 0) THEN
            MNI = MAXA(NI)
            MI  = I
            DO 10 J = 1 , NNODE
               NJ = NODE(J)
               IF (NJ .GT. 0) THEN
                  NIJ = NI - NJ
                  IF (NIJ .GE. 0) THEN
                     MIJ  = MNI + NIJ
                     MNIJ = MI
                     IF (J .GE. I) MNIJ = J + NDF
                     A(MIJ + NS) = A(MIJ + NS) + STIFF(MNIJ)
                  END IF
               END IF
               MI = MI + NNODE - J
   10       CONTINUE
         END IF
         NDF = NDF + NNODE - I
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE BOUNDC(NODE,NODES,MAXA,MHT,NHT,IBOUND,NNODE,NNODES,
     1                  NELEM,NDIM,NN,NNE,IPRINT,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO SET THE DIRICHLET AND NEUMANN BOUNDARY        .
C .                  CONDITIONS                                       .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION MAXA(1),MHT(1),NHT(1),NODE(NNODE,NELEM),
     1          NODES(NNODES,NELEM)
C
C      ESTABLISH BOUNDARY CONDITIONS
C
      GO TO (10,40,60,80),IBOUND
C
C      DIRICHLET - DIRICHLET BOUNDARY CONDITIONS
C
   10 IK    = 0
      NEL   = 1
      DO 20 IL = 1 , NDIM
         IK      = IK + 1
         NHT(IK) = NODES(IL,NEL)
   20 CONTINUE
      NEL  = NELEM
      DO 30 IL = 1 , NDIM
         IK       = IK + 1
         ILL      = NDIM * (NNODE - 1)
         NHT(IK)  = NODES(IL + ILL,NEL)
   30 CONTINUE
      GO TO 200
C
C      DIRICHLET - NEUMANN BOUNDARY CONDITIONS
C
   40 IK   = 0
      NEL  = 1
      DO 50 IL = 1 , NDIM
         IK       = IK + 1
         NHT(IK)  = NODES(IL,NEL)
   50 CONTINUE
      GO TO 200
C
C      NEUMANN - DIRICHLET BOUNDARY CONDITIONS
C
   60 IK    = 0
      NEL   = NELEM
      DO 70 IL = 1 , NDIM
         IK       = IK + 1
         ILL      = NDIM * (NNODE - 1)
         NHT(IK)  = NODES(IL + ILL,NEL)
   70 CONTINUE
      GO TO 200
C
C      NEUMANN - NEUMANN BOUNDARY CONDITION
C
   80 GO TO 300
C
  200 CONTINUE
C
C      ELIMINATE CORRESPONDING COLUMNS AND ROWS FOR THE
C      IMPLEMENTATION OF THE DIRICHLET BOUNDARY CONDITIONS
C      AND SET A NEW NUMERATION OF NODES
C
      DO 210 IL = 1 , NNE
         MAXA(IL) = IL
         MHT (IL) = 0
  210 CONTINUE
C
      DO 220 I = 1 , IK
         J        = NHT(I)
         MHT(J)   = 1
  220 CONTINUE
C
      DO 230 NEL = 1 , NELEM
      DO 230   I = 1 , NNODES
         K  = NODES(I,NEL)
         IF (K .EQ. 0) GO TO 230
         IF (MHT(K) .NE. 0) NODES(I,NEL) = 0
  230 CONTINUE
C
      DO 240 I = 1 , IK
         J       = NHT(I)
         MAXA(J) = 0
  240 CONTINUE
C
      NN  = 0
      DO 250 I = 1 , NNE
         MHT(I)   = 0
         IF (MAXA(I) .NE. 0) THEN
             NN       = NN + 1
             NHT(NN)  = I
             MHT(I)   = NN
         END IF
  250 CONTINUE
C
      DO 270 NEL = 1 , NELEM
         DO 260 I = 1 , NNODES
            K   = NODES(I,NEL)
            IF (K .NE. 0) THEN
               NODES(I,NEL) = MHT(K)
            END IF
  260    CONTINUE
  270 CONTINUE
C
  300 CONTINUE
C
      IF (IPRINT .GE. 3) THEN
         WRITE(IOUT,1020) NNODE,NELEM
         DO 320 NEL = 1 , NELEM
  320    WRITE(IOUT,1030) (NODES(IP,NEL),IP=1,NNODES)
         WRITE(IOUT,1040)
      END IF
      RETURN
 1020 FORMAT(//3X,'NODAL ARRAY NODES (',I3,',',I3,
     1       ') MODIFY ACCORDING TO DIRICHLET CONDITION :'/)
 1030 FORMAT(12I6)
 1040 FORMAT(/)
      END
C
      SUBROUTINE COLMHT(MHT,NNM,NNODE,NODE)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE COLUMN HEIGHTS                      .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION MHT(NNM),NODE(NNODE)
      NMAX = 1000000
      DO 10 I = 1 , NNODE
         NDI = NODE(I)
         IF (NDI .NE. 0) THEN
            IF ((NDI - NMAX) .LT. 0) NMAX = NDI
         END IF
   10 CONTINUE
      DO 20 I = 1 , NNODE
         NDI = NODE(I)
         IF (NDI .NE. 0) THEN
            NDIF = NDI - NMAX
            IF (NDIF .GT. MHT(NDI)) MHT(NDI) = NDIF
         END IF
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE FEGRID(RMESH,NMESH,RGRID,NGRID,NELEM,NPOL)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE NODAL POINTS RGRID(NGRID)           .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION RMESH(NMESH),RGRID(1)
      NG = 0
      NL = 0
      DO 20 I = 3 , NMESH, 2
          XM  = RMESH(I) - RMESH(I-2)
          XS  = XM / IDINT(RMESH(I-1)) / DFLOAT(NPOL)
          NM  = IDINT(XM / XS + 0.000000000001)
          NL  = NL + IDINT(RMESH(I-1))
          DO 10 J = 1 , NM
              NG            = NG + 1
              RGRID(NG)     = RMESH(I-2) + (J-1) * XS
   10     CONTINUE
   20 CONTINUE
      NELEM        = NL
      NGRID        = NG + 1
      RGRID(NGRID) = RMESH(NMESH)
      RETURN
      END
C
      SUBROUTINE MAXHT(MAXA,MHT,NODES,NNODES,NELEM,NN,NNM,NWK,NWM,
     1                 MK,IPRINT,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE ADDRESSES OF DIAGONAL ELEMENTS      .
C .                  IN BANDED MATRIX                                 .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION MAXA(NNM),MHT(NNM),NODES(NNODES,NELEM)
C
C      CALCULATE COLUMN HEIGHTS
C
      DO 10 IL = 1 , NN
         MHT (IL) = 0
         MAXA(IL) = 0
   10 CONTINUE
C
      DO 20 NEL = 1 , NELEM
         CALL COLMHT(MHT,NNM,NNODES,NODES(1,NEL))
   20 CONTINUE
C
C      SET ADDRESSES OF DIAGONAL ELEMENTS IN MAXA
C
      MAXA(1) = 1
      MAXA(2) = 2
      MK      = 0
      DO 30 IL = 1 , NN
         IF (MHT(IL) .GT. MK) MK = MHT(IL)
         MAXA(IL+1) = MAXA(IL) + MHT(IL) + 1
   30 CONTINUE
C
      MK   = MK + 1
      NNM  = NN + 1
      NWK  = MAXA(NNM) - 1
      NWM  = NWK
C
      IF (IPRINT .GE. 3) THEN
         WRITE(IOUT,1010) NNM
         WRITE(IOUT,1020) (MAXA(I),I=1,NNM)
         WRITE(IOUT,1030)
      END IF
C
      IF (IPRINT .GE. 3) THEN
 1010 FORMAT(/3X,'ADDRESSES OF DIAGONAL ELEMENTS SET IN PROFILE ',
     * 'MATRIX MAXA(',I5,') :'/)
 1020 FORMAT(10I7)
 1030 FORMAT(/)
      END IF
      RETURN
      END
C
      SUBROUTINE NODGEN(NODE,NNODE,NODES,NNODES,NDIM,NELEM,IPRINT,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE NODAL POINT DISTRIBUTION        .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION NODE(NNODE,NELEM),NODES(NNODES,NELEM)
C
C      NUMERATE THE NODES
C
      DO I = 1 , NNODES
        DO J = 1 , NELEM
          NODES(I,J) = 0
          IF (I .LE. NNODE) NODE(I,J) = 0
        END DO
      END DO
C
      DO 10 NI = 1 , NNODE
         NODE(NI,1) = NI
   10 CONTINUE
C
      DO 30 NEL = 2 , NELEM
         DO 20 ND = 1 , NNODE
            NODE(ND,NEL) = NODE(NNODE,NEL-1) + ND - 1
   20    CONTINUE
   30 CONTINUE
C
C      GENERATE NODAL NUMBERS FOR MULTICHANNEL CASE
C
      DO 50 NEL = 1 , NELEM
      DO 50 NL  = 1 , NNODE
         DO 40 KL = 1 , NDIM
            IL = (NL - 1) * NDIM + KL
            JL = NODE(NL,NEL) * NDIM + KL - NDIM
            NODES(IL,NEL) = JL
   40    CONTINUE
   50 CONTINUE
C
      IF (IPRINT .GE. 3) THEN
         WRITE(IOUT,1000)
         WRITE(IOUT,1010) NNODE,NELEM
         DO 160 NEL = 1 , NELEM
  160    WRITE(IOUT,1020) NEL,(NODE(IP,NEL),IP=1,NNODE)
         WRITE(IOUT,1030) NNODE,NELEM
         DO 170 NEL = 1 , NELEM
  170    WRITE(IOUT,1040) (NODES(IP,NEL),IP=1,NNODES)
         WRITE(IOUT,1050)
      END IF
      RETURN
 1000 FORMAT(//15X,'N O D A L   D A T A   I N F O R M A T I O N '/
     *         15X,'------------------------------------------- '/)
 1010 FORMAT(/20X,'NODAL ARRAY NODE (',I2,' , ',I4,') :'/)
 1020 FORMAT(5X,'ELEMENT NO ',I5,5X,10I5)
 1030 FORMAT(//20X,'NODAL ARRAY NODES (',I4,' , ',I4,') :'/)
 1040 FORMAT(12I6)
 1050 FORMAT(/)
      END
C
      SUBROUTINE SHAPEF(X,PSI,DPSI,XNODE,NPOL1)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE SHAPE FUNCTIONS PSI(I) AND THEIR    .
C .                  DERIVATIVES DPSI(I) WITH RESPECT TO THE MASTER   .
C .                  ELEMENT COORDINATE AT A SPECIFIED VALUE OF X     .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION PSI(NPOL1),DPSI(NPOL1),XNODE(NPOL1)
      DATA ZERO / 0.D0 /, ONE / 1.D0 /, TWO / 2.D0 /
      NPOL  = NPOL1 - 1
      HNODE = TWO / DFLOAT(NPOL)
      DO 10 IP = 1 , NPOL1
          XNODE(IP) = - ONE + HNODE * (IP - 1)
   10 CONTINUE
      DO 40 NP = 1 , NPOL1
          PSIN       = ONE
          PSI(NP)    = ONE
          DPSI(NP)   = ZERO
          DO 30 IP = 1 , NPOL1
              IF (IP .NE. NP) THEN
                  PSI(NP) = PSI(NP) * (X - XNODE(IP))
                  PSIN    = PSIN * (XNODE(NP) - XNODE(IP))
                  PSID    = ONE
                  DO 20 JP = 1 , NPOL1
                      IF (JP .NE. IP .AND. JP .NE. NP) THEN
                           PSID = PSID * (X - XNODE(JP))
                      END IF
   20             CONTINUE
                  DPSI(NP) = DPSI(NP) + PSID
              END IF
   30     CONTINUE
          PSI(NP)  = PSI(NP)  / PSIN
          DPSI(NP) = DPSI(NP) / PSIN
   40 CONTINUE
      RETURN
      END
C
C     PUBLISHED PROGRAMS
C
      SUBROUTINE GAULEG(X1,X2,X,W,N,EPS)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE NODES AND WEIGHTS OF THE        .
C .                  GAUSS-LEGENDRE QUADRATURE. PUBLISHED IN:         .
C .                  W.H. PRESS, B.F. FLANERY, S.A. TEUKOLSKY,        .
C .                  AND W.T.VETTERLEY, NUMERICAL RECIPES: THE        .
C .                  ART OF SCIENTIFIC COMPUTING (CAMBRIDGE           .
C .                  UNIVERSITY PRESS, CAMBRIDGE, 1986), CHAP. 4      .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION X(N),W(N)
      DATA ZERO/0.D0/, QUAT/0.25D0/, HALF/0.5D0/, ONE/1.D0/, TWO/2.D0/
      PI = DACOS(-ONE)
      M=(N+1)/2
      XM=HALF*(X2+X1)
      XL=HALF*(X2-X1)
      DO 12 I = 1 , M
        Z=DCOS(PI*(I-QUAT)/(N+HALF))
    1   CONTINUE
          P1=ONE
          P2=ZERO
          DO 11 J = 1 , N
            P3=P2
            P2=P1
            P1=((TWO*J-ONE)*Z*P2-(J-ONE)*P3)/J
   11     CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-ONE)
          Z1=Z
          Z=Z1-P1/PP
        IF (ABS(Z-Z1).GT.EPS) GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=TWO*XL/((ONE-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
   12 CONTINUE
      RETURN
      END
C
      SUBROUTINE REDBAC(ZA,ZV,MAXA,NN)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                 .
C .   P R O G R A M                                                 .
C .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS          .
C .                                                                 .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Y), COMPLEX(8) (Z)
      DIMENSION ZA(1),ZV(NN),MAXA(NN+1)
      DO 400 N = 1 , NN
        KL = MAXA(N) + 1
        KU = MAXA(N+1) - 1
        IF (KU - KL) 400,410,410
  410   K = N
        ZC = DCMPLX(0.D0,0.D0)
        DO 420 KK = KL , KU
          K = K - 1
  420   ZC = ZC + ZA(KK) * ZV(K)
        ZV(N) = ZV(N) - ZC
  400 CONTINUE
      DO 480 N = 1 , NN
        K = MAXA(N)
  480 ZV(N) = ZV(N) / ZA(K)
      IF (NN .EQ. 1) RETURN
      N = NN
      DO 500 L = 2 , NN
        KL = MAXA(N) + 1
        KU = MAXA(N+1) - 1
        IF (KU - KL) 500,510,510
  510   K = N
        DO 520 KK = KL , KU
          K = K - 1
  520   ZV(K) = ZV(K) - ZA(KK) * ZV(N)
  500 N = N - 1
      RETURN
      END
C
      SUBROUTINE DECOMC(ZA,MAXA,NN,ISH,IFPR,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .         TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF              .
C .         STIFFNESS MATRIX                                          .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Y), COMPLEX(8) (Z)
      DIMENSION ZA(1),MAXA(NN+1)
C
      IF (NN .EQ. 1) RETURN
      NLOW = 0
      DO 200 N = 1 , NN
        KN = MAXA(N)
        KL = KN+1
        KU = MAXA(N+1) - 1
        KH = KU - KL
        IF (KH) 200,240,210
  210   K = N - KH
        IC = 0
        KLT = KU
        DO 260 J = 1 , KH
          IC = IC + 1
          KLT = KLT - 1
          KI = MAXA(K)
          ND = MAXA(K+1) - KI - 1
          IF (ND) 260,260,270
  270     KK = MIN0(IC,ND)
          ZC = DCMPLX(0.D0,0.D0)
          DO 280 L = 1 , KK
  280     ZC = ZC + ZA(KI+L) * ZA(KLT+L)
          ZA(KLT) = ZA(KLT) - ZC
  260   K = K + 1
  240   K = N
        ZB = DCMPLX(0.D0,0.D0)
        DO 300 KK = KL , KU
          K = K - 1
          KI = MAXA(K)
          ZC = ZA(KK) / ZA(KI)
          IF (CDABS(ZC) .LT. 1.D16) GO TO 290
          WRITE (IOUT,2010) N,ZC
          STOP
  290     ZB = ZB + ZC * ZA(KK)
  300   ZA(KK) = ZC
        ZA(KN) = ZA(KN) - ZB
  200 CONTINUE
      RETURN
 2010 FORMAT (//47H STOP - STURM SEQUENSE CHECK FAILED BECAUSE OF    ,/
     18X,35HMULTIPLIER GROWTH FOR COLUMN NUMBER,I7,//
     212H MULTIPLIER=,2E20.8)
      END
C