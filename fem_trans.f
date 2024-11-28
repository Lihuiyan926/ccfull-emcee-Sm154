      SUBROUTINE FEM(CU,NLEVEL,RMIN00) 
      PARAMETER (NMESHM=111,MDIMM=111)
      PARAMETER (nlevelmax=100)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION RMESH(NMESHM),NDIL(MDIMM),THRSHL(MDIMM),THRSHR(MDIMM)
      COMPLEX(8) CU(NLEVEL,NLEVEL)
      CHARACTER*70 TITLE
      CHARACTER*55 FNOUT,POTEN,FMATR,EVWFN
      
	COMMON/DYN/E,RMASS
	COMMON/CONST/HBAR,PI
	COMMON/GRID/RMIN,RMAX,DR
      COMMON/COUP1/CPOT1(NLEVELMAX)
      COMMON /COUNTER / ICOUNT
      COMMON /PARA / IMODEL,RMAX_K,ngrid
C
C      SET INITIAL VALUES OF SOME CONSTANTS AND VARIABLES
C                  FOR SCATTERING PROBLEM
C
      NAMELIST /DIMS/  NMESH,MDIM,NDIR
      NAMELIST /PARAS/ TITLE,IPTYPE,ISC,ISCAT,IDIM,NPOL,SHIFT,IPRINT,
     1                 IPRSTP,RMESH,NDIL,NMDIL,THRSHL,IBOUND,
     2                 FNOUT,IOUT,POTEN,IOUP,FMATR,IOUM,EVWFN,IOUF
C
      IN = 110
      OPEN(UNIT=IN,FILE='TEST2.INP',STATUS='UNKNOWN',FORM='FORMATTED')
C
C     Read the first data block
C
      READ (IN,DIMS)
C
C     Read the second data block
C
      READ (IN,PARAS)
C
      MDIM = NLEVEL 
      NDIL = MDIM 
      RMESH(1) = 0D0
      RMESH(2) = 1500D0
      !RMESH(3) = 30D0 !c-wenpw2024.07.30
      RMESH(3) = RMAX_K
C      
      IF (NMESH .GT. NMESHM) THEN
         WRITE (*,2000) NMESH, NMESHM
         STOP
      ENDIF
      IF (NDIR .GT. MDIMM) THEN
         WRITE (*,2010) NDIR, NDIRR
         STOP
      ENDIF
      IF (MDIM .GT. MDIMM) THEN
         WRITE (*,2020) MDIM, MDIMM
         STOP
      ENDIF
C
C      MODIFIED BY CHUKA 20.08.2024      
CX      FAC = 2.D0 * RMASS / HBAR**2
      FAC = 1D0
      CALL CMAT1    
C
C      MODIFIED BY CHUKA 20.08.2024  
      DO I = 1 , MDIM
        THRSHR(I) = CPOT1(I) * FAC
      END DO
C
C      MODIFIED BY CHUKA 20.08.2024  
      SHIFT = E * FAC
      NMDIL = ICOUNT
      CALL KANTBP(TITLE,IPTYPE,ISC,ISCAT,NROOT,MDIM,IDIM,NPOL,
     1            RTOL,NITEM,SHIFT,IPRINT,IPRSTP,NMESH,RMESH,
     2            NDIR,NDIL,NMDIL,THRSHL,THRSHR,IBOUND,FNOUT,
     3            IOUT,POTEN,IOUP,FMATR,IOUM,EVWFN,IOUF)
C      
      CALL TRSQR(IOUF,CU,NLEVEL)
C
      CLOSE(IN)
      CLOSE(IOUT)
      CLOSE(IOUP)
      CLOSE(IOUM)
      CLOSE(IOUF)
      RETURN
C
 1000 FORMAT(/1X,'PROBLEM: ',A70/1X,8('*')/)
 2000 FORMAT(/"Requested size for array RMESH, ",I3,
     *        ", is larger the maximum size allowed, ", I3/)
 2010 FORMAT(/"Requested size for array NDIL, ",I3,
     *        ", is larger the maximum size allowed, ", I3/)
 2020 FORMAT(/"Requested size for array THRSHL, ",I3,
     *        ", is larger the maximum size allowed, ", I3/)
      END
C
      SUBROUTINE TRSQR(IOUF,CU,NLEVEL)
      IMPLICIT REAL(8) (A-H,O-Y)
      IMPLICIT COMPLEX(8) (Z)
      COMPLEX*16, ALLOCATABLE :: ZR(:,:),ZT(:,:)
      COMPLEX(8) CU(NLEVEL,NLEVEL)
      REWIND(IOUF)
      READ(IOUF) NDIM,NN,NOPEN1,NOPEN2,NGRID
      ALLOCATE (ZR(NOPEN2,NOPEN2),ZT(NOPEN1,NOPEN2))
      REWIND(IOUF)
      READ(IOUF) NDIM,NN,NOPEN1,NOPEN2,NGRID,
     1                              ((ZR(I,J),I=1,NOPEN2),J=1,NOPEN2)
      DO I = 1 , NOPEN2
          DO J = 1 , NOPEN2
             CU(I,J) = ZR(I,J)
          END DO   
      END DO

      DEALLOCATE (ZR,ZT)
      RETURN
      END
C
C     USER-SUPPLIED SUBROUTINES
C
      SUBROUTINE ASYMEV(ZMAX,NDIM,MDIM,SHIFT,DLAMBDA,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE INITIAL VALUE DLAMBDA FOR       .
C .                  THE HOMOGENEOUS THIRD TYPE BOUNDARY CONDITION    .
C .                  AT ZMAX                                          .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      RETURN
      END
C
      SUBROUTINE ASYMSL(ZMIN,ZMAX,NDIM,NOPENL,NOPENR,SHIFT,THRSHL,
     1                  THRSHR,EIGL,EIGR,PREGL,DREGL,PREGR,PIRRR,DREGR,
     2                  DIRRR,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE REGULAR, IRREGULAR              .
C .                  ASYMPTOTIC MATRIX SOLUTIONS PREGR, PIRRR         .
C .                  AND THEIR DERIVATIVES DREGR, DIRRR AT ZMAX,      .
C .                  THE REGULAR MATRIX SOLOTION PREGL AND ITS        .
C .                  DERIVATIVE DREGL AT ZMIN                         .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION  THRSHL(NDIM),THRSHR(NDIM),EIGL(NDIM,NDIM),
     1           EIGR(NDIM,NDIM)
      DIMENSION  FC1(0:3000),FCP1(0:3000),GC1(0:3000),GCP1(0:3000)
C
      COMPLEX(8) PREGL(NDIM,NDIM),DREGL(NDIM,NDIM),
     1           PREGR(NDIM,NOPENR),PIRRR(NDIM,NOPENR),
     1           DREGR(NDIM,NOPENR),DIRRR(NDIM,NOPENR),ZI
      COMMON/HION/AP,ZP,RP,AT,ZT,RT
      COMMON/CONST/HBAR,PI
CX      COMMON/DYN/E,RMASS
      COMMON/ANGULAR/L
C      ADDED BY CHUKA 20.08.2024 
      COMMON/TRANS/FTR,QTRANS,NTRANS
C
      ZI = (0.D0,1D0)
      EIGR(1:NDIM,1:NDIM) = 0D0
      DO I = 1 , NDIM
          EIGR(I,I) = 1D0
      END DO    
      ANMASS=938.D0
C
      DO J = 1 , NOPENR
C
C      MODIFIED BY CHUKA 20.08.2024
        IF (J .LE. NDIM-NTRANS) THEN 
          RMASS=AP*AT/(AP+AT)*ANMASS
        ELSE 
          RMASS=(AP+2)*(AT-2)/(AP+AT)*ANMASS
        END IF
        FAC = 2.D0 * RMASS / HBAR**2
        QR = (SHIFT - THRSHR(J)) * FAC 
        
        QR  = DSQRT(QR)
        SQR = DSQRT(QR)
        RHO = QR * ZMAX
        MINL=L
        MAXL=L    
        CHARGE = ZP * ZT * HBAR / 137.D0 * (RMASS/HBAR**2)
        ETA = CHARGE / QR
C
        ACCUR = 1D-14
        STEP = 100.D0
CX        ARGM = ARGGAMMA(L+1,ETA,1D-14)
        CALL RCWFNN(RHO,ETA,MINL,MAXL,FC1,FCP1,GC1,GCP1,ACCUR,STEP)
C       
C      MODIFIED BY CHUKA 20.08.2024
        DO I = 1 , NDIM
          PREGR(I,J) =  EIGR(I,J) * (GC1 (L) - ZI * FC1 (L)) 
     1               / SQR * SQRT(RMASS)
          PIRRR(I,J) = -EIGR(I,J) * (GC1 (L) + ZI * FC1 (L)) 
     1               / SQR * SQRT(RMASS)
          DREGR(I,J) =  EIGR(I,J) * (GCP1(L) - ZI * FCP1(L))
     1               * SQR * SQRT(RMASS)
          DIRRR(I,J) = -EIGR(I,J) * (GCP1(L) + ZI * FCP1(L))
     1               * SQR * SQRT(RMASS)
        END DO  
C
CX        PREGR(J,J) = PREGR(J,J) / CDEXP(-ZI * ARGM)
CX        PIRRR(J,J) = PIRRR(J,J) / CDEXP( ZI * ARGM)
CX        DREGR(J,J) = DREGR(J,J) / CDEXP(-ZI * ARGM)
CX        DIRRR(J,J) = DIRRR(J,J) / CDEXP( ZI * ARGM) 
C
      END DO


      RETURN
      END
C
      SUBROUTINE ASYMSR(ZMIN,ZMAX,NDIM,NOPENL,NOPENR,SHIFT,THRSHL,
     1                  THRSHR,EIGL,EIGR,PREGR,DREGR,PREGL,PIRRL,DREGL,
     2                  DIRRL,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE REGULAR, IRREGULAR              .
C .                  ASYMPTOTIC MATRIX SOLUTIONS PREGL, PIRRL         .
C .                  AND THEIR DERIVATIVES DREGL, DIRRL AT ZMIN,      .
C .                  THE REGULAR MATRIX SOLOTION PREGR AND ITS        .
C .                  DERIVATIVE DREGR AT ZMAX                         .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION  THRSHL(NDIM),THRSHR(NDIM),EIGL(NDIM,NDIM),
     1           EIGR(NDIM,NDIM)
      COMPLEX(8) PREGR(NDIM,NDIM),DREGR(NDIM,NDIM),
     1           PREGL(NDIM,NOPENL),PIRRL(NDIM,NOPENL),
     1           DREGL(NDIM,NOPENL),DIRRL(NDIM,NOPENL),ZI,ZP,ZM

      RETURN
      END
C
      SUBROUTINE POTCAL(RHO,VV,QQ,MDIM,IOUT)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  TO CALCULATE THE POTENTIAL MATRIX ELEMENTS       .
C .                  VV AND QQ OF DIMENSION MDIM X MDIM IN POINT      .
C .                  RHO                                              .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT REAL(8) (A-H,O-Z)
      COMPLEX(8) VV(MDIM,MDIM),QQ(MDIM,MDIM)
      PARAMETER (nlevelmax=100)
      DIMENSION CPOT0(NLEVELMAX,NLEVELMAX)
      DIMENSION CPOTW0(NLEVELMAX,NLEVELMAX)
      COMMON/CONST/HBAR,PI
      COMMON/DYN/E,RMASS
      COMMON/POT/V0,R0,A0,POW
      COMMON/IMAGP/W0,RW,AW
      COMMON/COUP1/CPOT1(NLEVELMAX)
      COMMON/COUP_R/CPOTR_REAL(NLEVELMAX,NLEVELMAX)
     1             ,CPOTR_IMG(NLEVELMAX,NLEVELMAX)   !ADDED-WENPW-20240326 
C
      DO I1 = 1 , MDIM 
        DO I2 = I1 , MDIM 
          QQ(I1,I2) = 0.D0
        ENDDO
      ENDDO
C      
      CALL CMAT(RHO,CPOT0,CPOTW0)  
C
CX      FAC = 2.D0 * RMASS / HBAR ** 2 
      FAC = 1D0
      DO I = 1 , MDIM
        DO J = I , MDIM
          VV(I,J) = DCMPLX(CPOTR_REAL(I,J),-CPOTR_IMG(I,J)) * FAC
        ENDDO
      END DO
C
      CALL CMAT1
CX          write(5432,'(7e15.6)') 
CX     1     rHO,VV(1,1),VV(2,2),VV(1,2)
      
      XX   = 1.D0 + EXP((RHO-R0)/A0)
      XX_W = 1.D0 + EXP((RHO-RW)/AW)
      VC0  = VC(RHO)

      DO I = 1 , MDIM
        VV(I,I) = VV(I,I) 
     1          + (-DCMPLX(V0/XX,W0/XX_W) + VC0 + CPOT1(I)) * FAC 
      END DO
C     
      RETURN
      END
C        
      SUBROUTINE RCWFNN(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEP)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
C .                  RCWFNN TO CALCULATE THE COULOMB REGULAR AND      .
C .                  IRREGULAR SOLUTIONS AND THEIR DERIVATIVES.       .
C .                  THIS PROGRAM IS THE MODIFIED VERSION OF THE      .
C .                  PROGRAM RCWFN FOR REAL(8) TYPE.         .
C .                  PUBLISHED IN: A.R. BARNETT, ET AL, COMPUTER      .
C .                  PHYSICS COMMUNICATIONS 8 (1974) 377--395.        .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT  REAL(8)   (A-H,O-Z)
      REAL(8)   K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(1),FCP(1),GC(1),GCP(1)
C *** COULOMB WAVEFUNCTIONS CALCULATED AT R = RHO BY THE
C *** CONTINUED-FRACTION METHOD OF STEED   MINL,MAXL ARE ACTUAL L-VALUES
C *** SEE BARNETT FENG STEED AND GOLDFARB COMPUTER PHYSICS COMMUN 1974
      PACE = STEP
      ACC  = ACCUR
      IF(PACE.LT.100.0D0) PACE = 100.0D0
      !the accuracy can be changed at here, Prof. Chuka
      !note-wenpw20190307
      IF(ACC.LT.1.0D-15.OR.ACC.GT.1.0D-6) ACC = 1.0D-10
      R    = RHO
      KTR  = 1
      LMAX = MAXL
      LMIN1= MINL + 1
      XLL1 = FLOAT(MINL*LMIN1)
      ETA2 = ETA*ETA
      TURN = ETA + SQRT(ETA2 + XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.0D-6) KTR = -1
      KTRP = KTR
      GO TO 2
1     R    = TURN
      TF   = F
      TFP  = FP
      LMAX = MINL
      KTRP = 1
2     ETAR = ETA*R
      RHO2 =   R*R
      PL   = FLOAT(LMAX + 1)
      PMX  = PL + 0.5D0
C *** CONTINUED FRACTION FOR FP(MAXL)/F(MAXL)  XL IS F  XLPRIME IS FP **
      FP  = ETA/PL + PL/R
      DK  = ETAR*2.0D0
      DEL = 0.0D0
      D   = 0.0D0
      F   = 1.0D0
      K   = (PL*PL - PL + ETAR)*(2.0D0*PL - 1.0D0)
      IF(PL*PL+PL+ETAR.NE.0.0D0) GO TO 3
      R   = R + 1.0D-6
      GO TO 2
3     H   = (PL*PL + ETA2)*(1.0D0 - PL*PL)*RHO2
      K   = K + DK + PL*PL*6.0D0
      D   =  1.0D0/(D*H + K)
      DEL =  DEL*(D*K - 1.0D0)
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL + 1.0D0)*D/PL
      PL  = PL + 1.0D0
      FP  = FP + DEL
      IF(D.LT.0.0D0) F = -F
      IF(PL.GT.20000.D0) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP  = F*FP
      IF( LMAX.EQ.MINL) GO TO 5
      FC (LMAX+1) = F
      FCP(LMAX+1) = FP
C *** DOWNWARD RECURSION TO MINL FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      L  = LMAX
      DO 4 LP  = LMIN1,LMAX
      PL = FLOAT(L)
      GC (L+1) = ETA/PL + PL/R
      GCP(L+1) = SQRT(ETA2 + PL*PL)/PL
      FC (L)   = (GC(L+1)*FC(L+1) + FCP(L+1))/GCP(L+1)
      FCP(L)   =  GC(L+1)*FC(L)   - GCP(L+1)*FC(L+1)
4     L  = L - 1
      F  = FC (LMIN1)
      FP = FCP(LMIN1)
5     IF(KTRP.EQ.-1) GO TO 1
C *** REPEAT FOR R = TURN IF RHO LT TURN
C *** NOW OBTAIN P + I.Q FOR MINL FROM CONTINUED FRACTION (32)
C *** REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING REAL(8)
      P  = 0.0D0
      Q  = R - ETA
      PL = 0.0D0
      AR = -(ETA2 + XLL1)
      AI =   ETA
      BR = 2.0D0*Q
      BI = 2.0D0
      WI = 2.0D0*ETA
      DR =   BR/(BR*BR + BI*BI)
      DI =  -BI/(BR*BR + BI*BI)
      DP = -(AR*DI + AI*DR)
      DQ =  (AR*DR - AI*DI)
6     P  =  P + DP
      Q  =  Q + DQ
      PL = PL + 2.0D0
      AR = AR + PL
      AI = AI + WI
      BI = BI + 2.0D0
      D  = AR*DR - AI*DI + BR
      DI = AI*DR + AR*DI + BI
      T  = 1.0D0/(D*D + DI*DI)
      DR =  T*D
      DI = -T*DI
      H  = BR*DR - BI*DI - 1.0D0
      K  = BI*DR + BR*DI
      T  = DP*H  - DQ*K
      DQ = DP*K  + DQ*H
      DP = T
      IF(PL.GT.46000.D0) GO TO 11
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P  = P/R
      Q  = Q/R
C *** SOLVE FOR FP,G,GP AND NORMALISE F  AT L=MINL
      G  = (FP - P*F)/Q
      GP = P*G - Q*F
      W  = 1.0D0/SQRT(FP*G - F*GP)
      G  = W*G
      GP = W*GP
      IF(KTR.EQ.1) GO TO 8
      F  = TF
      FP = TFP
      LMAX = MAXL
C *** RUNGE-KUTTA INTEGRATION OF G(MINL) AND GP(MINL) INWARDS FROM TURN
C ***             SEE FOX AND MAYERS 1968 PG 202
      I2=IDINT(PACE*ABS(RHO-TURN)+0.001D0)
      R3=1.D0/3.D0
      H=(RHO-TURN)/(1.0D0+DFLOAT(I2))
      H2=0.5D0*H
      ETAH=ETA*H
      H2LL=H2*XLL1
      S=(ETAH+H2LL/R)/R-H2
7     RH2= R + H2
      T  = (ETAH + H2LL/RH2)/RH2 - H2
      K1 = H2*GP
      M1 =  S*G
      K2 = H2*(GP + M1)
      M2 =  T*(G  + K1)
      K3 =  H*(GP + M2)
      M3 =  T*(G  + K2)
      M3 =     M3 + M3
      K4 = H2*(GP + M3)
      RH = R + H
      S  = (ETAH + H2LL/RH )/RH  - H2
      M4 =  S*(G + K3)
      G  = G  + (K1 + K2 + K2 + K3 + K4)*R3
      GP = GP + (M1 + M2 + M2 + M3 + M4)*R3
      R  = RH
      I2 = I2 - 1
      IF(ABS(GP).GT.1.0D300) GO TO 11
      IF(I2.GE.0) GO TO 7
      W  = 1.0D0/(FP*G - F*GP)
C *** UPWARD RECURSION FROM GC(MINL) AND GCP(MINL),STORED VALUES ARE R,S
C *** RENORMALISE FC,FCP FOR EACH L-VALUE
8     GC (LMIN1) = G
      GCP(LMIN1) = GP
      IF(LMAX.EQ.MINL) GO TO 10
      DO  9  L = LMIN1,LMAX
      T        = GC(L+1)
      GC (L+1) = (GC(L)*GC (L+1) - GCP(L))/GCP(L+1)
      GCP(L+1) =  GC(L)*GCP(L+1) - GC(L+1)*T
      FC (L+1) = W*FC (L+1)
9     FCP(L+1) = W*FCP(L+1)
      FC (LMIN1) = FC (LMIN1)*W
      FCP(LMIN1) = FCP(LMIN1)*W
      RETURN
10    FC (LMIN1) = W*F
      FCP(LMIN1) = W*FP
      RETURN
11    W  = 0.0D0
      G  = 0.0D0
      GP = 0.0D0
      GO TO 8
      END