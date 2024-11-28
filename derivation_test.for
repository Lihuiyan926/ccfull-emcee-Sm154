
!   second derivative of the excitation function....jason's program appended
!    here
!  modified by Lihuiyan 2024-3-22
!  theta =180 degree, ei=Eeff
!   **************************************************************************
      program DERIVATIVE
      IMPLICIT REAL*8 (A-H,O-Z)


      DIMENSION e(5000) , sig(5000), ei(5000), elab(5000), ecm(5000)
      OPEN(10,FILE='Ccfull-sc2.inp',STATUS='old')
      OPEN (8,FILE='AngDis.dat',STATUS='old')

      open(16,file='db_test.out',status='unknown')
      !barrier distribution DB=d^2(E*sigma)/d E^2
      !open(17,file='sfact.out',status='unknown')
      !s(e)= E* sgima(E) * exp(2*pi*eta), eta=z1*z2*e^2/hbar/v : sommerfiled parameter
c      open(18,file='slope.out',status='unknown')
      !L(e)=d(E* sgima(E)) / dE /(E* sgima(E))

      write(16,'(a10, 3a15)')"e","DB","S","L"


      Stp=2. !step size
      !read(10,*)
      READ(10,*)AP,ZP,AT,ZT

      HBAR=197.329D0
      PI=3.141592653D0
      ANMASS=938.D0
      RMASS=AP*AT/(AP+AT)*ANMASS

       READ (8,*)
      !calculate centroid barrier
      sum1 = 0.
      sum2 = 0.
      accur = .05
      i = 0
 910  i = i + 1
c      READ (8,*,END=977) e(i) , sig(i) , err , int
      READ (8,*,END=977) ei(i) ,a,a,a,a,a,a,a, sig(i)  !sig is only for kantbp calculation
      !theta = 180*3.14159/180
      !theta = 175*PI/180+dasin(AP*dsin(175*PI/180)/AT)
      theta = 180*PI/180
      !csc1 = 1.d0/dsin(0.5d0*theta)
      !e(i) = ei(i)*(csc1-1.d0)/(csc1+1.d0)
      ! ecm(i) = 2.d0*ei(i)*dsin(0.5d0*theta)/(dsin(0.5d0*theta)+1.d0)
      !e(i) = ecm(i)*(ap+at)/at
      e(i) = ei(i)

      GOTO 910
 977  n = i - 1

      DO i = 1 , n
         l = 0
         DO j = i , n
            diff = ABS(e(j)-e(i))

            IF ( diff.LE.(Stp+accur) .AND. diff.GE.(Stp-accur) ) l = j
            IF ( l.NE.0  .AND. i.NE.l  ) THEN

               e1 = e(i)
               e3 = e(l)
               e2=(e3+e1)/2

               s1 = sig(i)
               s3 = sig(l)

               !barrier distribution
               !three point formula
               deriv =(s3 - s1) / (e3 - e1)
               deriv_1= (1./(e3-e1))*(e3*s3-e1*s1)


               write(16,'(F10.2, 2E15.4)')e2,-deriv,deriv_1
               ! write(*,'(F10.2, 2E15.4)')e2,-deriv ,deriv_1

               l = 0
               k = 0
            ENDIF


         ENDDO
      ENDDO

      END
