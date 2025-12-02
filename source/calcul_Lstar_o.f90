!***************************************************************************************************
! Copyright 2003, 2004, D. Boscher, S. Bourdarie
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
! Boscher modifie pour la nieme fois le 4Feb2004
! A. C. Kellerman added code to exit search along NH FL when Nrebmax is reached Jan, 2017
!
!
       SUBROUTINE calcul_Lstar(t_resol,r_resol,  &
                lati,longi,alti,Lm,Lstar,leI0,B0,Bmin)
!
       IMPLICIT NONE
!
       INTEGER*4  t_resol,r_resol
       REAL*8     xx0(3)
       REAL*8     lati,longi,alti
       REAL*8     B0,Bmin
       REAL*8     Lm,Lstar
       REAL*8     leI0
!
       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
!
       call calcul_Lstar_opt(t_resol,r_resol,    &
                xx0,Lm,Lstar,leI0,B0,Bmin)

       RETURN
       END
       SUBROUTINE calcul_Lstar_opt(t_resol,r_resol,     &
                xx0,Lm,Lstar,leI0,B0,Bmin)
!
       USE fieldline_utils
       IMPLICIT NONE
       INCLUDE 'variables.inc'
!
       INTEGER*4  Nreb_def,Nder_def,Ntet_def
       PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
!
       INTEGER*4  Nder,Nreb,Ntet
       INTEGER*4  k_ext,k_l,kint,n_resol,t_resol,r_resol
       INTEGER*4  Nrebmax
       REAL*8     rr,rr2
       REAL*8     xx0(3),xx(3),x1(3),x2(3)
       REAL*8     xmin(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,Bmin,B1,B3
       REAL*8     dsreb,smin

       INTEGER*4  I,J,Iflag,Iflag_I,Ilflag,Ifail,Iflag3
       REAL*8     Lm,Lstar,Lb
       REAL*8     leI,leI0,leI1
       REAL*8     XY,YY
       REAL*8     aa,bb
!
       REAL*8     tt
       REAL*8     tet(Nder_def*r_resol),phi(Nder_def*r_resol)
       REAL*8     tetl,tet1,dtet, phil
       REAL*8     somme,BrR2
       REAL*8     Bmin3,x3(3),xmin3(3),Bl3,rr32,x13(3),Bl13,ds3
!
       REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
!
       INTEGER*4     FLAG_IN_EARTH
!
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /flag_L/Ilflag
       COMMON /magmod/k_ext,k_l,kint

       REAL*8     pi,rad
       common /rconst/rad,pi
!
       INTEGER*4 depth

       tet(:) = 0
       FLAG_IN_EARTH = 0
!
       Nder=Nder_def*r_resol
       Nreb=Nreb_def
       Ntet=Ntet_def*t_resol
       dtet = pi/Ntet
!
       Nrebmax = 20*Nreb
!
       Lm = baddata
       Lstar = baddata
       leI0 = baddata

       CALL COMPUTE_L_DIPOLE(xx0, Lb)

       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
	   B0=baddata
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
	   RETURN
       ENDIF
       Bmin = B0
!
       dsreb = Lb/Nreb

       CALL COMPUTE_FIELDLINE_FROM_MIR(xx0, dsreb, nrebmax, xmin, Bmin, leI0, Lm)
!      Ici, on a besoin de calculer l'invariant sur la ligne de champ
!      Mais aussi, le champ à l'équateur, et Lm

! derive
!
       ! write(6,*) xx0, b0, dsreb, nrebmax, nder, ntet
       CALL COMPUTE_DRIFTSHELL_FOOTPOINTS(xx0, leI0, B0, dsreb, nrebmax, tet, phi, nder, ntet)
       ! write(6,*) "phis", phi
       ! write(6,*) "tets", tet
! calcul de somme de BdS sur la calotte nord
! "calculation of sum of BdS on the northern cap"
!
       CALL INTEGRATE_FLUX_ON_POLAR_CAP(Nder, Ntet, tet, phi, somme)

       if (k_l .eq.1) Lstar = 2.D0*pi*Bo/somme
       if (k_l .eq.2) Lstar = somme  ! Phi and not Lstar
       IF (Lm.LT.0.D0) Lstar =-Lstar
       !IF(Lm.GE.0.D0 .and. flag_in_earth.eq.1) Lstar = -Lstar
       Ilflag = 1
!
       END
!
!***********************************************************************
!* RUNGE-KUTTA d'ordre 4
!***********************************************************************
        SUBROUTINE sksyst(h,xx,x2,Bl,Ifail)
!
        IMPLICIT NONE
!
        INTEGER*4  Ifail
        REAL*8 xx(3),x2(3)
        REAL*8 B(3),Bl
        REAL*8 h
        REAL*8 xwrk(4,3)
!
!-----------------------------------------------------------------------
!
!        write(6,*)'sksyst'
!        write(6,*)xx(1),xx(2),xx(3),h
        CALL CHAMP(xx,B,Bl,Ifail)
!	write(6,*)xx,B,Bl,Ifail
	IF (Ifail.LT.0) RETURN
!        write(6,*)'b',B(1),B(2),B(3),Bl
        xwrk(1,1) = h*B(1)/Bl
        xwrk(1,2) = h*B(2)/Bl
        xwrk(1,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(1,1)/2.D0
        x2(2) = xx(2)+xwrk(1,2)/2.D0
        x2(3) = xx(3)+xwrk(1,3)/2.D0
!        write(6,*)x2(1),x2(2),x2(3),Bl
!
        CALL CHAMP(x2,B,Bl,Ifail)
	IF (Ifail.LT.0) RETURN
        xwrk(2,1) = h*B(1)/Bl
        xwrk(2,2) = h*B(2)/Bl
        xwrk(2,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(2,1)/2.D0
        x2(2) = xx(2)+xwrk(2,2)/2.D0
        x2(3) = xx(3)+xwrk(2,3)/2.D0
!        write(6,*)x2(1),x2(2),x2(3),Bl
!
        CALL CHAMP(x2,B,Bl,Ifail)
	IF (Ifail.LT.0) RETURN
        xwrk(3,1) = h*B(1)/Bl
        xwrk(3,2) = h*B(2)/Bl
        xwrk(3,3) = h*B(3)/Bl
        x2(1) = xx(1)+xwrk(3,1)
        x2(2) = xx(2)+xwrk(3,2)
        x2(3) = xx(3)+xwrk(3,3)
!        write(6,*)x2(1),x2(2),x2(3),Bl
!
        CALL CHAMP(x2,B,Bl,Ifail)
	IF (Ifail.LT.0) RETURN
        xwrk(4,1) = h*B(1)/Bl
        xwrk(4,2) = h*B(2)/Bl
        xwrk(4,3) = h*B(3)/Bl
!
!-----------------------------------------------------------------------
!
        x2(1) = xx(1)+(   xwrk(1,1)+2.D0*xwrk(2,1)      &
                       + 2.D0*xwrk(3,1)+   xwrk(4,1))/6.D0
        x2(2) = xx(2)+(   xwrk(1,2)+2.D0*xwrk(2,2)      &
                       + 2.D0*xwrk(3,2)+   xwrk(4,2))/6.D0
        x2(3) = xx(3)+(   xwrk(1,3)+2.D0*xwrk(2,3)      &
                       + 2.D0*xwrk(3,3)+   xwrk(4,3))/6.D0
        CALL CHAMP(x2,B,Bl,Ifail)
	IF (Ifail.LT.0) RETURN
!        write(6,*)x2(1),x2(2),x2(3),Bl
!        read(5,*)
!
        RETURN
        END
        SUBROUTINE sksyst2(h,xx,x2,Bl,Ifail)
!
        IMPLICIT NONE
!
        INTEGER*4  II,JJ,Ifail
        REAL*8 xx(3),x2(3)
        REAL*8 B(3),Bl
        REAL*8 h, c1, c2, c3
        REAL*8 xw, xf(3)
!
!-----------------------------------------------------------------------
!
!        write(6,*)'sksyst'
!        write(6,*)xx(1),xx(2),xx(3),h
        CALL CHAMP(XX,B,Bl,Ifail)
!        write(6,*)xx,B,Bl,Ifail
        IF (Ifail.LT.0) RETURN
	C1 = 1.D0
	C2 = 0.5D0
	DO JJ=1,3
          C3 = h/Bl
	  IF (JJ.eq.3) C2 = 1.D0
          DO II=1,3
!            write(6,*)'b',B(1),B(2),B(3),Bl
            xw = C3*B(II)
            x2(II) = xx(II)+xw*C2
            XF(II) = XF(II)+xw*C1
          ENDDO
          CALL CHAMP(X2,B,Bl,Ifail)
!          write(6,*)xx,B,Bl,Ifail
          IF (Ifail.LT.0) RETURN
	  C1 = 2.D0
        ENDDO
        C3 = h/Bl
        DO II=1,3
!          write(6,*)'b',B(1),B(2),B(3),Bl
          xw = C3*B(II)
          x2(II) = xx(II) + xf(ii)/6.D0
        ENDDO
        CALL CHAMP(X2,B,Bl,Ifail)
!        write(6,*)xx,B,Bl,Ifail
        IF (Ifail.LT.0) RETURN
        END


