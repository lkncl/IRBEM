!***************************************************************************************************
! Copyright 2025 L. Nicolle, 2004 S. Bourdarie
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
!
SUBROUTINE trace_drift_shell_opt_constant_time(xx0,Lm,Lstar,leI0,Bposit,Bmin,posit,Nposit, t_resol, r_resol)
!
       IMPLICIT NONE
       INCLUDE 'variables.inc'
!
       INTEGER*4, PARAMETER :: NREB_DEF = 50
       INTEGER*4, PARAMETER :: NDER_DEF = 48
       INTEGER*4, PARAMETER :: NTET_DEF = 720
       INTEGER*4, PARAMETER :: NREBMAX = NREB_DEF*20


       REAL*8, intent(in) :: xx0(3)
       REAL*8, intent(inout) :: Lm, Lstar, leI0
       REAL*8, intent(out) :: Bmin
       INTEGER*4, intent(in) :: t_resol, r_resol

       ! Already allocated arrays from caller
       real*8, intent(inout) :: posit(3, 20*NREB_DEF, NDER_DEF*r_resol)
       real*8, intent(inout) :: Bposit(20*NREB_DEF, NDER_DEF*r_resol)
       integer*4, intent(inout) :: Nposit(NDER_DEF*r_resol)

       INTEGER*4  Nreb,Nder,Ntet
!
       INTEGER*4  k_ext,k_l,kint
       REAL*8     rr,rr2
       REAL*8     xx(3),x1(3),x2(3)
       REAL*8     xmin(3),Xsave(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,B1,B3
       REAL*8     dsreb,smin

       INTEGER*4  I,J,Iflag,Iflag_I,Ilflag,ind,II,Ifail
       REAL*8     Lb
       REAL*8     leI,leI1
       REAL*8     XY,YY
       REAL*8     aa,bb
       REAL*8     tt
       REAL*8     tetl,tet1,dtet
       REAL*8     somme,BrR2
!
       REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
!
       REAL*8, ALLOCATABLE :: tet(:), phi(:)
!
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
       COMMON /flag_L/Ilflag
       COMMON /magmod/k_ext,k_l,kint
       REAL*8     pi,rad
       common /rconst/rad,pi


       !===============================
       ! Résolutions dynamiques
       !===============================
       Nreb = NREB_DEF
       Nder = NDER_DEF * r_resol
       Ntet = NTET_DEF * t_resol


       ALLOCATE(tet(Nder))
       ALLOCATE(phi(Nder))

       dtet = pi/Ntet
!
       Lm = baddata
       Lstar = baddata
       leI0 = 0.D0
!
       CALL GEO_SM(xx0,xx)
       rr = SQRT(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
       tt = ACOS(xx(3)/rr)
       Lb  = rr/SIN(tt)/SIN(tt)
!       write(6,*)'L bete ',Lb
!
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
	  RETURN
       ENDIF
       Bmin = B0
!
       dsreb = Lb/Nreb
!
! calcul du sens du depart
!
       CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
	  RETURN
       ENDIF
       B1 = Bl
       CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          leI0 = baddata
          Bmin = baddata
          Ilflag = 0
	  RETURN
       ENDIF
       B3 = Bl
!
! attention cas equatorial
!
       IF(B1.GT.B0 .AND. B3.GT.B0)THEN
         aa = 0.5D0*(B3+B1-2.D0*B0)
         bb = 0.5D0*(B3-B1)
         smin = -0.5D0*bb/aa
         Bmin = B0 - aa*smin*smin
         leI0 = SQRT(1.D0-Bmin/B0)*2.D0*ABS(smin*dsreb)
         Lm = (Bo/Bmin)**(1.D0/3.D0)
!         write(6,*)'L McIlwain eq ',B0,leI0,Lm
         GOTO 100
       ENDIF
       IF (B3.GT.B1) THEN
         dsreb = -dsreb
       ENDIF
!
! calcul de la ligne de champ et de I
!
       Bmin = B0
       leI = 0.D0
       ind=0
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
!
!       write(6,*)dsreb
       DO J = 1,NREBMAX
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            leI0 = baddata
            Bmin = baddata
            Ilflag = 0
	    RETURN
	 ENDIF
!	 write(6,*)J,x1(1),x1(2),x1(3),Bl
         IF (Bl.LT.Bmin) THEN
           xmin(1) = x2(1)
           xmin(2) = x2(2)
           xmin(3) = x2(3)
           Bmin = Bl
         ENDIF
         IF (Bl.GT.B0) GOTO 20
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
         leI = leI + SQRT(1.D0-Bl/B0)
	 B1 = Bl
       ENDDO
20     CONTINUE
!	write(6,*)J,leI
!
       IF (J.GE.NREBMAX) THEN !open field line
        leI0 = baddata
        Bmin = baddata
        Ilflag = 0
	RETURN
       ENDIF
!
       leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)
       leI = leI*ABS(dsreb)
       leI0 = leI
!
! calcul de L Mc Ilwain (Mc Ilwain-Hilton)
!
       XY = leI*leI*leI*B0/Bo
       YY = 1.D0 + 1.35047D0*XY**(1.D0/3.D0) &
           + 0.465376D0*XY**(2.D0/3.D0) &
           + 0.0475455D0*XY
       Lm = (Bo*YY/B0)**(1.D0/3.D0)
!       write(6,*)'L McIlwain ',B0,leI0,Lm
!
! calcul de Bmin
!
       CALL sksyst(dsreb,xmin,x1,B3,Ifail)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
          Ilflag = 0
	  RETURN
       ENDIF
       CALL sksyst(-dsreb,xmin,x1,B1,Ifail)
       IF (Ifail.LT.0) THEN
          Bmin = baddata
          Ilflag = 0
	  RETURN
       ENDIF
       aa = 0.5D0*(B3+B1-2.D0*Bmin)
       bb = 0.5D0*(B3-B1)
       smin = -0.5D0*bb/aa
       Bmin = Bmin - aa*smin*smin
       IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.D0) THEN
        Lm = -Lm
       ENDIF
!
100    CONTINUE
       IF (ABS(Lm) .GT. 10.D0) THEN
        Ilflag = 0
        RETURN
       ENDIF
!
! derive
!
!
! calcul du point sur la ligne de champ a la surface de la terre du
! cote nord
!
       DO I = 1,3
         x1(I)  = xx0(I)
       ENDDO
       dsreb = ABS(dsreb)
       DO J = 1,NREBMAX
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
	    RETURN
	 ENDIF
	 rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
	 IF (rr.LT.1.D0) GOTO 102
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
102    CONTINUE
       smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
       smin = (1.D0-smin)/(rr-smin)
       CALL sksyst(smin*dsreb,x1,x2,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          Ilflag = 0
	  RETURN
       ENDIF
       rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
       tet(1) = ACOS(x2(3)/rr)
       phi(1) = ATAN2(x2(2),x2(1))
!       write(6,*)rr,tet(1),phi(1)
!
! trace la ligne de champ complete.
!
       DO I = 1,3
         xsave(I)  = x1(I)
       ENDDO
       x1(1) = rr*SIN(tet(1))*COS(phi(1))
       x1(2) = rr*SIN(tet(1))*SIN(phi(1))
       x1(3) = rr*COS(tet(1))
       ind=1
       Nposit(1)=ind
       posit(1,ind,1)=x1(1)
       posit(2,ind,1)=x1(2)
       posit(3,ind,1)=x1(3)
       Bposit(ind,1)=Bl
       DO J = 1,NREBMAX
         CALL sksyst(-dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
	    RETURN
	 ENDIF
	 ind=ind+1
	 posit(1,ind,1)=x2(1)
	 posit(2,ind,1)=x2(2)
	 posit(3,ind,1)=x2(3)
         Bposit(ind,1)=Bl
!	 write(6,*)J,x1(1),x1(2),x1(3),Bl
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
         IF (rr2.LT.1.) GOTO 201
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
       ENDDO
201    CONTINUE
       Nposit(1)=ind
        DO I = 1,3
         x1(I)  = xsave(I)
       ENDDO
       if (k_l .eq.0) then
          Ilflag = 0
	  RETURN
       endif
!
! et on tourne -> on se decale sur la surface en phi et on cherche teta
! pour avoir leI0 et B0 constants
!
       dsreb = -dsreb
       DO I = 2,Nder
        phi(I) = phi(I-1)+2.D0*pi/Nder
        Iflag_I = 0
!	write(6,*)tet(I)
	IF (Ilflag.EQ.0) THEN
         tetl = tet(I-1)
	 IF (I.GT.2) tetl = 2.D0*tet(I-1)-tet(I-2)
         tet1 = tetl
	ELSE
	 tetl = tet(I)
	 tet1 = tetl
	ENDIF
!	write(6,*)tetl
!	read(5,*)
	leI1 = -1.
!
107     CONTINUE
        x1(1) = SIN(tetl)*COS(phi(I))
        x1(2) = SIN(tetl)*SIN(phi(I))
        x1(3) = COS(tetl)
        Iflag = 0
        leI = 0.D0
!
        DO J = 1,NREBMAX
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
	    RETURN
	 ENDIF
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
	 IF (Bl.LT.B0) THEN
	  IF (Iflag .EQ. 0) THEN
	    CALL CHAMP(x1,B,B1,Ifail)
	    IF (Ifail.LT.0) THEN
               Ilflag = 0
	       RETURN
	    ENDIF
	    leI = 0.5D0*SQRT(1.D0-Bl/B0)*(1.D0+(Bl-B0)/(Bl-B1))
	    Iflag = 1
	  ELSE
	    leI = leI+SQRT(1.D0-Bl/B0)
	  ENDIF
	 ENDIF
         IF (Bl.GT.B0 .AND. Iflag.EQ.1) GOTO 103
	 IF (rr2.LT.1.D0) GOTO 103
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
!	 write(6,*)J,Bl,B0,leI*ABS(dsreb)
        ENDDO
103     CONTINUE
        IF (rr2.LT.1.D0) THEN
         leI = 0.D0
	ENDIF
        IF (J.LT.NREBMAX .AND. rr2.GE.1.D0) THEN
            CALL CHAMP(x1,B,B1,Ifail)
	    IF (Ifail.LT.0) THEN
               Ilflag = 0
	       RETURN
	    ENDIF
            leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)
            leI = leI*ABS(dsreb)
	ENDIF
!	write(6,*)I,tetl,leI,leI0,J
!
        IF (Iflag_I .EQ.0) THEN
	 IF (J.GE.NREBMAX) THEN
	  tetl = tetl-dtet
         ELSE
	  tetl = tetl+dtet
	 ENDIF
	 leI1 = leI
	 tet1 = tetl
	 Iflag_I = 1
	 GOTO 107
	ENDIF
	IF ((leI-leI0)*(leI1-leI0) .LT. 0.D0) GOTO 108
	leI1 = leI
	tet1 = tetl
	IF (leI.LT.leI0) THEN
	 tetl = tetl-dtet
	ElSE
	 tetl = tetl+dtet
	ENDIF
	IF (tetl.GT.pi .OR. tetl.LT.0.D0) GOTO 108
	GOTO 107
108     CONTINUE
        tet(I) = 0.5D0*(tetl+tet1)
!	read(5,*)
	IF (J.GE.NREBMAX .AND. leI.GT.0.) THEN
	 Ilflag = 0
	 RETURN
	ENDIF
!
        x1(1) = SIN(tet(I))*COS(phi(I))
        x1(2) = SIN(tet(I))*SIN(phi(I))
        x1(3) = COS(tet(I))
!
! trace la ligne de champ complete
!
        DO II = 1,3
         xsave(II)  = x1(II)
        ENDDO
	ind=1
        Nposit(I)=ind
	posit(1,ind,I)=x1(1)
	posit(2,ind,I)=x1(2)
	posit(3,ind,I)=x1(3)
        Bposit(ind,I)=Bl
        DO J = 1,NREBMAX
         CALL sksyst(dsreb,x1,x2,Bl,Ifail)
	 IF (Ifail.LT.0) THEN
            Ilflag = 0
	    RETURN
	 ENDIF
	 ind=ind+1
	 posit(1,ind,I)=x2(1)
	 posit(2,ind,I)=x2(2)
	 posit(3,ind,I)=x2(3)
         Bposit(ind,I)=Bl
!	 write(6,*)J,x1(1),x1(2),x1(3),Bl
         rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
         IF (rr2.LT.1.) GOTO 205
	 x1(1) = x2(1)
	 x1(2) = x2(2)
	 x1(3) = x2(3)
        ENDDO
205     CONTINUE
        Nposit(I)=ind
        DO II = 1,3
         x1(II)  = xsave(II)
        ENDDO
!
	CALL CHAMP(x1,B,Bl,Ifail)
	IF (Ifail.LT.0) THEN
           Ilflag = 0
	   RETURN
	ENDIF
	IF (Bl.LT.B0) THEN
	 Ilflag = 0
	 Lstar = baddata
	 RETURN
	ENDIF
       ENDDO
!       write(6,*)(tet(I),I=1,Nder)
!       read(5,*)
!
! calcul de somme de BdS sur la calotte nord
!
       x1(1) = 0.D0
       x1(2) = 0.D0
       x1(3) = 1.D0
       CALL CHAMP(x1,B,Bl,Ifail)
       IF (Ifail.LT.0) THEN
          Ilflag = 0
	  RETURN
       ENDIF
       BrR2 = abs((x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3))) ! phi integrates B dot dA, or Br*R^2dphi*dtheta, R=1
       somme = BrR2*pi*dtet*dtet/4.D0
       DO I = 1,Nder
         tetl = 0.D0
         DO J = 1,Ntet
	  tetl = tetl+dtet
	  IF (tetl .GT. tet(I)) GOTO 111
          x1(1) = SIN(tetl)*COS(phi(I))
          x1(2) = SIN(tetl)*SIN(phi(I))
          x1(3) = COS(tetl)
          CALL CHAMP(x1,B,Bl,Ifail)
          IF (Ifail.LT.0) THEN
             Ilflag = 0
	     RETURN
          ENDIF
          BrR2 = abs((x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3))) ! phi integrates B dot dA, or Br*R^2dphi*dtheta, R=1
	  somme = somme+BrR2*SIN(tetl)*dtet*2.D0*pi/Nder
	 ENDDO
111      CONTINUE
       ENDDO
       Lstar = 2.D0*pi*Bo/somme
!       IF (Lm.LT.0.D0) Lstar = -Lstar
       Ilflag = 0
!

       ! dealloc tableaux internes
       DEALLOCATE(tet, phi)
END SUBROUTINE
!