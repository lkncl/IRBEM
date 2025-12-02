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
SUBROUTINE trace_drift_shell_opt(xx0,Lm,Lstar,leI0,Bposit,Bmin,posit,Nposit, t_resol, r_resol)
!
       USE fieldline_utils
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
       real*8, intent(inout) :: posit(3, NREBMAX, NDER_DEF*r_resol)
       real*8, intent(inout) :: Bposit(NREBMAX, NDER_DEF*r_resol)
       integer*4, intent(inout) :: Nposit(NDER_DEF*r_resol)

       INTEGER*4  Nreb,Nder,Ntet
       REAL*8     rr,rr2
       REAL*8     xx(3),x1(3),x2(3)
       REAL*8     xmin(3),Xsave(3)
       REAL*8     lati,longi,alti
       REAL*8     B(3),Bl,B0,B1,B3
       REAL*8     dsreb,smin

       INTEGER*4  I,J,Iflag,Iflag_I,ind,II,Ifail
       REAL*8     Lb
       REAL*8     leI,leI1
       REAL*8     tt
       REAL*8     tetl,tet1,dtet
       REAL*8     somme,BrR2
!
       REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
!
       REAL*8, ALLOCATABLE :: tet(:), phi(:)
!
       COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
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
       leI0 = baddata
       Bmin = baddata
!
       CALL COMPUTE_L_DIPOLE(xx0,Lb)
!
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) RETURN

       dsreb = Lb/Nreb

       CALL COMPUTE_FIELDLINE_FROM_MIR(xx0, dsreb, nrebmax, xmin, Bmin, leI0, Lm)
!      Ici, on a besoin de calculer l'invariant sur la ligne de champ
!      Mais aussi, le champ à l'équateur, et Lm

! derive
!
! calcul du point sur la ligne de champ a la surface de la terre du
! cote nord
       CALL COMPUTE_DRIFTSHELL_FOOTPOINTS(xx0, leI0, B0, dsreb, nrebmax, tet, phi, nder, ntet)
!
! trace la ligne de champ complete
!
       DO I = 1, Nder
              write(6,*) I, phi(I), tet(I)
              x1(1) = SIN(tet(I))*COS(phi(I))
              x1(2) = SIN(tet(I))*SIN(phi(I))
              x1(3) = COS(tet(I))

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
                     CALL sksyst(-dsreb,x1,x2,Bl,Ifail)
	              IF (Ifail.LT.0) THEN
                            RETURN
                     ENDIF
                     ind=ind+1
                     posit(1,ind,I)=x2(1)
                     posit(2,ind,I)=x2(2)
                     posit(3,ind,I)=x2(3)
                     Bposit(ind,I)=Bl
                     rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
                     IF (rr2.LT.1.) GOTO 205
                     x1(1) = x2(1)
                     x1(2) = x2(2)
                     x1(3) = x2(3)
              ENDDO
205           CONTINUE
              Nposit(I)=ind
              write(6,*) I, ind
              DO II = 1,3
                     x1(II)  = xsave(II)
              ENDDO
!
	       CALL CHAMP(x1,B,Bl,Ifail)
              IF (Ifail.LT.0) RETURN
              IF (Bl.LT.B0) RETURN
       ENDDO

!
! calcul de somme de BdS sur la calotte nord
!
       CALL INTEGRATE_FLUX_ON_POLAR_CAP(Nder, Ntet, tet, phi, somme)
       Lstar = 2.D0*pi*Bo/somme
!       IF (Lm.LT.0.D0) Lstar = -Lstar
!
       ! dealloc tableaux internes
       DEALLOCATE(tet, phi)
END SUBROUTINE
!