!***************************************************************************************************
! Copyright  2004 S. Bourdarie
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
!

!     Legacy wrapper w/ R0=1.0
       SUBROUTINE field_line_tracing(lati,longi,alti,&
              Lm,leI0,Bposit,Bmin,posit,Nposit)
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     lati,longi,alti,R0
       REAL*8     Lm,leI0,Bmin,xx0(3)
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
!
       call field_line_tracing_opt2(xx0,1.D0, & 
              Lm,leI0,Bposit,Bmin,posit,Nposit)

       RETURN
       END
       SUBROUTINE field_line_tracing2(lati,longi,alti, &
              R0,Lm,leI0,Bposit,Bmin,posit,Nposit)
!
!    modified from field_line_tracing to add R0 parameter: radius (Re) of
!    reference surface
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     lati,longi,alti,R0
       REAL*8     Lm,leI0,Bmin,xx0(3)
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
!
       call field_line_tracing_opt2(xx0,R0,& 
              Lm,leI0,Bposit,Bmin,posit,Nposit)

       RETURN
       END
       
       SUBROUTINE field_line_tracing_opt(xx0,&
              Lm,leI0,Bposit,Bmin,posit,Nposit)
       IMPLICIT NONE
       INTEGER*4  Nreb,Ntet
       PARAMETER (Nreb = 150, Ntet = 720)
       INTEGER*4  Nposit
       REAL*8     xx0(3),R0
       REAL*8     Lm,leI0,Bmin
       REAL*8     posit(3,20*Nreb),Bposit(20*Nreb)

       call field_line_tracing_opt2(xx0,1.D0,& 
              Lm,leI0,Bposit,Bmin,posit,Nposit)
       RETURN
       END

       SUBROUTINE field_line_tracing_opt2(xx0,&
              R0,Lm,leI0,Bposit,Bmin,posit,Nposit)
       !    modified from field_line_tracing to add R0 parameter: radius (Re) of
       !    reference surface
              IMPLICIT NONE
              INTEGER*4  Nreb,Ntet,Nrebmax
              PARAMETER (Nreb = 150, Ntet = 720, Nrebmax = Nreb*20)

              REAL(8), INTENT(IN) :: xx0(3), R0
              REAL(8), INTENT(OUT) :: Lm, leI0, Bmin
              REAL(8), INTENT(OUT) :: posit(3,Nrebmax),Bposit(Nrebmax)
              INTEGER(4), INTENT(OUT) :: Nposit

              CALL field_line_tracing_flex(xx0, &
                     R0,Lm,leI0,Bposit,Bmin,posit,Nposit, &
                     Nreb)
       END SUBROUTINE

       SUBROUTINE field_line_tracing_flex(xx0,&
              R0,Lm,leI0,Bposit,Bmin,posit,Nposit, &
              Nreb)

              USE fieldline_utils
              IMPLICIT NONE
              INCLUDE 'variables.inc'
       !
              INTEGER(4), INTENT(IN) :: Nreb
              INTEGER(4) :: Nrebmax
              REAL*8     rr2
              REAL*8     xx0(3),x1(3),x2(3)
              REAL*8     xmin(3)
              REAL*8     B(3),Bl,B0,Bmin
              REAL*8     dsreb

              INTEGER*4  I, J,ind,Ifail
              INTEGER*4  Nposit
              REAL*8     Lm,Lb
              REAL*8     leI0
       !
              INTEGER(4) :: stop_type, hemi_flag
              REAL(8) :: stop_value, bfoot(3), bfootmag
       !
       !
              REAL*8     posit(3,Nrebmax),Bposit(Nrebmax)
              real*8     R0,R02 ! R0^2
       !
       !
       !
              R02 = R0*R0
       !
              Nrebmax = 20*Nreb
       !
              Lm = baddata
              leI0 = baddata
              Bmin = baddata
       !

              CALL COMPUTE_L_DIPOLE(xx0, Lb)
       !
              CALL CHAMP(xx0,B,B0,Ifail)
              IF (Ifail.LT.0) RETURN

              dsreb = Lb/Nreb

              CALL COMPUTE_FIELDLINE_FROM_MIR(xx0, dsreb, nrebmax, xmin, Bmin, leI0, Lm)
              
              stop_type = 0 ! Radius GEO stop criterion
              stop_value = 1.0D0 ! Re = 1 stop condition
              hemi_flag = 1 ! north
              CALL FIND_FIELDLINE_FOOT(xx0, dsreb, nrebmax, stop_type, stop_value, &
              hemi_flag, x1, bfoot, bfootmag)

       ! trace la ligne de champ complete.
       !
              !tracing

              ind=1
              Nposit=ind
              posit(1,ind)=x1(1)
              posit(2,ind)=x1(2)
              posit(3,ind)=x1(3)
              Bposit(ind)=Bl
              DO J = 1,Nrebmax-1 ! this is the corrected version
                                   !, gives 3000 values, A. Kellerman
                     CALL sksyst(-dsreb,x1,x2,Bl,Ifail)
                     IF (Ifail.LT.0) RETURN
                     ind=ind+1
                     posit(1,ind)=x2(1)
                     posit(2,ind)=x2(2)
                     posit(3,ind)=x2(3)
                     Bposit(ind)=Bl
                     rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
                     IF (rr2.LT.R02) GOTO 201
                     x1(1) = x2(1)
                     x1(2) = x2(2)
                     x1(3) = x2(3)
              ENDDO
201           CONTINUE
              Nposit=ind
       END
!
