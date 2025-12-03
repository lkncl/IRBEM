!***************************************************************************************************
! Copyright 2004 S. Bourdarie
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
! S. Bourdarie (June 2004)
! Modified S./ Bourdarie (July 2004)
!
! Routine to find mirror point of a trapped particle
!
SUBROUTINE find_bm(lati,longi,alti,alpha,Bposit,Bmir,xmin)

    IMPLICIT NONE
    REAL*8     xx0(3),xmin(3)
    REAL*8     lati,longi,alti
    REAL*8     Bposit,Bmir
    REAL*8     alpha

    CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
    call find_bm_nalpha (xx0,1,alpha,Bposit,Bmir,xmin)

    RETURN
END SUBROUTINE

! alpha : degrees
SUBROUTINE find_bm_nalpha(xx0,nalp,alpha,Bposit,Bmir,xmin)
    
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: Nalp
    REAL(8), INTENT(IN) :: xx0(3), alpha(Nalp)
    REAL(8), INTENT(OUT) :: Bposit, Bmir(Nalp), xmin(3, Nalp)

    INTEGER*4  Nreb
    PARAMETER (Nreb = 50)

    CALL find_bm_nalpha_flex(xx0,nreb,nalp,alpha,Bposit,Bmir,xmin)

END SUBROUTINE


! alpha : degrees
SUBROUTINE find_bm_nalpha_flex(xx0,nreb,nalp,alpha,Bposit,Bmir,xmin)

    USE fieldline_utils
    IMPLICIT NONE
    INCLUDE 'variables.inc'

    INTEGER(4), INTENT(IN) :: Nreb, Nalp
    REAL(8), INTENT(IN) :: xx0(3), alpha(Nalp)
    REAL(8), INTENT(OUT) :: Bposit, Bmir(Nalp), xmin(3, Nalp)

    INTEGER(4) :: Ifail, I, ii, Nrebmax
    REAL(8) :: x1(3),x2(3)
    REAL(8) :: B(3),Bl,B1,B3,xmir(3), bmirvect(3), bmirmag
    REAL(8) :: Lb, dsreb

    Nrebmax = 20*Nreb

    ! initialization
    Bposit=baddata
    Bmir = baddata
    xmin = baddata

    CALL COMPUTE_L_DIPOLE(xx0, Lb)

    CALL CHAMP(xx0,B,Bl,Ifail)
    IF (Ifail.LT.0) RETURN
    Bposit = Bl

    dsreb = Lb/Nreb
    
    do ii=1,nalp
    
        IF( alpha(ii).eq.90.d0 ) then
            Bmir(ii) = 0.d0 ! why zero ? should be B at xx0 -> Bl
            xmin(1,ii) = xx0(1)
            xmin(2,ii) = xx0(2)
            xmin(3,ii) = xx0(3)
            goto 100
	    ELSE
            ! compute direction of departure
            CALL sksyst(-dsreb,xx0,x1,Bl,Ifail)
            IF (Ifail.LT.0) GOTO 100
            B1 = Bl

            CALL sksyst(dsreb,xx0,x2,Bl,Ifail)
            IF (Ifail.LT.0) GOTO 100
            B3 = Bl

            IF (B1.GT.B3) dsreb = -dsreb

            ! computation of mirror point
            CALL FIND_MIRROR_POINT(xx0, alpha(ii), dsreb, nrebmax, xmir, bmirvect, bmirmag)
            DO I = 1,3
                xmin(I, ii) = xmir(I)
            ENDDO
            Bmir(ii) = bmirmag
        ENDIF

100 CONTINUE
    ENDDO
END SUBROUTINE