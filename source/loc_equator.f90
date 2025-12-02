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
!
       SUBROUTINE loc_equator(lati,longi,alti,Bmin,posit)
!
       IMPLICIT NONE
       REAL*8     lati,longi,alti
       REAL*8     Bmin
       REAL*8     posit(3)
       REAL*8     xx0(3)
       
       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
       
       call loc_equator_opt(xx0,Bmin,posit)

       RETURN
       END

       SUBROUTINE loc_equator_opt(xx0,Bmin,posit)
!
       USE fieldline_utils
       IMPLICIT NONE
       INCLUDE 'variables.inc'
!
       INTEGER*4  Nreb
       PARAMETER (Nreb = 50)
!
       INTEGER*4  Nrebmax
       REAL*8     xx0(3)
       REAL*8     B(3),B0,Bmin
       REAL*8     dsreb

       INTEGER*4  Ifail
       REAL*8     Lb
!
!
       REAL*8     posit(3)
!
!
       Nrebmax = 10*Nreb
!
       CALL COMPUTE_L_DIPOLE(xx0, Lb)
!
       CALL CHAMP(xx0,B,B0,Ifail)
       IF (Ifail.LT.0) THEN
              Bmin=baddata
              posit(1) = baddata
              posit(2) = baddata
              posit(3) = baddata
              RETURN
       ENDIF
       Bmin = B0

       dsreb = Lb/Nreb
       CALL FIND_MAGEQUATOR(xx0, B0, dsreb, nrebmax, posit, Bmin)

       END
       
