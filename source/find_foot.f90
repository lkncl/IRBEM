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
! Modified from find_bm and find_mirror_point by TPO (July 2008)
! Modified to prevent overshooting into the other hemisphere ACK (Nov 2016)
!
! Routine to find foot point of field line at specified altitude and hemi
! finds foot point at specified altitude to within 1 km
! --------------------------------------------------------------------
!
        SUBROUTINE find_foot_point1(kext,options,sysaxes,iyearsat,&
           idoy,UT,xIN1,xIN2,xIN3,stop_alt,hemi_flag,maginput,&
           XFOOT,BFOOT,BFOOTMAG)
!
!       INPUTS have the usual meaning, except:
!      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
!      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
!      0 - same Hemisphere as start point
!      +1 - Northern Hemisphere
!      -1 - Southern Hemisphere
!      2 - opposite Hemisphere as start point
!
!      OUTPUTS
!      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
!      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
!      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)

	IMPLICIT NONE
	INCLUDE 'variables.inc'
!
! declare inputs
        INTEGER*4    kext,k_ext,k_l,kint,options(5)
        INTEGER*4    sysaxes
	INTEGER*4    iyearsat
	integer*4    idoy
	real*8     UT
	real*8     xIN1,xIN2,xIN3
	real*8     stop_alt
        INTEGER*4  hemi_flag
	real*8     maginput(25)
!
! Declare internal variables
	INTEGER*4    isat,iyear,Iint
        INTEGER*4    Ndays,activ,Ifail
	REAL*8     psi,mlon,tilt
	REAL*8     xGEO(3)
	real*8     alti,lati,longi
	real*8     BxGEO(3)
!
! Declare output variables	
        REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG
!
	COMMON /magmod/k_ext,k_l,kint
      integer*4 int_field_select, ext_field_select
!

        if ((stop_alt.lt.0).or.(stop_alt.ge.6378.0*500.0)) then
           goto 999 ! fail, stop_alt out of range 0 to 500 Re
        endif
	kint = int_field_select ( options(5) )
	k_ext = ext_field_select ( kext )
!
        CALL INITIZE
	
	call init_fields ( kint, iyearsat, idoy, ut, options(2) )
	
	call get_coordinates ( sysaxes, xIN1, xIN2, xIN3, &
       alti, lati, longi, xGEO )
	    
	call set_magfield_inputs ( k_ext, maginput, ifail )
 
	if ( ifail.lt.0 ) goto 999

        if (k_ext .eq. 13 .or. k_ext .eq. 14) then 
            call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
            call INIT_TS07D_TLPR
	        if ( ifail.lt.0 ) goto 999
        end if
!
!
        CALL find_foot(lati,longi,alti,stop_alt,&
               hemi_flag,XFOOT,BFOOT,BFOOTMAG)
        return
 999    continue ! bad data
        XFOOT(1) = baddata
        XFOOT(2) = baddata
        XFOOT(3) = baddata
        BFOOT(1) = baddata
        BFOOT(2) = baddata
        BFOOT(3) = baddata
        BFOOTMAG = baddata
	END
      
       SUBROUTINE find_foot(lati,longi,alti,stop_alt,hemi_flag,&     
         XFOOT,BFOOT,BFOOTMAG)
!
!      inputs: 
!      REAL*8 lati - geodetic latitude of start point (gdz), degrees
!      REAL*8 longi - geodetic longitude of start point  (gdz), degrees
!      REAL*8 alti - geodetic altitude of start point  (gdz), km
!      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
!      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
!      0 - same Hemisphere as start point
!      +1 - Northern Hemisphere
!      -1 - Southern Hemisphere
!      2 - opposite Hemisphere as start point
!
!      outputs:
!      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
!      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
!      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)
       IMPLICIT NONE

       INTEGER*4  k_ext,k_l,kint,Ifail
       REAL*8     xx0(3)
       REAL*8     lati,longi,alti
       REAL*8     stop_alt
       INTEGER*4  hemi_flag
       REAL*8     XFOOT(3),BFOOT(3),BFOOTMAG

       CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
!
       call find_foot_opt ( xx0,stop_alt,hemi_flag,&
            XFOOT,BFOOT,BFOOTMAG)

       RETURN
       END

       SUBROUTINE find_foot_opt(xx0,stop_alt,hemi_flag,&
            XFOOT,BFOOT,BFOOTMAG)
!
!      inputs: 
!      REAL*8 xx0(3) - GEO cartesian coordinates
!      REAL*8 stop_alt - geodetic altitude of desired foot point (gdz), km
!      integer*4 hemi_flag - hemishere flag, specifies hemisphere of foot point
!      0 - same Hemisphere as start point
!      +1 - Northern Hemisphere
!      -1 - Southern Hemisphere
!      2 - opposite Hemisphere as start point
!
!      outputs:
!      REAL*8 XFOOT(3) - GDZ position of foot point (alt, lat, lon)
!      REAL*8 BFOOT(3) - Magnetic field at foot point (nT, GEO)
!      REAL*8 BFOOTMAG - Magnetic field at foot point (nT)

       USE fieldline_utils
       IMPLICIT NONE
       INCLUDE 'variables.inc'
!
       REAL(8), INTENT(IN) :: stop_alt
       INTEGER(4), INTENT(IN) :: hemi_flag
       REAL(8), INTENT(OUT) :: xx0(3), XFOOT(3),BFOOT(3),BFOOTMAG
       INTEGER(4)  Nreb
       PARAMETER (Nreb = 50)

       INTEGER(4) :: stop_type, nrebmax
       REAL(8) :: dsreb, Lb

       CALL COMPUTE_L_DIPOLE(xx0, Lb)

       dsreb = Lb/(Nreb*1.d0) ! step size
       nrebmax = 10*Nreb ! legacy 500 steps
       stop_type = 1
       CALL FIND_FIELDLINE_FOOT(xx0, dsreb, nrebmax, stop_type, stop_alt, &
         hemi_flag, xfoot, bfoot, bfootmag)

       END SUBROUTINE
