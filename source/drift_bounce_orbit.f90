!***************************************************************************************************
! Copyright 2010 T.P. O'Brien
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
!-----------------------------------------------------------------------------
      SUBROUTINE drift_bounce_orbit1(kext,options,sysaxes,&
           iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,&
           Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind)
      ! inputs
      INTEGER*4 kext, options(5),sysaxes,iyearsat,idoy
      REAL*8 UT,xIN1,xIN2,xIN3,alpha,maginput(25)
      ! outputs
      REAL*8 Lm,Lstar
      REAL*8 BLOCAL(1000,25),BMIN,BMIR,XJ,posit(3,1000,25)
      INTEGER*4 ind(25)
      ! locals
      REAL*8 R0,hmin,hmin_lon
      R0 = 1.D0 ! stop search at surface of Earth
      call drift_bounce_orbit2_1(kext,options,sysaxes,&
           iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,&
           R0,&
           Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind,&
           hmin,hmin_lon)
      end

      
      SUBROUTINE drift_bounce_orbit2_1(kext,options,sysaxes,&
           iyearsat,idoy,UT,xIN1,xIN2,xIN3,alpha,maginput,&
           R0,&
           Lm,Lstar,BLOCAL,BMIN,BMIR,XJ,posit,ind,&
           hmin,hmin_lon)
      
!     computes posit(3,1000,25), BLOCAL(1000,25) and ind(25) in same
!     format as drift_shell1 (note 25 not 48 azimuths)
!     also provides Bmin, Bmirror and usual stuff
!     allows user to specify R0 = surface below which a particle is considered
!     lost (R0 is in RE and the usual value is 1 but for the drift loss cone, one
!     may wish to allow R0<1)
!     also provides hmin and hmin_lon, the altitude and longitude (GDZ) of the
!     minimum altitude along the drift orbit
      
      
      IMPLICIT NONE
      INCLUDE 'variables.inc'
!     
!     declare inputs
      INTEGER*4    kext,k_ext,k_l,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      integer*4    idoy
      real*8     UT
      real*8     xIN1,xIN2,xIN3
      real*8     alpha
      real*8     maginput(25),R0
      REAL*8 hmin,hmin_lon
!     
!     
!     Declare internal variables
      INTEGER*4    kint,i,j,k
      INTEGER*4    ifail,Ilflag,t_resol,r_resol
      REAL*8     xGEO(3),xGEOmir(3), BL
      REAL*8     alti,lati,longi
!     
!     Declare output variables
      INTEGER*4  ind(25)
      REAL*8     posit(3,1000,25)
      REAL*8     BLOCAL(1000,25)
      REAL*8     BMIN,BMIR
      REAL*8     XJ
      REAL*8     Lm,Lstar
!     
      COMMON /magmod/k_ext,k_l,kint
      COMMON /flag_L/Ilflag
      integer*4 int_field_select, ext_field_select
!     
!     initialize outputs
      Lm=baddata
      Lstar=baddata
      XJ=baddata
      BMIN=baddata
      BMIR=baddata
      hmin = baddata
      hmin_lon = baddata
      do i=1,25
         ind(i)=0
         do j=1,1000
            BLOCAL(j,i)=baddata
            do k=1,3
               posit(k,j,i)=baddata
            enddo
         enddo
      enddo

      Ilflag=0
      k_ext=kext
      if (options(1).eq.0) options(1)=1
      if (options(3).lt.0 .or. options(3).gt.9) options(3)=0
      t_resol=options(3)+1
      r_resol=options(4)+1
      k_l=options(1)
!
	kint = int_field_select ( options(5) )
	k_ext = ext_field_select ( kext )
!
        CALL INITIZE
	
	call init_fields ( kint, iyearsat, idoy, ut, options(2) )
	
	call get_coordinates ( sysaxes, xIN1, xIN2, xIN3, & 
      alti, lati, longi, xGEO )
	    
	call set_magfield_inputs ( kext, maginput, ifail )
	
	if ( ifail.lt.0 ) RETURN
      if (kext .eq. 13 .or. kext .eq. 14) then !special script to read files and
         call INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
         call INIT_TS07D_TLPR
	     if ( ifail.lt.0 ) RETURN
      end if
!
!    

      CALL find_bm_nalpha(xGEO,1,alpha,BL,BMIR,xGEOmir)
      IF (Bmir.NE.baddata) THEN
         Ilflag=0
         call trace_drift_bounce_orbit_opt(t_resol,r_resol,&
                 xGEOmir,R0,Lm,Lstar,XJ,&
                 BLOCAL,Bmin,Bmir,posit,ind,hmin,hmin_lon)
         if (Lstar.eq.baddata) then
            hmin = baddata
            hmin_lon = baddata
         endif
         
      ENDIF
 99   continue
      
      end                       ! end subroutine drift_bounce_orbit1


!     --------------------------------------------------------------------
!     

      SUBROUTINE trace_drift_bounce_orbit(t_resol,r_resol,&
           lati,longi,alti,R0,Lm,Lstar,leI0,Bposit,&
           Bmin,Bmir,posit,Nposit,hmin,hmin_lon)
!     
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      
!     Parameters
      INTEGER*4  Nreb_def,Nder_def,Ntet_def
      PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
      
!     Input Variables
      INTEGER*4 t_resol,r_resol
      REAL*8     lati,longi,alti,R0
      REAL*8     xx0(3)
!     Output Variables       
      REAL*8     Lm,Lstar,leI0,Bmin,Bmir
      REAL*8     posit(3,20*Nreb_def,Nder_def)
      REAL*8     Bposit(20*Nreb_def,Nder_def)
      REAL*8     hmin,hmin_lon
      INTEGER*4  Nposit(Nder_def)

      CALL GDZ_GEO(lati,longi,alti,xx0(1),xx0(2),xx0(3))
!     
      call trace_drift_bounce_orbit_opt(t_resol,r_resol,&
           xx0,R0,Lm,Lstar,leI0,Bposit,&
           Bmin,Bmir,posit,Nposit,hmin,hmin_lon)

      RETURN
      END      

      SUBROUTINE trace_drift_bounce_orbit_opt(t_resol,r_resol,&
           xx0,R0,Lm,Lstar,leI0,Bposit,&
           Bmin,Bmir,posit,Nposit,hmin,hmin_lon)
!     
      USE fieldline_utils
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      
!     Parameters
      INTEGER*4  Nreb_def,Nder_def,Ntet_def
      PARAMETER (Nreb_def = 50, Nder_def = 25, Ntet_def = 720)
      
!     Input Variables
      INTEGER*4 t_resol,r_resol
      REAL*8     R0
      
!     Internal Variables
      INTEGER*4  Nder,Nreb,Ntet
      INTEGER*4  k_ext,k_l,kint,n_resol
      INTEGER*4  Nrebmax
      REAL*8     rr,rr2
      REAL*8     xx0(3),xx(3),x1(3),x2(3)
      REAL*8     xmin(3)
      REAL*8     B(3),Bl,B0,B1,B3, bfoot(3), bfootmag
      REAL*8     dsreb0,dsreb,smin
      
      INTEGER*4  I,J,K,Iflag,Iflag_I,Ilflag,Ifail
      INTEGER*4  Ibounce_flag,hemi_flag,stop_type 
      !istore is azimuth cursor
      INTEGER*4  istore
      REAL*8     Lb, stop_value
      REAL*8     leI,leI1
      REAL*8     XY,YY
      REAL*8     aa,bb
      
!     
      REAL*8     tt
      REAL*8     tet(10*Nder_def),phi(10*Nder_def)
      REAL*8     tetl,tet1,dtet,lasttet
      REAL*8     somme,BrR2
        REAL*8     pi,rad
        common /rconst/rad,pi

! variables to deal with leI~0 case
      REAL*8     x1old(3)
      INTEGER*4  I0flag
!     
      REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
      
!     Output Variables       
      REAL*8     Lm,Lstar,leI0,Bmin,Bmir
      REAL*8     posit(3,20*Nreb_def,Nder_def)
      REAL*8     Bposit(20*Nreb_def,Nder_def)
      REAL*8     hmin,hmin_lon
      INTEGER*4  Nposit(Nder_def), tet_count
!     
      COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
      COMMON /calotte/tet
      COMMON /flag_L/Ilflag
      COMMON /magmod/k_ext,k_l,kint
!     
!     


      Nder=Nder_def*r_resol     ! longitude steps
      Nreb=Nreb_def             ! steps along field line
      Ntet=Ntet_def*t_resol     ! latitude steps
      dtet = pi/Ntet            ! theta (latitude) step
!     
      Nrebmax = 20*Nreb         ! maximum steps along field line
!     
!
!     initialize hmin,hmin_lon to starting point
      hmin = baddata
      call check_hmin(xx0(1),xx0(2),xx0(3),hmin,hmin_lon)

      CALL COMPUTE_L_DIPOLE(xx0, Lb)
!     
      CALL CHAMP(xx0,B,B0,Ifail)
      Bmir = B0 ! local field at starting point - for Bmir=0 case when alpha=90
      IF (Ifail.LT.0) THEN
         Ilflag = 0
         RETURN
      ENDIF
      Bmin = B0

      

!     
      dsreb0 = Lb/Nreb           ! step size dipole L / Nsteps
      dsreb = dsreb0
      CALL COMPUTE_FIELDLINE_FROM_MIR(xx0, dsreb, nrebmax, xmin, Bmin, leI0, Lm)

      !     calcul du point sur la ligne de champ a la surface de la terre du
      !     cote nord
      !     (compute the point nothern footpoint at Earth's surface)

      stop_type = 0 ! GEO radius stop cdt
      stop_value = R0 ! Rstop = R0
      hemi_flag = 1 ! northern hemi
      CALL FIND_FIELDLINE_FOOT(xx0, dsreb, nrebmax, stop_type, stop_value, hemi_flag, &
                                    x2, bfoot, bfootmag)

      if (k_l .eq.0) then
         Ilflag = 0
         RETURN
      endif
      IF (ABS(Lm) .GT. 10.D0) THEN
         Ilflag = 0
         RETURN
      ENDIF

      rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
      tet(1) = ACOS(x2(3)/rr)
      phi(1) = ATAN2(x2(2),x2(1))

!     et on tourne -> on se decale sur la surface en phi et on cherche teta
!     pour avoir leI0 et B0 constants
!     (find the thetta/phi contour on the surface of the earth that conserves
!     B0= Bmirror and leI0=I)
      call trace_bounce_orbit(x2,B0,1,dsreb0,R0,&
           Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
      if (Ibounce_flag.ne.1) then ! try again from starting point
         call trace_bounce_orbit(xx0,B0,1,dsreb0,R0,&
                 Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
      endif
      if (Ibounce_flag.ne.1) then
         Ilflag = 0
         RETURN
      endif

      istore = 2 ! istore=1 already stored in first bounce orbit trace
      dsreb = -dsreb
      DO I = 2,Nder
         phi(I) = phi(I-1)+2.D0*pi/Nder
         Iflag_I = 0
         I0flag = 0
         IF (Ilflag.EQ.0) THEN
            tetl = tet(I-1)
            IF (I.GT.2) tetl = 2.D0*tet(I-1)-tet(I-2)
            tet1 = tetl
         ELSE
            tetl = tet(I)
            tet1 = tetl
         ENDIF
         leI1 = baddata
!     
         tet_count = 0          ! number of tries of different theta
 107     CONTINUE
         tet_count = tet_count+1
         leI = baddata
         IF (tetl.GT.pi .OR. tetl.LT.0.D0) GOTO 108
         x1(1) = R0*SIN(tetl)*COS(phi(I))
         x1(2) = R0*SIN(tetl)*SIN(phi(I))
         x1(3) = R0*COS(tetl)
         Iflag = 0
         CALL CHAMP(x1,B,B1,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         if (B1.lt.Bmir) then
            if (tet_count.gt.Ntet*2) goto 108 ! abort
            ! starting point is in middle of field line. Try again at lower colat
            tetl = tetl-dtet
            goto 107
         endif
!     
         dsreb = dsreb/abs(dsreb)*dsreb0
         DO J = 1,Nrebmax
 109        CALL sksyst_var(dsreb,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) THEN
               Ilflag = 0
               RETURN
            ENDIF
            rr2 = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
            IF (Bl.LT.B0) THEN
               lasttet = tetl
               IF (Iflag .EQ. 0) THEN
                  CALL CHAMP(x1,B,B1,Ifail)
                  IF (Ifail.LT.0) THEN
                     Ilflag = 0
                     RETURN
                  ENDIF
                  leI = 0.5D0*SQRT(1.D0-Bl/B0)*(1.D0+(Bl-B0)/(Bl-B1))&
                                   *abs(dsreb)
                  Iflag = 1
               ELSE
                  leI = leI+SQRT(1.D0-Bl/B0)*abs(dsreb)
               ENDIF
            ENDIF
            IF (Bl.GT.B0 .AND. Iflag.EQ.1) GOTO 103
            IF (rr2.LT.R0) GOTO 103
            x1(1) = x2(1)
            x1(2) = x2(2)
            x1(3) = x2(3)
            B1 = Bl
         ENDDO
 103     CONTINUE
         if ((Iflag.eq.0).and.(J.LT.Nrebmax)) then ! never got to B < Bmirror
            leI=0 ! assume particle is equatorially mirroring
                  ! Not technically true, but effectively gives desired result
                  ! (i.e., leI=0 for field lines way inside trajectory)
            I0flag = 1
         else
            I0flag = 0
!     Pourquoi? (why? I don't know!)
            IF (rr2.LT.R0) THEN
               leI = baddata
            ENDIF
            IF (J.LT.Nrebmax .AND. rr2.GE.R0) THEN
               CALL CHAMP(x1,B,B1,Ifail)
               IF (Ifail.LT.0) THEN
                  Ilflag = 0
                  RETURN
               ENDIF
               lasttet = tetl
               leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)*&
                             abs(dsreb)
            ENDIF
         endif
!     
         IF (Iflag_I .EQ.0) THEN
            IF (J.GE.Nrebmax) THEN
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
         GOTO 107
 108     CONTINUE
         IF ((J.GE.Nrebmax .AND. leI.GT.0.D0).or.(leI.LT.0.D0)) THEN
            Ilflag = 0
            RETURN
         ENDIF

         tetl = 0.5D0*(tetl+tet1) ! both tetl and tet1 work, so average
         tet(I) = tetl
         x1(1) = R0*SIN(tet(I))*COS(phi(I))
         x1(2) = R0*SIN(tet(I))*SIN(phi(I))
         x1(3) = R0*COS(tet(I))
         CALL CHAMP(x1,B,Bl,Ifail)
         IF (Ifail.LT.0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         IF (Bl.LT.B0) THEN
            Ilflag = 0
            RETURN
         ENDIF
         
! Trace bounce orbit on field line from x1
         call trace_bounce_orbit(x1,Bmir,istore,dsreb0,R0,&
                 Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
         if (Ibounce_flag.ne.1) then
            ! try again from lasttet, which worked
            tet(I) = lasttet
            x1(1) = R0*SIN(tet(I))*COS(phi(I))
            x1(2) = R0*SIN(tet(I))*SIN(phi(I))
            x1(3) = R0*COS(tet(I))
            call trace_bounce_orbit(x1,Bmir,istore,dsreb0,R0,&
                       Bposit,posit,Nposit,Ibounce_flag,hmin,hmin_lon)
         endif
         if (Ibounce_flag.ne.1) then
            Ilflag = 0
            RETURN
         endif
         if (MOD(I-1,r_resol).eq.0) then
            istore = istore+1 ! next time, store in next array step
         endif                  ! mode(I,r_resol)
      ENDDO                     ! end of do I = 2,Nder

      CALL INTEGRATE_FLUX_ON_POLAR_CAP(Nder, Ntet, tet, phi, somme, R0)

      if (k_l .eq.1) Lstar = 2.D0*pi*Bo/somme
      if (k_l .eq.2) Lstar = somme ! Phi and not Lstar
      IF (Lm.LT.0.D0) Lstar = -Lstar
      Ilflag = 1
!     
      END
      
      SUBROUTINE trace_bounce_orbit(xstart,Bmirror,istore,dsreb0,&
           R0,Bposit,posit,Nposit,Iflag,hmin,hmin_lon)

      IMPLICIT NONE

!     Declare input variables
      !xstart in GEO, Bmirror is particle mirror field strength
      !istore stores bounce orbit
      REAL*8 xstart(3)
      REAL*8 Bmirror
      INTEGER*4 istore
      ! hmin,hmin_lon are also inputs

!     Declare output variables
      INTEGER*4  Nposit(25),Iflag
      REAL*8     posit(3,1000,25)
      REAL*8     Bposit(1000,25),R0,hmin,hmin_lon

!     Declare internal variables
      INTEGER*4 i,j,k,Ifail
      REAL*8 lati,longi,alti,Bl,Bmir
      REAL*8 xmir(3)
      REAL*8 alpha
      REAL*8 dsreb0,dsreb,x1(3),x2(3),B1,B2,Bvec(3)
      ! store or not?
      INTEGER*4 store

      j = 1 ! TPO 5 Jan 2020 seems like j was never initialized

      ! Iflag=1 means success, set to fail until succeeded
      Iflag = 0

      if ((istore.gt.25).or.(istore.lt.1)) then
         store = 0
      else
         store = 1
      endif

!     trace up to mirror point, leave starting pointin Bmir, xmir

      call champ(xstart,Bvec,B1,Ifail)

      dsreb = -dsreb0 ! assume starting at northern R=1 foot point

      do i = 1,3
         x1(i) = xstart(i)
      enddo

      call sksyst_var(dsreb,x1,x2,B2,Ifail)
      if (B2.GT.B1) then
         dsreb = -dsreb ! going wrong way
      endif

      call sksyst_var(dsreb,x1,x2,B2,Ifail)
      if (B2.GT.B1) then ! already at local min w/in resolution of dsreb
         ! store & return
         if (store.eq.1) then
            Bposit(1,istore) = B2
            Nposit(istore) = 1
            do k = 1,3
               posit(k,j,istore) = x1(k)
            enddo
         endif
         call check_hmin(x1(1),x1(2),x1(3),hmin,hmin_lon)
         Iflag=1
         return
      endif

      do j = 1,999
         call sksyst_var(dsreb,x1,x2,B2,Ifail)
         if (B2.GT.Bmirror) then ! still haven't crossed Bmirror
            do i = 1,3
               x1(i) = x2(i)
            enddo
            B1 = B2
         else
            do i = 1,3
               xmir(i) = x2(i)
            enddo
            Bmir = B2
            if (abs((Bmir/Bmirror-1.0D0)).lt.1.0D-6) then
               goto 101         ! good convergence
            endif
            dsreb = dsreb/2.D0 ! try again with smaller step (converges logarithmically)
         endif
      enddo
      return ! did not converge (Iflag=0)
 101  continue

      call check_hmin(xmir(1),xmir(2),xmir(3),hmin,hmin_lon)
      dsreb = -dsreb0

!     trace bounce trajectory
      call sksyst_var(dsreb,xmir,x2,B2,Ifail)
      if (B2.GT.Bmir) then
         dsreb = -dsreb ! going wrong way
      endif

!     trace to opposite mirror point
!     set up trace
      B2 = BMIR
      do k=1,3
         x1(k)=xmir(k)
      enddo
!     do trace
      do j=1,999
         if (store.eq.1) then
            Bposit(j,istore) = B2
            Nposit(istore) = j
            do k = 1,3
               posit(k,j,istore) = x1(k)
            enddo
         endif
         call sksyst_var(dsreb,x1,x2,B2,Ifail)
         if (Ifail.LT.0) then
            Iflag = 0
            return
         endif
         if (B2.ge.BMIR) goto 201
         do k=1,3
            x1(k)=x2(k)
         enddo
      enddo
 201  continue
      if (B2.lt.BMIR) then
         Iflag = 0 ! failed to get past BMIR
         return
      elseif (B2.gt.BMIR) then
!     finish trace to Bm: Bm between x1 and x2
         BL = B2 ! closest stored point so far (x1)
         j = j+1                ! prepare to store in next element
         do i = 1,10            ! this converges logarithmically, so reduces step size by up to 2^10
            dsreb = dsreb/2.0D0 ! restart with half step size
            call sksyst_var(dsreb,x1,x2,B2,Ifail)
            if (B2.LT.Bmirror) then ! still inside bounce orbit
               do k=1,3
                  x1(k) = x2(k)
               enddo
               if (store.eq.1) then ! got closer, so store x1
                  Bposit(j,istore) = B2
                  Nposit(istore) = j
                  do k=1,3
                     posit(k,j,istore) = x1(k)
                  enddo
               endif
            endif
         enddo
         call check_hmin(x1(1),x1(2),x1(3),hmin,hmin_lon)
      endif

      Iflag = 1
      END ! end sub trace_bounce_orbit

      subroutine check_hmin(x,y,z,hmin,hmin_lon)
      ! check alt/lat/lon at x,y,z and replace hmin,hmin_lon if needed
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      ! inputs - x,y,z GEO
      real*8 x,y,z
      
      ! outputs
      real*8 hmin,hmin_lon

      ! locals
      real*8 lati,longi,alti

      call geo_gdz(x,y,z,lati,longi,alti)
      if ((hmin.eq.baddata).or.(alti .lt. hmin)) then
         hmin = alti
         hmin_lon = longi
      endif
      end


      subroutine sksyst_var(dsreb,x0,x1,B,Ifail)
      ! sksyt but takes 1/10 size steps for R<1
      IMPLICIT NONE

      REAL*8 dsreb,x0(3),x1(3),B,xtmp(3)
      INTEGER*4 Ifail,i,j

      if ((x0(1)**2+x0(2)**2+x0(3)**2).ge.1.0) then
         ! R>= 1.0, so use regular dsreb
         call sksyst(dsreb,x0,x1,B,Ifail)
         return
      endif
      ! else, take 10 small steps
      do i = 1,3
         xtmp(i) = x0(i)
      enddo
      do j = 1,10
         call sksyst(dsreb/10.D0,xtmp,x1,B,Ifail)
         if (Ifail.LT.0) return
         do i = 1,3
            xtmp(i) = x1(i)
         enddo
      enddo
      end
