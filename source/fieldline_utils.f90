!   Module containing functions relative to field lines
!   Created by : L. Nicolle, 2025

MODULE fieldline_utils
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: FIND_MAGEQUATOR
    PUBLIC :: FIND_MIRROR_POINT
    PUBLIC :: FIND_FIELDLINE_FOOT
    PUBLIC :: COMPUTE_L_DIPOLE
    PUBLIC :: COMPUTE_L_MCILWAIN_HILTON
    PUBLIC :: COMPUTE_FIELDLINE_FROM_MIR
    PUBLIC :: INTEGRATE_FLUX_ON_POLAR_CAP
    PUBLIC :: COMPUTE_DRIFTSHELL_FOOTPOINTS
    

    CONTAINS

    !=======================================================================
    ! Purpose : compute L shell assuming a perfect dipole
    ! Input : 
    !   posGEO(3) - position in GEO coordinates
    ! Output :
    !   Ld - L shell (perfect dipole)
    !=======================================================================
    SUBROUTINE COMPUTE_L_DIPOLE(posGEO, Ld)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: posGEO(3)
        REAL(8), INTENT(OUT) :: Ld

        REAL(8) :: r, theta, xx(3)

        CALL GEO_SM(posGEO,xx)

        ! Distance from Earth center
        r = SQRT(xx(1)**2 + xx(2)**2 + xx(3)**2)

        ! Colatitude (0 = north pole, pi/2 = equator)
        theta = ACOS(xx(3) / r)

        ! Dipole L parameter
        Ld = r / (SIN(theta)**2)
    
    END SUBROUTINE

    !=======================================================================
    ! Purpose : compute L shell using the Mc Ilwain-Hilton approximation
    ! See : https://doi.org/10.1029/JA076i028p06952
    ! Input : 
    !   leI - Normalized second adiabatic invariant
    !   B0 - Magnitude of the magnetic field at mirror point
    ! Output :
    !   Lmh - L shell (Mc Ilwain Hilton)
    !=======================================================================
    SUBROUTINE COMPUTE_L_MCILWAIN_HILTON(leI, Bmirr, Lmh)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: leI, Bmirr
        REAL(8), INTENT(OUT) :: Lmh

        REAL(8) :: XY, YY
        REAL(8)     Bo,xc,yc,zc,ct,st,cp,sp
        COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp

        XY = leI*leI*leI*Bmirr/Bo
        YY = 1.D0 + 1.35047D0*XY**(1.D0/3.D0) &
                + 0.465376D0*XY**(2.D0/3.D0) &
                + 0.0475455D0*XY
        Lmh = (Bo*YY/Bmirr)**(1.D0/3.D0)
    END SUBROUTINE

    !=======================================================================
    ! Purpose : Returns field line properties from a mirror point
    ! Input : 
    !   xmir(3) - position (GEO) of mirror point
    !   dsreb - integration arc-length step along the magnetic field line
    !   rebmax - maximum number of dsreb steps along the magnetic field line
    ! Output :
    !   xeq(3) - position (GEO) of the magnetic equator
    !   beq - magnitude of the magnetic field at xeq
    !   lei - normalized second adiabatic invariant
    !   Lmh - Lshell computed using the McIlwain-Hilton approximation
    !=======================================================================
    SUBROUTINE COMPUTE_FIELDLINE_FROM_MIR(xmir, dsreb, nrebmax, xeq, beq, lei, Lmh)

        IMPLICIT NONE
        INCLUDE 'variables.inc'

        REAL(8), INTENT(IN) :: xmir(3)
        REAL(8), INTENT(IN) :: dsreb
        INTEGER(4), INTENT(IN) :: nrebmax
        REAL(8), INTENT(OUT) :: xeq(3)
        REAL(8), INTENT(OUT) :: beq, lei, Lmh

        REAL(8) :: xx0(3), x1(3), x2(3), B(3)
        REAL(8) :: B0, Bl, B1, B3, XY, YY, aa, bb, smin, leI0, dsreb_loc
        INTEGER(4) :: Ifail, Ilflag, I, J

        REAL*8     Bo,xc,yc,zc,ct,st,cp,sp
        INTEGER*4  k_ext,k_l,kint

        COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
        COMMON /flag_L/Ilflag
        COMMON /magmod/k_ext,k_l,kint

        DO I = 1, 3
            xx0(I) = xmir(I)
        ENDDO

        dsreb_loc = dsreb

        CALL CHAMP(xx0,B,B0,Ifail)

! calcul du sens du depart
!"calculation of the direction of the departure"
!
        CALL sksyst(-dsreb_loc,xx0,x1,Bl,Ifail)
            IF (Ifail.LT.0) THEN
                leI0 = baddata
                beq = baddata
                Ilflag = 0
	            RETURN
            ENDIF
        B1 = Bl
        CALL sksyst(dsreb_loc,xx0,x2,Bl,Ifail)
        IF (Ifail.LT.0) THEN
            leI0 = baddata
            beq = baddata
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
            beq = B0 - aa*smin*smin
            leI0 = SQRT(1.D0-beq/B0)*2.D0*ABS(smin*dsreb_loc)
            Lmh = (Bo/beq)**(1.D0/3.D0)
!         write(6,*)'L McIlwain eq ',B0,leI0,Lm
        GOTO 100
        ENDIF
        IF (B3.GT.B1) THEN
            dsreb_loc = -dsreb_loc
        ENDIF

!
! calcul de la ligne de champ et de I
! "calculation of the field line and I"
!
        beq = B0
        leI = 0.D0
        DO I = 1,3
            x1(I)  = xx0(I)
        ENDDO
!
!       write(6,*)dsreb_loc
        DO J = 1,Nrebmax
            CALL sksyst(dsreb_loc,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) THEN
                leI0 = baddata
                beq = baddata
                Ilflag = 0
                RETURN
            ENDIF
!	 write(6,*)J,x1(1),x1(2),x1(3),Bl
            IF (Bl.LT.beq) THEN
                xeq(1) = x2(1)
                xeq(2) = x2(2)
                xeq(3) = x2(3)
                beq = Bl
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
        IF (J.GE.Nrebmax) THEN !open field line
            leI0 = baddata
            beq = baddata
            Ilflag = 0
            RETURN
        ENDIF
!
        leI = leI+0.5D0*SQRT(1.D0-B1/B0)*(B0-Bl)/(Bl-B1)
        leI = leI*ABS(dsreb_loc)
        leI0 = leI

        CALL COMPUTE_L_MCILWAIN_HILTON(leI, B0, Lmh)

!
! calcul de Bmin
!

        CALL sksyst(dsreb_loc,xeq,x1,B3,Ifail)
        IF (Ifail.LT.0) THEN
            beq = baddata
            Ilflag = 0
            RETURN
        ENDIF
        CALL sksyst(-dsreb_loc,xeq,x1,B1,Ifail)
        IF (Ifail.LT.0) THEN
            beq = baddata
            Ilflag = 0
            RETURN
        ENDIF
        aa = 0.5D0*(B3+B1-2.D0*beq)
        bb = 0.5D0*(B3-B1)
        smin = -0.5D0*bb/aa
        beq = beq - aa*smin*smin
        IF (x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3).LT.1.D0) THEN
            Lmh = -Lmh
        ENDIF

        ! write(6,*) "Beq", beq

!
100     CONTINUE
        if (k_l .eq.0) then
            Ilflag = 0
            RETURN
        endif
        IF (ABS(Lmh) .GT. 10.D0) THEN
            Ilflag = 0
            RETURN
        ENDIF

    END SUBROUTINE




    SUBROUTINE FIND_MAGEQUATOR(xx0, blocal, dsreb, nrebmax, xmin, bmin)

        IMPLICIT NONE
        INCLUDE 'variables.inc'
        REAL(8), INTENT(IN) :: xx0(3), blocal, dsreb
        INTEGER(4), INTENT(IN) :: nrebmax
        REAL(8), INTENT(OUT) :: xmin(3), Bmin

        REAL(8) :: dsreb_loc, x1(3), x2(3), Bl, B1, B3, B0
        REAL(8) :: aa, bb, smin
        REAL(8) :: posit(3)
        INTEGER(4) :: Ifail, I, J


        B0 = blocal
        dsreb_loc = dsreb

        ! calcul du sens du depart 

        CALL sksyst(-dsreb_loc,xx0,x1,Bl,Ifail)
        IF (Ifail.LT.0) THEN
            Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	        RETURN
        ENDIF
        B1 = Bl
        CALL sksyst(dsreb_loc,xx0,x2,Bl,Ifail)
        IF (Ifail.LT.0) THEN
	        Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	        RETURN
        ENDIF
        B3 = Bl

        ! attention cas equatorial
        IF(B1.GT.B0 .AND. B3.GT.B0)THEN
            aa = 0.5D0*(B3+B1-2.D0*B0)
            bb = 0.5D0*(B3-B1)
            smin = -0.5D0*bb/aa
            Bmin = B0 - aa*smin*smin
            CALL sksyst(smin*dsreb_loc,xx0,posit,Bl,Ifail)
            IF (Ifail.LT.0) THEN
                Bmin=baddata
                posit(1) = baddata
                posit(2) = baddata
                posit(3) = baddata
            RETURN
            ENDIF
	        RETURN
        ENDIF
        IF (B3.GT.B1) THEN
            dsreb_loc = -dsreb_loc
        ENDIF

        ! calcul de la ligne de champ et de I
        Bmin = B0
        DO I = 1,3
            x1(I)  = xx0(I)
        ENDDO

        DO J = 1,Nrebmax
            CALL sksyst(dsreb_loc,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) THEN
	            Bmin=baddata
                posit(1) = baddata
                posit(2) = baddata
                posit(3) = baddata
	            RETURN
            ENDIF
            !  bmin, xmin update
            IF (Bl.LT.Bmin) THEN
                xmin(1) = x2(1)
                xmin(2) = x2(2)
                xmin(3) = x2(3)
                Bmin = Bl
            ENDIF

            !found lowest
            IF (Bl.GT.B0) GOTO 20
	        x1(1) = x2(1)
	        x1(2) = x2(2)
	        x1(3) = x2(3)
	        B1 = Bl
        ENDDO
20      CONTINUE

        IF (J.GE.Nrebmax) THEN !open field line
            Bmin = baddata
            posit(1)=baddata
            posit(2)=baddata
            posit(3)=baddata
            RETURN
        ENDIF

! calcul de Bmin
        CALL sksyst(dsreb_loc,xmin,x1,B3,Ifail)
        IF (Ifail.LT.0) THEN
	        Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	        RETURN
        ENDIF
        CALL sksyst(-dsreb_loc,xmin,x1,B1,Ifail)
        IF (Ifail.LT.0) THEN
	        Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	        RETURN
        ENDIF
        aa = 0.5D0*(B3+B1-2.D0*Bmin)
        bb = 0.5D0*(B3-B1)
        smin = -0.5D0*bb/aa
        Bmin = Bmin - aa*smin*smin
        CALL sksyst(smin*dsreb_loc,xmin,posit,Bl,Ifail)
        IF (Ifail.LT.0) THEN
	        Bmin=baddata
            posit(1) = baddata
            posit(2) = baddata
            posit(3) = baddata
	        RETURN
        ENDIF
    END SUBROUTINE

    SUBROUTINE FIND_MIRROR_POINT(pos)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: pos(3)
    END SUBROUTINE

    !=======================================================================
    ! Purpose    : Compute the magnetic footpoint (Earth intersection
    !              points) of the field line corresponding to a given
    !              GEO coordinate
    !
    ! Inputs     :
    !   xGEO(3) - position (GEO) on the field line
    !   dsreb - integration arc-length step along the magnetic field line
    !   nrebmax - maximum number of dsreb steps along the magnetic field line
    !   stop_type :
    !       - 0 : stop at a given GEO radius
    !       - 1 : stop at a given GDZ alti
    !   stop_value :
    !       - stop_value = stop_radius (Re) if stop_type = 0
    !       - stop_value = stop_alt (km) if stop_type = 1
    !   hemi_flag : hemishere flag, specifies hemisphere of foot point
    !      0 - same Hemisphere as start point
    !      +1 - Northern Hemisphere
    !      -1 - Southern Hemisphere
    !      2 - opposite Hemisphere as start point
    !
    ! Outputs    :
    !    xfoot - footpoint (GEO (Re) if stop_type = 0, GDZ (alt, lat, long) if stop_type = 1)
    !    bfoot - magnetic field vector at foot point(nT, GEO)
    !    bfootmag - Magnetic field at foot point (nT)
    !=======================================================================
    SUBROUTINE FIND_FIELDLINE_FOOT(xGEO, dsreb, nrebmax, stop_type, stop_value, hemi_flag, &
                                    xfoot, bfoot, bfootmag)
        IMPLICIT NONE
        INCLUDE 'variables.inc'
        REAL(8), INTENT(IN) :: xGEO(3), dsreb, stop_value
        INTEGER(4), INTENT(IN) :: stop_type, nrebmax, hemi_flag
        REAL(8), INTENT(OUT) :: xfoot(3), bfoot(3), bfootmag

        REAL(8) :: x1(3), x2(3), x3(3), dsreb_loc, rr, Bl, B1, B3, smin
        INTEGER(4) :: I, J, Ifail
        LOGICAL :: UNDER_LIMIT

        ! init
        DO I = 1,3
            x1(I) = xGEO(I)
            xfoot(I) = baddata
            bfoot(I) = baddata
        ENDDO
        bfootmag = baddata
        dsreb_loc = abs(dsreb)

        ! check if initial point is already under limit
        IF (stop_type .eq. 0) THEN
            rr = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
            IF (rr .le. stop_value) RETURN
        ELSE IF (stop_type .eq. 1) THEN
            CALL GEO_GDZ(x1(1), x1(2), x1(3), &
                        x2(2), x2(3), x2(1))
            IF (x1(1) .le. stop_value) RETURN
        ENDIF

        dsreb_loc = min(dsreb_loc, 1.0D0)

        ! compute direction of departure
        ! get field one step southward (B1)
        CALL SKSYST(-dsreb_loc, x1, x2, B1, Ifail)
        IF (Ifail .lt. 0) RETURN

        ! get field one step northward (B3)
        CALL SKSYST(dsreb_loc, x1, x2, B3, Ifail)
        IF (Ifail .lt. 0) RETURN
        
        ! compute direction considering hemi_flag
        ! if hemi = 0 and already in northern hemi -> go north -> positive dsreb
        ! if hemi = 0 and in southern hemi -> go south -> negative dsreb
        IF ((hemi_flag .eq. 0) .and. (B1 .gt. B3)) dsreb_loc = -dsreb_loc

        ! if hemi = 1, northern hemi -> go north -> positive dsreb
        ! if hemi = -1, southern hemi -> go south -> negative dsreb
        IF (hemi_flag .eq. -1) dsreb_loc = -dsreb_loc

        ! if hemi = 2, and in southern hemi -> go north -> positive dsreb
        ! if hemi = 2, and in northern hemi -> go south -> negative dsreb
        IF ((hemi_flag .eq. 2) .and. (B1 .lt. B3)) dsreb_loc = -dsreb_loc

        ! follow field line until we reach stopping condition
15      UNDER_LIMIT = .false.
        DO J = 1, nrebmax
            CALL SKSYST(dsreb_loc, x1, x2, Bl, Ifail)
            IF (Ifail .lt. 0) RETURN

            ! Check stop condition
            IF (stop_type .eq. 0) THEN
                rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
                IF (rr.LT.stop_value) UNDER_LIMIT = .true.
            ELSE IF (stop_type .eq. 1) THEN
                call geo_gdz(x2(1),x2(2),x2(3),x3(2),x3(3),x3(1))
                if (x3(1).LE.stop_value) UNDER_LIMIT = .true.
            ENDIF
            ! exit loop
            IF (UNDER_LIMIT) GOTO 20
            DO I = 1, 3
                x1(I) = x2(I)
            ENDDO
        ENDDO
20      CONTINUE

        IF (J .ge. nrebmax) RETURN ! open field line

        ! correct foot position
        ! legacy correction for GEO stop condition
        IF (stop_type .eq. 0) THEN
            ! compute fraction of step to reach the earth (linear approx)
            smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
            smin = (1.D0-smin)/(rr-smin)
            CALL sksyst(smin*dsreb_loc,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) RETURN
            call CHAMP(x2,bfoot,bfootmag,Ifail)
            IF (Ifail.LT.0) RETURN

            DO I = 1, 3
                xfoot(I) = x2(I)
            ENDDO
            RETURN
            ! rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
            ! tet = ACOS(x2(3)/rr)
            ! phi = ATAN2(x2(2),x2(1))
        
        ELSE IF (stop_type .eq. 1) THEN
            ! footpoint is between x1 and x2
            IF(abs(x3(1)-stop_value).le.1.0) then
                !get B field at x2
                call CHAMP(x2,bfoot,bfootmag,Ifail)
                IF(Ifail.LT.0) RETURN
                DO I = 1, 3
                    xfoot(I) = x2(I)
                ENDDO
                RETURN
            ELSE  ! try loop again with smaller step
                dsreb_loc = dsreb_loc/100.0
                goto 15
            ENDIF
        ENDIF

    END SUBROUTINE

!     !=======================================================================
!     ! Purpose : Returns the position (colatitude, longitude GEO) of the footpoint (north)
!     ! Input : 
!     !   pos(3) - position (GEO) on the field line
!     !   dsreb - integration arc-length step along the magnetic field line
!     !   nrebmax - maximum number of dsreb steps along the magnetic field line
!     ! Output :
!     !   tet - colatitude (GEO) of the footpoint (North)
!     !   phi - longitude (GEO) of the footpoint (North)
!     !=======================================================================
!     SUBROUTINE FIND_FIELDLINE_FOOT_FROM_POS(pos, dsreb, nrebmax, tet, phi)
        
!         IMPLICIT NONE
!         REAL(8), INTENT(IN) :: pos(3)
!         REAL(8), INTENT(IN) :: dsreb
!         INTEGER(4), INTENT(IN) :: nrebmax
!         REAL(8), INTENT(OUT) :: tet, phi

!         REAL(8) :: x1(3), x2(3), rr, Bl, smin, dsreb_loc

!         INTEGER(4) :: Ilflag, Ifail
!         INTEGER(4) :: I, J
!         COMMON /flag_L/Ilflag
!         ! "calculation of the point of the field line on the surface 
!         !    of the earth of the northern slope(?)"
!         ! calcul du point sur la ligne de champ a la surface de la terre du
!         ! cote nord
!         !
!         DO I = 1,3
!                 x1(I)  = pos(I)
!         ENDDO
!         dsreb_loc = ABS(dsreb)
!         DO J = 1,Nrebmax
!             CALL sksyst(dsreb_loc,x1,x2,Bl,Ifail)
!             IF (Ifail.LT.0) THEN
!                 Ilflag = 0
!                 RETURN
!             ENDIF
!             rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
!             IF (rr.LT.1.D0) GOTO 102

!         ! need to return if we past here at Nrebmax, as we will end
!         ! up with /0 error at smin below

!             IF (J.EQ.Nrebmax) THEN
!                 Ilflag = 0
!                 RETURN
!             ENDIF
!             x1(1) = x2(1)
!             x1(2) = x2(2)
!             x1(3) = x2(3)
!         ENDDO
! 102     CONTINUE

!         ! compute fraction of step to reach the earth (linear approx)
!         smin = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
!         smin = (1.D0-smin)/(rr-smin)
!         CALL sksyst(smin*dsreb_loc,x1,x2,Bl,Ifail)
!         IF (Ifail.LT.0) THEN
!                 Ilflag = 0
!                 RETURN
!         ENDIF
!         rr = sqrt(x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3))
!         tet = ACOS(x2(3)/rr)
!         phi = ATAN2(x2(2),x2(1))

!     END SUBROUTINE


    SUBROUTINE INTEGRATE_J_ON_FIELDLINE_FROM_FOOTPOINT(xfoot, bmirr, dsreb, &
            nrebmax, rr2, leI, Ifail, exceeded_nmaxreb)
        
        IMPLICIT NONE
        INCLUDE 'variables.inc'

        REAL(8), INTENT(IN) :: xfoot(3), bmirr, dsreb
        INTEGER(4), INTENT(IN) :: nrebmax
        REAL(8), INTENT(OUT) :: rr2, leI
        INTEGER(4), INTENT(OUT) :: Ifail
        LOGICAL, INTENT(OUT) :: EXCEEDED_NMAXREB

        INTEGER(4) :: I, J, Ilflag, Iflag
        REAL(8) :: x1(3), x2(3), B(3), Bl, B1

        COMMON /flag_L/Ilflag

        Iflag = 0
        leI = baddata
        EXCEEDED_NMAXREB = .false.

        DO I = 1, 3
            x1(I) = xfoot(I)
        ENDDO

        DO J = 1,Nrebmax
            ! move on the field line
            CALL sksyst(dsreb,x1,x2,Bl,Ifail)
            IF (Ifail.LT.0) THEN
                Ilflag = 0
                RETURN
            ENDIF
            rr2 = x2(1)*x2(1)+x2(2)*x2(2)+x2(3)*x2(3)
            ! means we are between mirror points (so integrate second invar)
            IF (Bl.LT.bmirr) THEN
                IF (Iflag .EQ. 0) THEN
                    CALL CHAMP(x1,B,B1,Ifail)
                    IF (Ifail.LT.0) THEN
                        Ilflag = 0
                        RETURN
                    ENDIF
                    leI = 0.5D0*SQRT(1.D0-Bl/bmirr)*(1.D0+(Bl-bmirr)/(Bl-B1))
                    Iflag = 1
                ELSE
                    leI = leI+SQRT(1.D0-Bl/bmirr)
                ENDIF
            ENDIF

            ! end integration
            IF (Bl.GT.bmirr .AND. Iflag.EQ.1) GOTO 103

            ! inside earth fail condition
            IF (rr2.LT.1.D0) GOTO 103
            x1(1) = x2(1)
            x1(2) = x2(2)
            x1(3) = x2(3)
        ENDDO
        IF (J.GE.Nrebmax) EXCEEDED_NMAXREB = .true.


103     CONTINUE
! Pourquoi?  "why?"
        IF (rr2.LT.1.D0) THEN
            leI = baddata
        ENDIF

        ! integration worked, divide by the arc-length step
        IF (J.LT.Nrebmax .AND. rr2.GE.1.D0) THEN
            CALL CHAMP(x1,B,B1,Ifail)
            IF (Ifail.LT.0) THEN
                Ilflag = 0
                RETURN
            ENDIF
            leI = leI+0.5D0*SQRT(1.D0-B1/bmirr)*(bmirr-Bl)/(Bl-B1)
            leI = leI*ABS(dsreb)
        ENDIF

    END SUBROUTINE


    SUBROUTINE COMPUTE_DRIFTSHELL_FOOTPOINTS(pos, leI0, bmirr, dsreb, nrebmax, tet, phi, nder, ntet)
        
        IMPLICIT NONE
        INCLUDE 'variables.inc'

        REAL(8), INTENT(IN) :: pos(3), lei0, bmirr, dsreb
        INTEGER(4), INTENT(IN) :: nrebmax, nder, ntet
        REAL(8), INTENT(OUT) :: tet(nder), phi(nder)

        REAL(8) :: dsreb_loc
        REAL(8) :: phil, tetl, tet1, rr2, rr32, dtet
        REAL(8) :: leI, leI1, Bl, B(3)
        REAL(8) :: x1(3), x2(3)
        REAL(8) :: stop_value

        INTEGER(4) :: Ilflag, Ifail, Iflag_I
        INTEGER(4) :: stop_type, hemi_flag
        INTEGER(4) :: I, J

        LOGICAL :: FLAG_IN_EARTH, EXCEEDED_NMAXREB
        INTEGER*4 depth

        REAL(8) :: pi,rad
        COMMON /rconst/rad,pi
        COMMON /flag_L/Ilflag

        dtet = pi/ntet

        ! "and one turns" -> "one shifts on surface in phi and one seeks teta 
        !    to have the constant IO and BO"
        ! et on tourne -> on se decale sur la surface en phi et on cherche teta
        ! pour avoir leI0 et B0 constants
        !

        ! compute the first foot position from input position (north side)
        hemi_flag = 1
        stop_type = 0
        stop_value = 1.0D0
        CALL FIND_FIELDLINE_FOOT(pos, dsreb, nrebmax, stop_type, stop_value, &
            hemi_flag, x1, B, Bl)

        rr2 = sqrt(x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
        tet(1) = ACOS(x1(3)/rr2)
        phi(1) = ATAN2(x1(2),x1(1))

        ! instead of going from point to foot, start at a test footpoint
        dsreb_loc = -dsreb
        DO I = 2,Nder
            
            phi(I) = phi(I-1)+2.D0*pi/Nder
            ! write(6,*) "phi(", I, ") = ", phi(I)


            ! first guess 
            Iflag_I = 0
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
            depth = 0
!
107         CONTINUE

            x1(1) = SIN(tetl)*COS(phi(I))
            x1(2) = SIN(tetl)*SIN(phi(I))
            x1(3) = COS(tetl)

            CALL INTEGRATE_J_ON_FIELDLINE_FROM_FOOTPOINT(x1, bmirr, dsreb_loc, nrebmax, rr2, leI, Ifail, EXCEEDED_NMAXREB)

            ! write(6,*) "tet", tetl, "leI-leI0", leI - leI0, ">nrebmax", exceeded_nmaxreb, "ifail", ifail

            ! first try for phi(I)
            IF (Iflag_I .EQ.0) THEN
                ! open field line -> try a lower tet value (so more north)
                IF (EXCEEDED_NMAXREB) THEN
                    tetl = tetl-dtet
                ELSE
                    tetl = tetl+dtet
                ENDIF
                leI1 = leI
                tet1 = tetl
                Iflag_I = 1
                GOTO 107
            ENDIF

            ! root finding exit condition
            IF ((leI-leI0)*(leI1-leI0) .LT. 0.D0) THEN
            !
            !       IF ((leI.eq.baddata.or.leI1.eq.baddata).and.(depth<10)) THEN
            !              dtet = dtet/2.D0
            !              depth = depth + 1
            !       ELSE
                            GOTO 108
            !       ENDIF
            ENDIF
    
            
            leI1 = leI
            tet1 = tetl
            ! leI < target -> field line too close -> try more north
            IF (leI.LT.leI0) THEN
                tetl = tetl-dtet
            ElSE
                tetl = tetl+dtet
            ENDIF

            ! exit condition
            IF (tetl.GT.pi .OR. tetl.LT.0.D0) GOTO 108
            GOTO 107
108         CONTINUE
    
            ! 
            dtet = pi/Ntet
            IF (leI.eq.baddata .or. leI1.eq.baddata) FLAG_IN_EARTH = 1 ! tet failed
            tet(I) = 0.5D0*(tetl+tet1)

            IF (EXCEEDED_NMAXREB .AND. leI.GT.0.D0) THEN
                Ilflag = 0
                RETURN
            ENDIF

            x1(1) = SIN(tet(I))*COS(phi(I))
            x1(2) = SIN(tet(I))*SIN(phi(I))
            x1(3) = COS(tet(I))
            CALL CHAMP(x1,B,Bl,Ifail)
            IF (Ifail.LT.0) THEN
                Ilflag = 0
                RETURN
            ENDIF
            IF (Bl.LT.bmirr) THEN
                Ilflag = 0
                RETURN
            ENDIF

!! Ne semble pas avoir d'incidence sur le calcul de L*
!! rajout D. Boscher 25 Mars 2009
!            Bmin3 = Bl
!            xmin3(1) = x1(1)
!            xmin3(2) = x1(2)
!            xmin3(3) = x1(3)
!            x13(1) = x1(1)
!            x13(2) = x1(2)
!            x13(3) = x1(3)
!            Iflag3 = 0
!            DO J = 1,Nrebmax
!                CALL sksyst(dsreb,x13,x3,Bl3,Ifail)
!                IF (Bl3 .GT. Bl) GOTO 203
!                IF (Bl3.LT.Bmin3) THEN
!                    xmin3(1) = x3(1)
!                    xmin3(2) = x3(2)
!                    xmin3(3) = x3(3)
!                    Bmin3 = Bl3
!                ENDIF
!                IF (Bl3 .LT. B0) Iflag3 = 1    !test pas a l'equateur
!                x13(1) = x3(1)
!                x13(2) = x3(2)
!                x13(3) = x3(3)
!            ENDDO
!203         CONTINUE
!            rr32 = xmin3(1)*xmin3(1)+xmin3(2)*xmin3(2)+xmin3(3)*xmin3(3)
!            IF (rr32 .LT. 1.03d0) GOTO 303
!            IF (Iflag3 .EQ. 1) THEN
!                x13(1) = x1(1)
!                x13(2) = x1(2)
!                x13(3) = x1(3)
!                Bl13 = Bl
!                DO J=1,Nrebmax
!                    CALL sksyst(dsreb,x13,x3,Bl3,Ifail)
!                    IF (Bl3 .LT. B0) GOTO 204
!                    x13(1) = x3(1)
!                    x13(2) = x3(2)
!                    x13(3) = x3(3)
!                    Bl13 = Bl3
!                ENDDO
!204             CONTINUE
!                ds3 = dsreb*(B0-Bl13)/(Bl3-Bl13)
!                CALL sksyst(ds3,x13,x3,Bl3,Ifail)
!                rr32  = x3(1)*x3(1)+x3(2)*x3(2)+x3(3)*x3(3)
!                IF (rr32 .LT. 1.03D0) GOTO 303
!        !
!                x13(1) = x1(1)
!                x13(2) = x1(2)
!                x13(3) = x1(3)
!                Bl13 = Bl
!                Iflag3 = 0
!                DO J=1,Nrebmax
!                    CALL sksyst(dsreb,x13,x3,Bl3,Ifail)
!                    IF (Bl3 .GT. Bl) GOTO 205
!                    IF (Bl3 .GT. B0 .AND. Iflag3 .EQ. 1) GOTO 205
!                    IF (Bl3 .LT. B0) Iflag3 = 1
!                    x13(1) = x3(1)
!                    x13(2) = x3(2)
!                    x13(3) = x3(3)
!                    Bl13 = Bl3
!                ENDDO
!205             CONTINUE
!                IF (J .EQ. Nrebmax+1) GOTO 303
!                ds3 = dsreb*(B0-Bl13)/(Bl3-Bl13)
!                CALL sksyst(ds3,x13,x3,Bl3,Ifail)
!                rr32  = x3(1)*x3(1)+x3(2)*x3(2)+x3(3)*x3(3)
!                IF (rr32 .LT. 1.03D0) GOTO 303
!            ENDIF
!! end rajout
!! Ne semble pas avoir d'incidence sur le calcul de L
!
        ENDDO
303     CONTINUE

        IF (rr32 .LT. 1.03d0) THEN
            Ilflag = 0
            RETURN
        ENDIF



    END SUBROUTINE

    SUBROUTINE INTEGRATE_FLUX_ON_POLAR_CAP(Nder, Ntet, tet, phi, flux)
        IMPLICIT NONE
        INTEGER(4), INTENT(IN) :: Nder, Ntet
        REAL(8), INTENT(IN) ::  tet(Nder), phi(Nder)
        REAL(8), INTENT(OUT) :: flux

        REAL(8) :: x1(3), B(3)
        REAL(8) :: Bl, BrR2, tetl, dtet
        INTEGER(4) :: Ifail, I, J

        INTEGER(4) :: Ilflag
        REAL(8) ::   Bo,xc,yc,zc,ct,st,cp,sp
        INTEGER(4) ::  k_ext,k_l,kint
        REAL*8     pi,rad
        
        COMMON /dipigrf/Bo,xc,yc,zc,ct,st,cp,sp
        COMMON /flag_L/Ilflag
        COMMON /magmod/k_ext,k_l,kint
        COMMON /rconst/rad,pi

        dtet = pi/Ntet

        x1(1) = 0.D0
        x1(2) = 0.D0
        x1(3) = 1.D0
        CALL CHAMP(x1,B,Bl,Ifail)
        IF (Ifail.LT.0)THEN
            Ilflag = 0
            RETURN
        ENDIF
        BrR2 = abs((x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3))) ! phi integrates B dot dA, or Br*R^2dphi*dtheta, R=1
        flux = BrR2*pi*dtet*dtet/4.D0

        DO I = 1,Nder
            tetl = 0.D0
            DO J = 1,Ntet
                tetl = tetl+dtet
                IF (tetl .GT. tet(I)) GOTO 111
                x1(1) = SIN(tetl)*COS(phi(I))
                x1(2) = SIN(tetl)*SIN(phi(I))
                x1(3) = COS(tetl)
                CALL CHAMP(x1,B,Bl,Ifail)
                IF (Ifail.LT.0)THEN
                    Ilflag = 0
                    RETURN
                ENDIF
                BrR2 = abs((x1(1)*B(1)+x1(2)*B(2)+x1(3)*B(3))) ! phi integrates B dot dA, or Br*R^2dphi*dtheta, R=1
                flux = flux + BrR2*SIN(tetl)*dtet*2.D0*pi/Nder
            ENDDO
111     CONTINUE
        ENDDO

    END SUBROUTINE


END MODULE