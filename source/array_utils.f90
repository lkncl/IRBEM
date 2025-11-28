!=====================================================
! Created by L. Nicolle 2025
! check_array_sizes.f90
! Utility module to check ensure shape matching for F90 dynamic array inputs
!=====================================================

MODULE array_utils
  IMPLICIT NONE
CONTAINS

  !-----------------------------------------
  ! Vérifie un tableau REAL*8 jusqu'à 3D
  !-----------------------------------------
  SUBROUTINE check_dims_real8_1d(arr, dims_expected, name)
    REAL*8, INTENT(IN) :: arr(:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 1) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 1"
       STOP
    END IF

    IF (SIZE(arr) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr)
      STOP
    END IF
  END SUBROUTINE check_dims_real8_1d

  SUBROUTINE check_dims_real8_2d(arr, dims_expected, name)
    REAL*8, INTENT(IN) :: arr(:,:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 2) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 2"
       STOP
    END IF

    IF (SIZE(arr,1) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr,1)
      STOP
    END IF
    IF (SIZE(arr,2) /= dims_expected(2)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim2 wrong. Expected ", dims_expected(2), " got ", SIZE(arr,2)
      STOP
    END IF
  END SUBROUTINE check_dims_real8_2d

  SUBROUTINE check_dims_real8_3d(arr, dims_expected, name)
    REAL*8, INTENT(IN) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 3) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 3"
       STOP
    END IF

    write(6,*) dims_expected
    write(6,*) SIZE(arr, 1), SIZE(arr,2), SIZE(arr,3)

    IF (SIZE(arr,1) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr,1)
      STOP
    END IF
    IF (SIZE(arr,2) /= dims_expected(2)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim2 wrong. Expected ", dims_expected(2), " got ", SIZE(arr,2)
      STOP
    END IF
    IF (SIZE(arr,3) /= dims_expected(3)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim3 wrong. Expected ", dims_expected(3), " got ", SIZE(arr,3)
      STOP
    END IF
  END SUBROUTINE check_dims_real8_3d


  !-----------------------------------------
  ! Vérifie un tableau INTEGER jusqu'à 3D
  !-----------------------------------------
  SUBROUTINE check_dims_int_1d(arr, dims_expected, name)
    INTEGER, INTENT(IN) :: arr(:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 1) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 1"
       STOP
    END IF

    IF (SIZE(arr) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr)
      STOP
    END IF
  END SUBROUTINE check_dims_int_1d

  SUBROUTINE check_dims_int_2d(arr, dims_expected, name)
    INTEGER, INTENT(IN) :: arr(:,:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 2) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 2"
       STOP
    END IF

    IF (SIZE(arr,1) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr,1)
      STOP
    END IF
    IF (SIZE(arr,2) /= dims_expected(2)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim2 wrong. Expected ", dims_expected(2), " got ", SIZE(arr,2)
      STOP
    END IF
  END SUBROUTINE check_dims_int_2d

  SUBROUTINE check_dims_int_3d(arr, dims_expected, name)
    INTEGER, INTENT(IN) :: arr(:,:,:)
    INTEGER, INTENT(IN) :: dims_expected(:)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: ndims

    ndims = SIZE(dims_expected)
    IF (ndims /= 3) THEN
       WRITE(*,*) "F90: dims_expected must be an array of size 3"
       STOP
    END IF

    IF (SIZE(arr,1) /= dims_expected(1)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim1 wrong. Expected ", dims_expected(1), " got ", SIZE(arr,1)
      STOP
    END IF
    IF (SIZE(arr,2) /= dims_expected(2)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim2 wrong. Expected ", dims_expected(2), " got ", SIZE(arr,2)
      STOP
    END IF
    IF (SIZE(arr,3) /= dims_expected(3)) THEN
      WRITE(*,*) "F90: ", TRIM(name), " dim3 wrong. Expected ", dims_expected(3), " got ", SIZE(arr,3)
      STOP
    END IF
  END SUBROUTINE check_dims_int_3d

END MODULE array_utils