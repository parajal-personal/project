module subs_extrudate_swell_m

  use tfem_elem_m
  use math_defs_m

  implicit none

contains


! Special element routine for the surface tension force at the boundary point
! Used together with add_boundary_elements_point

  subroutine surface_tension_boundary_point1 ( mesh, problem, point, matrix, &
    vector, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: point
    logical, intent(in) :: matrix, vector
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer, parameter :: ndim = 2
    integer, parameter :: ndf = 3
    integer, parameter :: nodalp = 3
    integer, parameter :: ninti = 1

    real(dp) :: curvel(ninti), normal(ninti,ndim)
    real(dp) :: xig(ninti,1), x(nodalp,ndim)
    real(dp) :: phi(ninti,ndf), dphi(ninti,ndf,1), dxdxi(ninti,ndim)
    real(dp) :: normalc(ndim)

    integer :: curve, elem, coorsys
    real(dp) :: gammac


    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_boundary_point1:', &
        'No matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if

    xig(1,1) = -1._dp  ! assume first point

    ! surface tension
    curve = coefficients%i(501)
    elem = coefficients%i(502)
    coorsys = coefficients%i(23)
    gammac = coefficients%r(19)

    ! print *, 'curve =', curve, 'elem =', elem, 'coorsys =', coorsys, 'gammac =', gammac
    call shape_line_P2 ( xig, phi, dphi )

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )
    ! print *, 'x =', x
    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

!   outside normal to the outside "curve"

    normalc(1) =  normal(1,2)
    normalc(2) = -normal(1,1)
    ! print *, 'normalc1 =', normalc(1)
    ! print *, 'normalc2 =', normalc(2)
    ! print *, 'x(1, 2) =', x(1,2)
    elemvec = normalc * gammac
    ! print *, 'elemvec =', elemvec

    if ( coorsys == 1 ) then
!     axisymmetric
      elemvec = 2*pi*x(1,2)*elemvec
      ! print *, 'elemvec axisymmetric =', elemvec
    end if

  end subroutine surface_tension_boundary_point1


! Special element routine for the surface tension force at the boundary point
! Used together with buildsystem and object/onobjectnodes

  subroutine surface_tension_boundary_point2 ( mesh, problem, eleminfo, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: eleminfo
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer, parameter :: ndf = 9
    integer, parameter :: ninti = 1
    integer, parameter :: ndim = 2

    integer :: object

    real(dp), dimension(ndim) :: x1, x2, normalc
    real(dp) :: phi(ninti,ndf), xi(1,ndim)
    real(dp) :: d, gammac

!   test

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_boundary_point2:', &
        'No matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if

!   set object

    object = eleminfo%object

!   set coordinates (input)

    x1 = mesh%objects(object)%coor(1,:)
    x2 = coefficients%r(501:502)

!   shape function of fluid element

    xi(1,:) = mesh%objects(object)%refcoor(1,:)

    call shape_quad_Q2 ( xi, phi )

!   outside normal to the outside "curve"

    d = sqrt ( dot_product ( x1-x2, x1-x2 ) )

    normalc = ( x1 - x2 ) / d

!   surface tension coefficient

    gammac = coefficients%r(19)

!   element vector

    elemvec(    1:ndf  ) = gammac * normalc(1) * phi(1,:)
    elemvec(ndf+1:2*ndf) = gammac * normalc(2) * phi(1,:)

  end subroutine surface_tension_boundary_point2

end module subs_extrudate_swell_m
