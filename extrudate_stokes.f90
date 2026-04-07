! Newtonian extrudate swell problem
! Planar flow
! second-order Gear time integration for free surface
!
! Stripped-down version of extrudate_swell2d_c.f90 with all viscoelastic
! (constitutive equation, conformation tensor, DEVSS-G, log-conformation)
! parts removed.  Only Stokes + ALE free-surface remains.

program extrudate_swell2d_stokes

  use tfem_m
  use hsl_ma41_m
  use io_utils_m
  use surface_advection_elements_m
  use figplot_m
  use update_mesh_nodes_bc_m
  use stokes_elements_m
  use math_defs_m
  use metis5_m
  use subs_extrudate_swell_m

  implicit none

! constants flow problem

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gauss = 6,      & ! 6-point Gauss integration of triangles
    gaussb = 3,     & ! 3-point integration of boundary elements
    ndim = 2        ! dimension of space

! constants surface advection

  integer, parameter :: &
    ninti_sf_adv = 5,   & ! number of Gauss points
    vintpl = 6,         & ! P2 velocities
    method = 1            ! discretization method 0: Galerkin, 1: SUPG

! definitions – main problem

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln, solm1, sol_n
  type(sysvector_t) :: rhsd
  type(sample_t) :: sample_v
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients, coefficients_st
  type(solver_options_ma41_t) :: solver_options_u
  type(refinement_fields_t) :: refinement_fields
  type(vector_t) :: velocity, pressure

! 1D height function for surface advection

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh_sf_adv
  type(input_probdef_t) :: input_probdef_sf_adv
  type(problem_t) :: problem_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t), target :: rhsd_sf_adv, sol_sf_adv_pred, sol_sf_adv_pred_n
  type(oldvectors_t) :: oldvectors_sf_adv, oldvectors_sample_H, oldvectors_sf_adv_deriv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv, height_sf_adv
  type(subscript_t) :: hgt, hgt_end
  type(subscriptvec_t) :: velx_sf, vely_sf
  type(subscript_t) :: vel_in(2)
  type(subscript_t) :: velx7, vely7
  integer :: hintpl

! ALE mesh motion problem

  type(problem_t), target :: problem_lapl
  type(vector_t), target :: meshvel
  type(subscript_t) :: subsh
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial, xc
  real(dp), allocatable, dimension(:) :: Hhat, Hhatn

! Inlet problem (periodic Stokes channel)

  type(mesh_t) :: mesh_inlet
  type(input_probdef_t) :: input_probdef_inlet
  type(problem_t), target :: problem_inlet
  type(sysmatrix_t) :: sysmatrix_inlet
  type(sysvector_t), target :: sol_inlet, soln_inlet, solm1_inlet
  type(sysvector_t) :: rhsd_inlet
  type(oldvectors_t) :: oldvectors_ve_inlet, oldvectors_st
  type(coefficients_t) :: coefficients_inlet
  type(solver_options_ma41_t) :: solver_options_u_in
  real(dp), allocatable, dimension(:,:) :: vel_inlet

  integer :: physqvel, physqpress

! variables

  integer :: &
    timeint1 = 1,          & ! (first-order) Euler time integration (first step)
    timeint2 = 7,          & ! (second-order) semi-implicit Gear time integration
    numtimesteps = 1000,   & ! number of time steps
    coorsys = 1,           &      ! Cartesian coordinate system
    step0 = 0,             & ! initial step number
    vtkevery = 100           ! plot vtk every .. steps

  real(dp) :: &
    eta_s = 1.0_dp           ! total (Newtonian) viscosity

  real(dp), parameter :: &
    time0 = 0._dp,         & ! initial time
    beta = 1.0_dp,         & ! upwinding parameter in the SUPG method
    betai = 0.5_dp,        & ! upwinding parameter in the SUPG method free surface
    H = 1.0_dp,            & ! half height domain
    rmin = 0.1_dp,         & ! used in refinement fields
    rs_up = 2.5_dp,        & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up = 2.5_dp           ! integer_storage gradient-velocity-pressure LU (HSL)

  real(dp) :: U_avg = 1.0_dp       ! average velocity at the entry
  real(dp) :: flowrate             ! flowrate in half the channel (computed)

  real(dp) :: &
    dx_box = 0.4_dp,      & ! element spacing on the external boundaries
    dx_wall = 0.08_dp,    & ! element spacing on the upper boundary and die exit
    dx_inlet = 0.1_dp,    & ! element spacing on the inlet
    deltat = 1.e-2_dp       ! time step

  ! inlet geometry variables

  real(dp) :: &
    lx_in = 1.0_dp,      & ! size in x-direction
    ly_in = 1.0_dp,      & ! size in y-direction
    ox_in = -4.0_dp,     & ! x coordinate lower left corner
    oy_in = 0.0_dp,      & ! y coordinate lower left corner
    rs_up_in = 2.5_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up_in = 2.5_dp      ! integer_storage gradient-velocity-pressure LU (HSL)

  integer :: step, i, obj_ob, obj_sh, obj_inlet, surfimplct
  integer :: nnodes_inlet, ipost=0, ipost_in=0, iter
  integer :: ninner = 1          ! number of inner defect correction iterations

  real(dp), allocatable, dimension(:,:) :: coor
  real(dp) :: initial_h
  real(dp) :: ox=3.0_dp, L2=5.0_dp, oy=1.0_dp, gammac

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=300) :: filename
  character(len=300) :: curve4_fname
  character(len=32) :: s_eta, s_uavg

  logical :: lin_elem=.false.
  logical :: surface_tension = .false.

  namelist /comppar/ eta_s, U_avg, deltat, numtimesteps, &
    dx_box, dx_wall, dx_inlet, L2, gammac, vtkevery, surfimplct, coorsys

  call execute_command_line ( 'rm *.vtk' )

  ! call execute_command_line('rm -f mesh_inlet.geo mesh_inlet.msh mesh.geo mesh.msh' )

  read ( unit=*, nml=comppar )
  if ( gammac > 0.0_dp ) surface_tension = .true.
  print *, 'gammac =', gammac, '  surface_tension =', surface_tension
  print *, 'Newtonian extrudate swell'
  print *, 'eta_s =', eta_s, '  U_avg =', U_avg
  print *, 'Input mesh params: dx_box=', dx_box, ' dx_wall=', dx_wall, &
           ' dx_inlet=', dx_inlet

  if ( coorsys == 0 ) then
    flowrate = H*U_avg
  else
    flowrate = pi*H**2*U_avg
  end if
  print *, 'flowrate =', flowrate
! set physical quantity numbers (no DEVSS-G → only velocity + pressure)

  physqvel = 1
  physqpress = 2

! fill coefficients for main problem

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=600 )

  coefficients%i = 0

  coefficients%i(1) = uintpl      ! velocity interpolation
  coefficients%i(2) = pintpl      ! pressure interpolation
  coefficients%i(6) = physqvel    ! physical quantity velocity
  coefficients%i(7) = physqpress  ! physical quantity pressure
  coefficients%i(10) = gauss      ! Gauss quadrature
  coefficients%i(11) = gaussb     ! boundary Gauss quadrature
  coefficients%i(23) = coorsys    ! coordinate system
  coefficients%i(48) = 1          ! use mesh velocity for ALE formulation
  coefficients%i(60) = 1          ! use mesh velocity for ALE formulation

  coefficients%r = 0
  coefficients%r(1) = eta_s       ! viscosity
  coefficients%r(8) = deltat      ! time step
  coefficients%r(9) = beta        ! SUPG parameter
  coefficients%r(19) = gammac   ! surface tension coefficient

! fill coefficients for inlet channel

  call create_coefficients ( coefficients_inlet, ncoefi=450, ncoefr=400 )

  coefficients_inlet%i = 0

  coefficients_inlet%i(1)  = uintpl
  coefficients_inlet%i(2)  = pintpl
  coefficients_inlet%i(6)  = physqvel
  coefficients_inlet%i(7)  = physqpress
  coefficients_inlet%i(10) = gaussb
  coefficients_inlet%i(11) = gaussb
  coefficients_inlet%i(23) = coorsys
  coefficients_inlet%i(60) = 1 ! Separate vector for velocity gradient

  coefficients_inlet%r = 0

  coefficients_inlet%r(1:10) = &
    [ eta_s, 0._dp,   0._dp,   0._dp, 0._dp, &
       0._dp, 0._dp,  deltat,   beta,  0._dp &
      ]

  call create_coefficients ( coefficients_st, ncoefi=600, ncoefr=600 )

  coefficients_st%i = coefficients%i
  coefficients_st%r = 0._dp
  coefficients_st%r(1) = max(deltat, 1.0/gammac)
  coefficients_st%r(2) = gammac
  coefficients_st%i(11) = 0
  coefficients_st%i(12) = surfimplct
  print *, 'stab deltat', coefficients_st%r(1)
  call generate_read_mesh

  print *, 'nnodes', mesh%nnodes
  print *, 'nelem', mesh%nelem

  coefficients_inlet%r(6) = flowrate

  call define_inlet_problem
  print *, 'defined inlet problem'
  call define_problems
  print   *, 'defined main problem'
  call define_free_surface_problem
  print   *, 'defined free surface problem'
! create the structure oldvectors_ve for main problem

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nprob=1, nvec=1 )
  call create_vector ( problem, velocity, physq=physqvel )
  call create_vector ( problem, pressure, vec=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%v(1)%p => meshvel

  coefficients%i(22) = timeint1       ! first step integration scheme
  coefficients_inlet%i(22) = timeint1
  coefficients_sf_adv%i(5) = 1        ! start with first-order scheme

  call create_sysvector ( problem, sol_n )
  
  sol_n%u = 0._dp
  
  call create_oldvectors ( oldvectors_st, nsysvec=1, nprob=1 )
  
  oldvectors_st%s(1)%p => sol_n
  
  oldvectors_st%p(1)%p => problem

! ──────────────────────────── time loop ────────────────────────────

  do step = 1, numtimesteps

    if ( step == 2 ) then
      coefficients%i(22) = timeint2
      coefficients_inlet%i(22) = timeint2
      coefficients_sf_adv%i(5) = 2
    end if

    call solve_inlet_problem

    if ( step >= 2 ) then

      coefficients_sf_adv%i(5) = 2  ! second-order scheme

!     predict position of the surface and adapt mesh accordingly
      Hhatn = Hhat
      if ( .not. lin_elem ) then
        Hhat = 2._dp*sol_sf_adv_n%u(subsh%s) - sol_sf_adv_nm1%u(subsh%s)
      end if
      sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u(subsh%s) - sol_sf_adv_nm1%u(subsh%s)

      if ( lin_elem ) then

        call derive_vector ( mesh_sf_adv, problem_sf_adv, height_sf_adv, &
          elemsub=height_function_deriv, &
          coefficients=coefficients_sf_adv, oldvectors=oldvectors_sf_adv_deriv )

        Hhat = height_sf_adv%u

      end if

      meshcoor_nm1 = meshcoor_n
      meshcoor_n = mesh%coor

      call update_mesh_nodes_2D ( mesh, problem_lapl, disp=Hhat-Hhatn )

      call find_bounds_blocks ( mesh )

!     mesh velocity

      meshvel%u = reshape ( transpose ( &
          ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
                           [2*mesh%nnodes] )

    end if

!   sample velocity from inlet solution and use as boundary condition
!   at the inlet of the main domain

    nnodes_inlet = mesh%curves(6)%nnodes
    mesh_inlet%objects(obj_inlet)%coor(1:nnodes_inlet,1) = &
      mesh%coor(mesh%curves(6)%nodes,1)
    mesh_inlet%objects(obj_inlet)%coor(1:nnodes_inlet,2) = &
      mesh%coor(mesh%curves(6)%nodes,2)
    call fill_sample ( mesh_inlet, problem_inlet, sample_v, &
      ndegfd=ndim, object=obj_inlet, &
      elemsub=stokes_sample_velocity, coefficients=coefficients_inlet, &
      oldvectors=oldvectors_ve_inlet )
    vel_inlet(nnodes_inlet,:) = sample_v%u(nnodes_inlet,:)

!   fill solution vector with essential boundary conditions

    do i = 1, ndim
      sol%u(vel_in(i)%s) = sample_v%u(:,i)
    end do

!   inner iteration loop for defect correction convergence

    do iter = 1, ninner

!   build (assemble) matrix/vector for velocity/pressure problem (Stokes)

    call build_vpG

    if (surface_tension) then

      if ( coorsys == 1 .and. surface_tension ) then
        coefficients%r(65) = 0
        coefficients%r(24) = -gammac / sol_sf_adv_n%u(hgt_end%s(1))
        call add_boundary_elements ( mesh, problem, rhsd, curve=3, &
          physq=[physqvel], elemsub=stokes_natboun_normal, &
          coefficients=coefficients )
      end if

      if (surfimplct == 1) then
  !     build surface integral (implicit: adds matrix + vector)
        call add_boundary_elements ( mesh, problem,  rhsd, &
         sysmatrix = sysmatrix, elemsub=surface_tension_implicit_curve, &
          curve=4, coefficients=coefficients_st, physq=[physqvel], &
          oldvectors=oldvectors_st )
        else
        call add_boundary_elements ( mesh, problem, rhsd, &
          elemsub=surface_tension_curve, &
          curve=4, coefficients=coefficients, physq=[physqvel] )
        end if
      
      ! (die exit) to point 6 (outlet).
      coefficients%i(501) = 4
      coefficients%i(502) = 1

      call add_boundary_elements_point ( mesh, problem, rhsd, point=6, &
        elemsub=surface_tension_boundary_point1, &
        physq=(/physqvel/), coefficients=coefficients, buildmatrix=.false. )

    end if 
    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    if (step == 1 .and. iter == 1) then
      solver_options_u%pivot_order = 1
      solver_options_u%scaling = 1
      call renumber_metis_sysmatrix(sysmatrix)
    end if

    solver_options_u%real_storage=rs_up
    solver_options_u%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u )

!   update sol_n for next inner iteration
    sol_n%u = sol%u

    end do  ! inner iteration loop

    call extract_physvector ( mesh, problem, sol, velocity )

    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_ve )

!   solve surface advection (corrector)

    call solve_surface_height_corrector

    if (step == 1 ) then
      call postprocessing
      ipost = ipost + 1
    end if

!   write VTK files
    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing
        ipost = ipost + 1
      end if
    end if

!   copy old values

    call copy ( soln, solm1 )
    call copy ( sol, soln )

  end do

! ──────────────────────── end time loop ────────────────────────

!   Save final steady-state free-surface shape (curve 4)
  write(s_eta,  '(F8.4)') eta_s
  write(s_uavg, '(F8.1)') U_avg

  curve4_fname = 'curve4_stokes_' // trim(adjustl(s_eta)) // '_' // &
                 trim(adjustl(s_uavg)) // '.txt'
  open(unit=16, file=trim(curve4_fname), status='replace')
  write(16,'(A)') '# curve 4 final free-surface coordinates: x  y'
  write(16,'(A)') '# eta_s  = ' // trim(adjustl(s_eta))
  write(16,'(A)') '# U_avg  = ' // trim(adjustl(s_uavg))
  do i = 1, mesh%curves(4)%nnodes
    write(16,'(2ES24.16)') mesh%coor(mesh%curves(4)%nodes(i),1), &
                           mesh%coor(mesh%curves(4)%nodes(i),2)
  end do
  close(16)

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, soln, solm1, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors_ve )
  call delete ( meshvel )

  deallocate ( meshcoor_initial, meshcoor_n, meshcoor_nm1 )

  call delete ( mesh_sf_adv )
  call delete ( problem_sf_adv )
  call delete ( input_probdef_sf_adv )
  call delete ( sol_sf_adv, rhsd_sf_adv )
  call delete ( sol_sf_adv_pred, sol_sf_adv_n, sol_sf_adv_nm1 )
  call delete ( sysmatrix_sf_adv )
  call delete ( oldvectors_sf_adv )
  call delete ( coefficients_sf_adv )
  call delete ( hgt, hgt_end )
  call delete ( problem_lapl )

  call delete ( mesh_inlet )
  call delete ( problem_inlet )
  call delete ( input_probdef_inlet )
  call delete ( sol_inlet, soln_inlet, solm1_inlet, rhsd_inlet )
  call delete ( sysmatrix_inlet )
  call delete ( coefficients_inlet )
  call delete ( oldvectors_ve_inlet )

contains

! ═══════════════════════════════════════════════════════════════════
! generate and read mesh
! ═══════════════════════════════════════════════════════════════════

 subroutine generate_read_mesh

    integer :: i

    call write_gmsh_parameters ( ox, oy, L2, dx_box, dx_wall, dx_inlet )

    call execute_command_line ( 'gmsh -2 -order 2 -algo front2d -o mesh.msh &
                  &mesh.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=2, &
      physgeom=.true. )

    call add_to_mesh ( mesh, object='curve', objectcurve=6, topology=.true., &
      intrule=3 )
    obj_ob = mesh%nobjects

    call add_to_mesh ( mesh, curve=[4,5] ) ! curve7
    
    call fill_mesh_parts ( mesh )

!   define some arrays for the ALE mesh position at old times
    allocate ( meshcoor_n(mesh%nnodes,mesh%ndim), &
      meshcoor_nm1(mesh%nnodes,mesh%ndim) )

    allocate ( meshcoor_initial(mesh%nnodes,2) )
      meshcoor_initial = mesh%coor

    meshcoor_n = mesh%coor
    meshcoor_nm1 = mesh%coor

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  end subroutine generate_read_mesh

! write the parameters in gmsh format

  subroutine write_gmsh_parameters ( ox, oy, L2, dx_box, dx_wall, dx_inlet )

    real(dp) :: ox, oy, L2, dx_box, dx_wall, dx_inlet
    integer, parameter :: ncurve4_points = 200
    integer :: icurve
    real(dp), dimension(ncurve4_points,2) :: refinement_coor

!   points on curve 4 from the die exit (x=0) to the outlet (x=L2)
!   choose 50 uniformly spaced points in [0,1]
    do icurve = 1, ncurve4_points
      refinement_coor(icurve,1) = L2 * real(icurve - 1, dp) / real(ncurve4_points - 1, dp)
    end do

    refinement_coor(:,2) = H

    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,H],[1,2]), &
      distmin=2*rmin, distmax=4*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=0.05_dp, distmax=4*rmin, dx_fine=dx_wall, dx_coarse=dx_box )

    open ( unit=25, file='mesh.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', ox, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy = ', oy, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'L2 = ', L2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'

    write ( 25, '(1X,A,F18.14,A)' ) 'dx_wall = ', dx_wall, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_inlet = ', dx_inlet, ';'

!   write the refinement fields (the file needs to be open for writing)
    call write_refinement_fields ( refinement_fields, 'mesh.geo' )

    write ( 25, '(/1x,a)' ) 'Include "mesh_2D_bc.igo";'
    close ( 25 )

!   refinement points served their purpose
    call delete_refinement_fields ( refinement_fields )

  end subroutine write_gmsh_parameters

! ═══════════════════════════════════════════════════════════════════
! build Stokes velocity/pressure system (main domain)
! ═══════════════════════════════════════════════════════════════════

  subroutine build_vpG

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[physqvel,physqpress], &
      physqcol=[physqvel,physqpress], coefficients=coefficients )

  end subroutine build_vpG

! ═══════════════════════════════════════════════════════════════════
! solve convection equation for the surface height (corrector)
! ═══════════════════════════════════════════════════════════════════

  subroutine solve_surface_height_corrector

    use postprocessing_m

    real(dp) :: max_height, end_height
    type(solver_options_ma41_t) :: solver_options_h

    velocity_sf_adv%u(velx_sf%s) = sol%u(velx7%s)
    velocity_sf_adv%u(vely_sf%s) = sol%u(vely7%s)

    call build_system ( mesh_sf_adv, problem_sf_adv, sysmatrix_sf_adv, &
      rhsd_sf_adv, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_sf_adv, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_sf_adv )

    call add_effect_of_essential_to_rhs ( problem_sf_adv, sysmatrix_sf_adv, &
      sol_sf_adv, rhsd_sf_adv )

!   MA41 solver storage

    if (step == 1) then
      solver_options_h%pivot_order = 1
      solver_options_h%scaling = 1
      call renumber_metis_sysmatrix(sysmatrix_sf_adv)
    end if

    solver_options_h%integer_storage = 2.5
    solver_options_h%real_storage    = 2.5

!   solve system

    call solve_system_ma41 ( sysmatrix_sf_adv, rhsd_sf_adv, sol_sf_adv, &
                             solver_options=solver_options_h )

    call copy ( sol_sf_adv_n, sol_sf_adv_nm1 )
    call copy ( sol_sf_adv, sol_sf_adv_n )

!   compute the maximum, minimum and end radius

    max_height = maxval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1))

    print '(i6,4es16.8)', step0+step, time0+step*deltat, &
                            max_height, end_height, max_height/initial_h

  end subroutine solve_surface_height_corrector

! ═══════════════════════════════════════════════════════════════════
! postprocessing – VTK output (main domain)
! ═══════════════════════════════════════════════════════════════════

  subroutine postprocessing

    type(oldvectors_t) :: oldvectors_dve
    type(vector_t) :: pressure_pp, eff_shear

    if ( .not. mesh%meshparts ) call fill_mesh_parts ( mesh )

    call create_oldvectors ( oldvectors_dve, nsysvec=1 )
    oldvectors_dve%s(1)%p => sol

    call create_vector ( problem, pressure_pp, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure_pp, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    call create_vector ( problem, eff_shear, vec=3 )

    coefficients%i(13) = 8

!   derive the effective shear rate in all nodes
    call derive_vector ( mesh, problem, eff_shear, &
      elemsub=stokes_deriv, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    write(filename,'(a,i4.4,a)') 'flow_stokes', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure_pp, &
      dataname='pressure', filename=filename )

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    ! write(filename,'(a,i4.4,a)') 'shearrate_stokes', ipost, '.vtk'
    ! call write_scalar_vtk ( mesh, problem, vector=eff_shear, &
    !   dataname='effective_shearrate', filename=filename )

    call delete ( pressure_pp )
    call delete ( eff_shear )
    call delete ( oldvectors_dve )

  end subroutine postprocessing

! ═══════════════════════════════════════════════════════════════════
! define inlet problem (periodic Stokes channel)
! ═══════════════════════════════════════════════════════════════════

  subroutine define_inlet_problem

!   create mesh for channel

    open ( unit=25, file='mesh_inlet.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'ox_in = ', ox_in, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy_in = ', oy_in, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'lx_in = ', lx_in, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'ly_in = ', ly_in, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_wall = ', dx_wall, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_inlet = ', dx_inlet, ';'

    write ( 25, '(/1x,a)' ) 'Include "mesh_inlet_2D_bc.igo";'
    close ( 25 )

    call execute_command_line ( 'gmsh -2 -format msh2 -order 2 -algo front2d &
                &-o mesh_inlet.msh mesh_inlet.geo -v 0' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh_inlet, filename='mesh_inlet.msh', ndim=2, &
      physgeom=.true. )

    call add_to_mesh ( mesh_inlet, matchingcurve=[4,2], replace=4, &
      displacement=[-1._dp,0._dp] )

    call add_to_mesh ( mesh_inlet, curve=(/-4/) )     ! curve 5

    call fill_mesh_parts ( mesh_inlet )

!   problem definition of velocity/pressure for the channel

    call create_input_probdef ( mesh_inlet, input_probdef_inlet, nvec=3, nphysq=2 )

    input_probdef_inlet%vec_elementdof(1)%a =   &
        reshape ( (/ 2,2,2,2,2,2,    &  ! velocity
                     1,0,1,0,1,0,    &  ! pressure
                     1,1,1,1,1,1 /), &  ! scalar
                     (/6,3/) )

    input_probdef_inlet%physq = (/physqvel, physqpress/)
    input_probdef_inlet%probnr = 1

!   Dirichlet boundary conditions

!   velocity on lower boundary
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=1, physq=physqvel, &
      degfd=[0,1], excludecurves=(/4/) )
!   velocity on upper boundary
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=3, physq=physqvel, &
      excludecurves=(/4/) )
!   pressure in point 1
    call define_essential ( mesh_inlet, input_probdef_inlet, point=1, physq=physqpress )

!   constraints for periodical boundary conditions

!   velocities
    call define_constraint ( mesh_inlet, input_probdef_inlet, &
      physq=physqvel, curve1=2, curve2=4, discretization='collocation')

!   constraint for the flow rate
    call define_constraint ( mesh_inlet, input_probdef_inlet, &
      physq=physqvel, curve1=4, nglobalc=1 )

    call problem_definition ( input_probdef_inlet, mesh_inlet, problem_inlet )

!   create system vectors for velocity/pressure
!   (solution and right-hand side)

    call create_sysvector ( problem_inlet, sol_inlet, soln_inlet, solm1_inlet )
    call create_sysvector ( problem_inlet, rhsd_inlet )

    sol_inlet%u = 0._dp
    solm1_inlet%u = 0._dp

!   fill solution vector with essential boundary conditions

!   set y velocity = 0 on axis
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      curve1=1, physq=physqvel, degfd=2, value=0._dp )
!   set velocity = 0 on upper boundary
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      curve1=3, physq=physqvel, degfd=1, value=0._dp )
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      curve1=3, physq=physqvel, degfd=2, value=0._dp )
!   set pressure level = 0 in lower left corner
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      point=1, physq=physqpress, value=0._dp )

!   create system matrix of velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix_inlet, mesh_inlet, problem_inlet )
    call create_sysmatrix_structure_constraint ( sysmatrix_inlet, mesh_inlet, problem_inlet )
    call finalize_sysmatrix_structure ( sysmatrix_inlet )

    call create_sysmatrix_data ( sysmatrix_inlet )

!   create oldvectors for the inlet

    call create_oldvectors ( oldvectors_ve_inlet, nsysvec=1, nprob=1 )

    oldvectors_ve_inlet%s(1)%p => sol_inlet
    oldvectors_ve_inlet%p(1)%p => problem_inlet

  end subroutine define_inlet_problem

! ═══════════════════════════════════════════════════════════════════
! solve inlet problem (Stokes with periodic BC)
! ═══════════════════════════════════════════════════════════════════

  subroutine solve_inlet_problem

!   build (assemble) matrix/vector for velocity/pressure problem

    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, &
      elemsub=stokes_elem, physqrow=(/physqvel,physqpress/), &
      physqcol=(/physqvel,physqpress/), coefficients=coefficients_inlet )

!   periodical condition on velocities

    call build_system_constraint ( mesh_inlet, problem_inlet, &
      sysmatrix_inlet, rhsd_inlet, &
      constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients_inlet )

!   imposed flow rate

    call build_system_constraint ( mesh_inlet, problem_inlet, &
      sysmatrix_inlet, rhsd_inlet, &
      constraint1=2, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients_inlet )

    call add_effect_of_essential_to_rhs ( problem_inlet, sysmatrix_inlet, &
      sol_inlet, rhsd_inlet )

!   solve velocity/pressure problem for channel

    solver_options_u_in%real_storage=rs_up_in
    solver_options_u_in%integer_storage=is_up_in

    call solve_system_ma41 ( sysmatrix_inlet, rhsd_inlet, sol_inlet, &
      solver_options=solver_options_u_in )

    if (step == 1 ) then
      call postprocessing_inlet
      ipost_in = ipost_in + 1
    end if

!   write VTK files
    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing_inlet
        ipost_in = ipost_in + 1
      end if
    end if

!   copy old values

    call copy ( soln_inlet, solm1_inlet )
    call copy ( sol_inlet, soln_inlet )

  end subroutine solve_inlet_problem

! ═══════════════════════════════════════════════════════════════════
! postprocessing – VTK output (inlet)
! ═══════════════════════════════════════════════════════════════════

  subroutine postprocessing_inlet

    type(oldvectors_t) :: oldvectors_dve
    type(vector_t) :: pressure_in

    if ( .not. mesh_inlet%meshparts ) call fill_mesh_parts ( mesh_inlet )

    call create_oldvectors ( oldvectors_dve, nsysvec=1 )
    oldvectors_dve%s(1)%p => sol_inlet

    call create_vector ( problem_inlet, pressure_in, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh_inlet, problem_inlet, pressure_in, &
      elemsub=stokes_pressure, coefficients=coefficients_inlet, &
      oldvectors=oldvectors_dve )

    write(filename,'(a,i4.4,a)') 'flow_in', ipost_in, '.vtk'
    call write_scalar_vtk ( mesh_inlet, problem_inlet, vector=pressure_in, &
      dataname='pressure', filename=filename )

    call write_vector_vtk ( mesh_inlet, problem_inlet, filename=filename, &
      dataname='velocity', sysvector=sol_inlet, physq=physqvel, &
      append=.true. )

    call delete ( pressure_in )
    call delete ( oldvectors_dve )

  end subroutine postprocessing_inlet

! ═══════════════════════════════════════════════════════════════════
! define main problem (velocity + pressure)
! ═══════════════════════════════════════════════════════════════════

  subroutine define_problems

!   add an object in the channel mesh to sample the velocity
!   at the inlet of the extrudate swell problem

    allocate ( xc(mesh%curves(6)%nnodes,2) )

    xc(:,1) = mesh%coor(mesh%curves(6)%nodes,1)
    xc(:,2) = mesh%coor(mesh%curves(6)%nodes,2)

    warn_add_to_mesh_after_meshgen_parts = .false.
    call add_to_mesh ( mesh_inlet, object='coordinates', coor=xc )
    obj_inlet = mesh_inlet%nobjects
    call fill_mesh_parts_objects ( mesh_inlet, object1=obj_inlet )
    warn_add_to_mesh_after_meshgen_parts = .true.

    deallocate ( xc )

!   problem definition of velocity/pressure

    call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

    input_probdef%vec_elementdof(1)%a =  &
        reshape ( [ 2,2,2,2,2,2,   &  ! velocity
                    1,0,1,0,1,0,   &  ! pressure
                    1,1,1,1,1,1 ], &  ! scalar for postprocessing
                  [6,3] )

    input_probdef%physq = [physqvel, physqpress]
    input_probdef%probnr = 7

!   Dirichlet boundary conditions

!   inlet
    call define_essential ( mesh, input_probdef, &
      curve1=6, physq=physqvel )
!   center line
    call define_essential ( mesh, input_probdef, &
      curve1=1, physq=physqvel, degfd=[0,1] )
    call define_essential ( mesh, input_probdef, &
      curve1=2, physq=physqvel, degfd=[0,1] )
!   outflow: uy=0
    call define_essential ( mesh, input_probdef, &
      curve1=3, physq=physqvel, degfd=[0,1] )
!   wall
    call define_essential ( mesh, input_probdef, &
      curve1=5, physq=physqvel )

    call problem_definition ( input_probdef, mesh, problem )

!   subscript for the inlet BC

    do i = 1, ndim
      call create_subscript ( mesh, problem, vel_in(i), curves=[6], &
        physqarr=[physqvel], degfd=i )
    end do

    call create_subscript ( mesh, problem, velx7, physqarr=[physqvel], &
       degfd=1, curves=(/4/), fillnodes=.true. )
    call create_subscript ( mesh, problem, vely7, physqarr=[physqvel], &
       degfd=2, curves=(/4/), fillnodes=.true. )

    allocate ( vel_inlet(mesh%curves(6)%nnodes,2) )

!   create system vectors for velocity/pressure
!   (solution and right-hand side)

    call create_sysvector ( problem, sol, soln, solm1, rhsd )

!   fill solution vector with essential boundary conditions

    sol%u = 0._dp
    solm1%u = 0._dp

!   create system matrix of velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   define some arrays and vector for the ALE mesh movement

    call create_vector ( problem, meshvel, physq=physqvel )
    meshvel%u = 0._dp  ! initialize

  end subroutine define_problems

! ═══════════════════════════════════════════════════════════════════
! define free surface problem (identical to viscoelastic version)
! ═══════════════════════════════════════════════════════════════════

  subroutine define_free_surface_problem

    if ( lin_elem ) then
      hintpl = 2
    else
      hintpl = 6
    end if

!   fill coefficients for surface advection problem

    call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

    coefficients_sf_adv%i(1:4) = [ ninti_sf_adv, 2, 0, method ]
    coefficients_sf_adv%i(5:) = 0
    coefficients_sf_adv%i(5) = 1 ! time-integration
    coefficients_sf_adv%i(6) = hintpl ! height interpolation
    coefficients_sf_adv%i(9) = 3 ! numerical table for Gauss

    coefficients_sf_adv%r = 0

    coefficients_sf_adv%r(4) = deltat
    coefficients_sf_adv%r(5) = betai

    if ( lin_elem ) then
      coefficients_sf_adv%i(10) = vintpl
      coefficients_sf_adv%i(12) = 1 ! use shapefunc for velocity for element shape
    end if

!   create mesh for surface advection

    mesh_options%elshape = 2   ! three-node line elements

    mesh_options%nx = mesh%curves(4)%nelem

    call line1d ( mesh_sf_adv, mesh_options )

!   set coordinates

    do i = 1, mesh%curves(4)%nnodes
      mesh_sf_adv%coor(i,1) = mesh%coor(mesh%curves(4)%nodes(i),1)
    end do

    call fill_mesh_parts ( mesh_sf_adv )

!   add object for sampling velocity in height function advection

    initial_h = mesh%coor(mesh%curves(4)%nodes(1),2)

    allocate ( coor(mesh_sf_adv%nnodes,2) )

    coor(:,1) = mesh_sf_adv%coor(1:mesh_sf_adv%nnodes,1)
    coor(:,2) = initial_h

    warn_add_to_mesh_after_meshgen_parts = .false.
    call add_to_mesh ( mesh, object='coordinates', coor=coor )
    obj_sh = mesh%nobjects
    call fill_mesh_parts_objects ( mesh, object1=obj_sh )
    warn_add_to_mesh_after_meshgen_parts = .true.

    deallocate ( coor )

!   allocate arrays for the ALE displacement problem

    allocate ( Hhat(mesh_sf_adv%nnodes), Hhatn(mesh_sf_adv%nnodes) )
    Hhat = initial_h

!   problem definition for surface advection

    if ( lin_elem ) then
      call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=3, &
        nphysq=1 )

      input_probdef_sf_adv%vec_elementdof(1)%a = &
         reshape ( [ 1,0,1,    &  ! height function
                  2,2,2,    &  ! velocity
                  1,1,1 ],  &  ! height in all nodes
                 [3,3] )

      input_probdef_sf_adv%physq = [1]
      input_probdef_sf_adv%probnr = 10

    else
      call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=2, &
        nphysq=1 )

      input_probdef_sf_adv%vec_elementdof(1)%a = &
          reshape ( [ 1,1,1,    &  ! height function
                      2,2,2 ],  &  ! velocity
                     [3,2] )

      input_probdef_sf_adv%physq = [1]
      input_probdef_sf_adv%probnr = 9
    end if

    call define_essential ( mesh_sf_adv, input_probdef_sf_adv, point=2, physq=1 )

    print *, 'P1 BC',  mesh_sf_adv%coor(mesh_sf_adv%points(1),:)
    print *, 'P2 BC',  mesh_sf_adv%coor(mesh_sf_adv%points(2),:)

!   define problem

    call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

    call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
    call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )
    call create_sysvector ( problem_sf_adv, sol_sf_adv_pred, sol_sf_adv_pred_n )

!   fill solution vector with essential boundary conditions

    sol_sf_adv%u = initial_h
    sol_sf_adv_n%u = sol_sf_adv%u
    sol_sf_adv_nm1%u = sol_sf_adv_n%u

    call fill_sysvector ( mesh_sf_adv, problem_sf_adv, sol_sf_adv, point=2, &
      physq=1, value=initial_h )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_sf_adv, mesh_sf_adv, &
      problem_sf_adv )

    call create_sysmatrix_data ( sysmatrix_sf_adv )

!   create vectors
    call create_vector ( problem_sf_adv, velocity_sf_adv, vec=2 )
    if ( lin_elem ) then
      call create_vector ( problem_sf_adv, height_sf_adv, vec=3 )
    end if

!   subscripts for velocity

    call create_subscript ( mesh_sf_adv, problem_sf_adv, velx_sf, vec=2, degfd=1 )
    call create_subscript ( mesh_sf_adv, problem_sf_adv, vely_sf, vec=2, degfd=2 )

    if ( lin_elem ) then
      call create_oldvectors ( oldvectors_sf_adv_deriv, nsysvec=1 )
      oldvectors_sf_adv_deriv%s(1)%p => sol_sf_adv_pred ! prediction
    end if

    call create_oldvectors ( oldvectors_sf_adv, nsysvec=2, nvec=1 )

    oldvectors_sf_adv%s(1)%p => sol_sf_adv_n    ! corrector at n
    oldvectors_sf_adv%s(2)%p => sol_sf_adv_nm1  ! corrector at nm1

    oldvectors_sf_adv%v(1)%p => velocity_sf_adv ! advection velocity at np1

    call create_oldvectors ( oldvectors_sample_H, nsysvec=1 )
    oldvectors_sample_H%s(1)%p => sol_sf_adv_pred

!   create subscript for the height values

    call create ( mesh_sf_adv, problem_sf_adv, subsh )
    call create ( mesh_sf_adv, problem_sf_adv, hgt )
    call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[1] )

  end subroutine define_free_surface_problem

 subroutine surface_tension_implicit_curve ( mesh, problem, curve, elem, &
  matrix, vector, first, last, coef_st, oldvectors, elemmat, elemvec )

    use stokes_globals_m
! 
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coef_st
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: N, i, j, ip, nc, surfimpl
    integer :: nodes(size(mesh%curves(curve)%topology, 1))
    real(dp) :: deltat, gammac
    real(dp), save, allocatable :: xold(:,:), xng(:,:), dxdxi_n(:,:), g11_up(:)
    real(dp), save, allocatable :: vel_old(:), curvel_n(:), g11_up_n(:)

    if ( first ) then

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve, &
        maxvel3D=1 )

      allocate ( wg(ninti), curvel(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )
      allocate ( g1_up(ninti,ndim) )
      allocate ( work(ninti) )

      allocate ( xold(nodalp,ndim), xng(ninti,ndim) )
      allocate ( dxdxi_n(ninti,ndim), g11_up(ninti) )
      allocate ( vel_old(size(elemvec)) )
      allocate ( curvel_n(ninti), g11_up_n(ninti) )

      call set_Gauss_integration ( gauss, xig, wg )
      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      g1_up=g1_up )

    call isoparametric_coordinates ( x, phi, xg )

    g11_up = sum ( g1_up**2, dim=2 )

    nodes = mesh%curves(curve)%topology(:,elem,2)
    xold = meshcoor_n(nodes,:)

    call isoparametric_deformation_curve ( xold, dphi(:,:,1), dxdxi_n, work )
    call isoparametric_coordinates ( xold, phi, xng )

    curvel_n = work
    g11_up_n = 1._dp / (curvel_n**2)

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
      curvel_n = 2 * pi * xng(:,2) * curvel_n
    end if

    deltat = coef_st%r(1)
    surfimpl = coef_st%i(12)
    gammac = coef_st%r(2)

    if ( matrix ) then

      elemmat = 0._dp

      if ( surfimpl == 1 ) then
        do i = 1, ndf
          do j = i, ndf
            work = dphi(:,i,1) * g11_up_n * dphi(:,j,1)
            elemmat(i,j) = sum ( work * curvel_n * wg )
            elemmat(j,i) = elemmat(i,j)
          end do
        end do

        do i = 1 + ndf, 2 * ndf
          do j = i, 2 * ndf
            work = dphi(:,i-ndf,1) * g11_up_n * dphi(:,j-ndf,1)
            if ( coorsys == 1 ) then
              work = work + phi(:,i-ndf) * phi(:,j-ndf) / xng(:,2)**2
            end if
            elemmat(i,j) = sum ( work * curvel_n * wg )
            elemmat(j,i) = elemmat(i,j)
          end do
        end do

        elemmat = gammac * deltat * elemmat

      end if

    end if

    if ( vector ) then

      do i = 1, ndim
        do N = 1, ndf
          if ( surfimpl == 1 ) then
            work = dphi(:,N,1) * g11_up_n * dxdxi_n(:,i)
            tmp(N,i) = sum ( work * curvel_n * wg )
          else
            work = dphi(:,N,1) * g1_up(:,i)
            tmp(N,i) = sum ( work * curvel * wg )
          end if
        end do
      end do

      if ( coorsys == 1 ) then
        do N = 1, ndf
          if ( surfimpl == 1 ) then
            work = phi(:,N) / xng(:,2)
            tmp(N,2) = tmp(N,2) + sum ( work * curvel_n * wg )
          else
            work = phi(:,N) / xg(:,2)
            tmp(N,2) = tmp(N,2) + sum ( work * curvel * wg )
          end if
        end do
      end if

      tmp = -gammac * tmp

      nc = ndf * ndim

      if ( vel3D == 0 ) then
        elemvec = reshape ( tmp, [nc] )
      else
        elemvec(1:nc) = reshape ( tmp, [nc] )
        elemvec(nc+1:) = 0._dp
      end if

    end if

!   defect correction: add + gamma*dt*K*v_n to cancel O(dt) steady-state error

    if ( vector .and. surfimpl == 1 ) then
      call get_sysvector_geometry ( mesh, oldvectors%p(1)%p, &
        oldvectors%s(1)%p, elem, vel_old, curve=curve, physq=[physqvel] )
      elemvec = elemvec + matmul ( elemmat, vel_old )
    end if

    if ( last ) then

      deallocate ( wg, curvel )
      deallocate ( xig, phi, x )
      deallocate ( tmp )
      deallocate ( xg, xng )
      deallocate ( dphi, dxdxi, dxdxi_n )
      deallocate ( g1_up, g11_up )
      deallocate ( work )
      deallocate ( xold )
      deallocate ( vel_old )
      deallocate ( curvel_n, g11_up_n )

    end if

  end subroutine surface_tension_implicit_curve

end program extrudate_swell2d_stokes
