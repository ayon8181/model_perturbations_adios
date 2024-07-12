program single_gaussian_perturb_for_psf

! Ebru
! Originally written in Nice, 2015
! Updated in Santa Barbara, July 2018
! Ridvan
! Updated in CSM, 2022

  use adios_helpers_mod
  use manager_adios

  implicit none

  include 'constants.h'
  include 'values_from_mesher.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NGLOB=NGLOB_CRUST_MANTLE
  integer,parameter:: NKERNEL=7    !vpv,vph,vsv,vsh,rho,eta,qmu

  integer, parameter :: NSPEC_MAX = NSPEC_CRUST_MANTLE_ADJOINT
  integer, parameter :: NGLOB_MAX = NGLOB_CRUST_MANTLE
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x, y, z
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: solver_point
  real(kind=CUSTOM_REAL) :: max_new
  character(len=512):: output_file,solver_file
  character(len=512):: output_bp_file,solver_bp_file
  character(len=150):: solver_name(3)
  character(len=512):: latt,lonn,rr,ssigma, update_params
  character(len=150):: kernel_name(NKERNEL)

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL):: new_kernel

  real(kind=CUSTOM_REAL) :: xloc,yloc,zloc,sigma,exp_val
  real(kind=CUSTOM_REAL) :: lat,lon,r,phi,theta

  ! ADIOS variables
  integer(kind=8) :: local_dim, group_size_inc

  ! MPI
  integer:: myrank, sizeprocs

  ! kernels
  integer:: i, j, k, ispec, iglob, iker

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  call initialize_adios()

! read in parameters / directories of mesh files & output kernel
  call getarg(1,latt)
  call getarg(2,lonn)
  call getarg(3,rr)
  call getarg(4,ssigma)
  call getarg(5,solver_file)
  call getarg(6,output_file)
  call getarg(7,update_params)
  if (trim(latt) == '' &
    .or. trim(lonn) == '' &
    .or. trim(rr) == '' &
    .or. trim(ssigma) == '' &
    .or. trim(solver_file) == '' &
    .or. trim(output_file) == '' &
    .or. trim(update_params) == '') then
      call exit_MPI(myrank, 'Usage: lat lon r sigma solver_file output_file update_param')
  end if

  ! read the location and the sigma of the Gaussian perturbation
  read(latt,*) lat
  read(lonn,*) lon
  read(rr,*) r
  read(ssigma,*) sigma
  if (myrank == 0) print*,"input parameters:", lat,lon, r, ssigma, trim(update_params)

! convert lat, lon, r and sigma to unit sphere

  phi = lon * DEGREES_TO_RADIANS;
  theta = (90.0 - lat) * DEGREES_TO_RADIANS
  r =(EARTH_R_KM - r) / EARTH_R_KM
  xloc = r * sin(theta) * cos(phi)
  yloc = r * sin(theta) * sin(phi)
  zloc = r * cos(theta)
  sigma = sigma / EARTH_R_KM

  if (myrank == 0) print*,r,xloc,yloc,zloc,sigma
  if (myrank == 0) print*,DEGREES_TO_RADIANS,EARTH_R_KM, NGLLX, NGLLY, NGLLZ, NSPEC

  write(solver_bp_file,'(a,i6.6,a)') trim(solver_file)//'/solver_data.bp'
  write(output_bp_file,'(a,i6.6,a)') trim(output_file)//'/model_pert.bp'
  if (myrank == 0) print*,"solver_file",trim(solver_bp_file)
  if (myrank == 0) print*,"output_file",trim(output_bp_file)

  

   kernel_name=(/character(len=150) :: &
       "reg1/dvph", &
       "reg1/dvpv", &
       "reg1/dvsh", &
       "reg1/dvsv", &
       "reg1/deta", &
       "reg1/drho", &
       "reg1/dqmu" /)

  solver_name=(/character(len=350) :: &
       "reg1/x_global", "reg1/y_global", "reg1/z_global"/)

  ! Open the mesh file
  call init_adios_group(myadios_group, "MeshReader")
  call open_file_adios_read_and_init_method(myadios_file, &
          myadios_group, trim(solver_bp_file))
  ! Read ibool
  ibool = 0.0
  call read_adios_array(myadios_file, myadios_group, &
       myrank, NSPEC, "reg1/ibool", ibool)
  !  if (myrank == 0) print*,ibool
  
  ! Read global positioning
  do iker=1,3
    solver_point = 0.0
    call read_adios_array(myadios_file, myadios_group, &
         myrank, NGLOB, trim(solver_name(iker)), solver_point)
    if (iker==1) then
       x=solver_point
    elseif (iker==2) then
       y=solver_point
    elseif (iker==3) then
       z=solver_point
    endif
    
  
  end do

  ! end of reading parameters from the mesh

   do iker=1,NKERNEL

      ! perturb bulk_betav_kl453368

      ! if (myrank == 0) print *, "kernel_name:", iker, kernel_name(iker)(6:9)
      if (trim(kernel_name(iker)(7:9)) == update_params) then
         do ispec = 1,NSPEC
            do k = 1,NGLLZ
               do j = 1,NGLLY
                  do i = 1,NGLLX
                     iglob=ibool(i,j,k,ispec)
                     !if (myrank == 0) print*,'iglob=',iglob,ispec,i,j,k

                     !dist=sqrt( (x(iglob)-xloc)**2 + (y(iglob)-yloc)**2 + (z(iglob)-zloc)**2 ) 
                     iglob = ibool(i,j,k,ispec)
                     !if (myrank == 45) print*, dist

                     exp_val = exp( - (x(iglob)-xloc)**2/sigma**2 - (y(iglob)-yloc)**2/sigma**2 - (z(iglob)-zloc)**2/sigma**2)
                     
                     new_kernel(i,j,k,ispec, iker) = exp_val
                     
                  end do
               end do
            end do
         end do
      else
         new_kernel(:,:,:,:,iker) = 0.0_CUSTOM_REAL
      end if

      call max_allreduce_cr(maxval(abs(new_kernel(:,:,:,:,iker))), max_new)
      if (myrank == 0) then
         print *, 'Max kernel (', trim(kernel_name(iker)), ") = ", max_new
      end if
   end do


  call init_adios_group(myadios_val_group, "Writer")
  group_size_inc = 0

  call open_file_adios_write(myadios_val_file, &
       myadios_val_group, trim(output_bp_file), "Writer")

  call set_adios_group_size(myadios_val_file, group_size_inc)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC
  do iker=1,NKERNEL
     !------------------------.
     ! Define ADIOS Variables |
     !------------------------'
     call define_adios_global_array1D(myadios_val_group, &
          group_size_inc, local_dim, '', &
          trim(kernel_name(iker)), new_kernel)
  enddo

  do iker=1,7
     !------------------------------------------.
     ! Write previously defined ADIOS variables |
     !------------------------------------------'
     local_dim = NGLLX * NGLLY * NGLLZ * NSPEC
     if (myrank == 0) print*, local_dim
     call write_adios_global_1d_array(myadios_val_file, &
          myadios_val_group, myrank, sizeprocs_adios, &
          local_dim, trim(kernel_name(iker)), &
          new_kernel(:, :, :, :, iker))
     call synchronize_all()
  end do
  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call write_adios_perform(myadios_val_file)

  call close_file_adios(myadios_val_file)
  call close_file_adios_read_and_finalize_method(myadios_file)
  if (myrank==0) write(*,*) 'done'

  call finalize_adios()
  call finalize_mpi()
  end program single_gaussian_perturb_for_psf
