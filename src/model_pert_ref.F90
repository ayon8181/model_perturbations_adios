program model_pert_ref

! ADIOS implementation: Ebru & Matthieu
! Princeton, September 2013
! Updated adios lib
! Ridvan Orsvuran, 2021

  use adios_helpers_mod
  use manager_adios

  implicit none

  include 'constants.h'
  include 'values_from_mesher.h'

#if defined(USE_AZIMUTHAL_ANISOTROPY)
  integer, parameter :: NPARAMS=7  ! vpv, vph, vsv, vsh, eta, rho, Gc_prime, Gs_prime
#elif defined(USE_QMU)
  integer, parameter :: NPARAMS=7  ! qmu
#else
  integer, parameter :: NPARAMS=7  ! vpv, vph, vsv, vsh, eta, rho
#endif

  integer,parameter :: NSPEC = NSPEC_CRUST_MANTLE
  !integer,parameter :: NSPEC_OC = NSPEC_OUTER_CORE
  !integer,parameter :: NSPEC_IC = NSPEC_INNER_CORE
  

  integer :: myrank, sizeprocs, iker
  character(len=512):: dir_model_new,dir_model_old,dir_model_pert
  character(len=256) :: model_name(NPARAMS)
  character(len=256) :: model_name_pert(NPARAMS)


  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NPARAMS) :: &
       model_new, model_old, model_pert


  ! adios
  integer(kind=8) :: group_size_inc, local_dim

  call init_mpi()

  call world_size(sizeprocs)
  call world_rank(myrank)

  call initialize_adios()
  call synchronize_all()

! read in input directory
  call getarg(1,dir_model_new)
  call getarg(2,dir_model_old)
  call getarg(3,dir_model_pert)
  
  if (myrank == 0) then
    print*,trim(dir_model_new)
    print*,trim(dir_model_old)
    print*,trim(dir_model_pert)
  end if
  if (trim(dir_model_new) == '' &
       .or. trim(dir_model_old) == '' &
       .or. trim(dir_model_pert) == '' ) then
     call exit_MPI(myrank, 'Usage: xmodel_pert_ref dir_model_new dir_model_old dir_model_pert')
  end if

  ! initial array
  model_new = 0.0_CUSTOM_REAL
  model_old = 0.0_CUSTOM_REAL
  model_pert = 0.0_CUSTOM_REAL

#if defined(USE_AZIMUTHAL_ANISOTROPY)
  model_name=(/character(len=150) :: &
       "reg1/vpv","reg1/vph", &
       "reg1/vsv","reg1/vsh", &
       "reg1/eta","reg1/rho", &
       "reg1/Gc_prime","reg1/Gs_prime"/)

  model_name_pert=(/character(len=150) :: &
       "vpv_crust_mantle","vph_crust_mantle", &
       "vsv_crust_mantle","vsh_crust_mantle", &
       "reg1/eta","reg1/rho", &
       "reg1/Gc_prime","reg1/Gs_prime"/)
#elif defined(USE_QMU)
  model_name=(/character(len=150) :: &
       "reg1/vpv","reg1/vph", &
       "reg1/vsv","reg1/vsh",&
       "reg1/eta","reg1/rho", &
       "reg1/qmu"/)
  model_name_pert=(/character(len=150) :: &
       "reg1/dvpv","reg1/dvph", &
       "reg1/dvsv","reg1/dvsh",&
       "reg1/deta","reg1/drho", &
       "reg1/dqmu"/)
#else
  model_name_pert=(/character(len=150) :: &
       "reg1/dvpv","reg1/dvph", &
       "reg1/dvsv","reg1/dvsh",&
       "reg1/deta","reg1/drho", &
       "reg1/dqmu"/)

  model_name=(/character(len=150) :: &
       "reg1/vpv","reg1/vph", &
       "reg1/vsv","reg1/vsh",&
       "reg1/eta","reg1/rho", &
       "reg1/qmu"/)
#endif
#if defined(USE_QMU)
  ! Read old model
  call init_adios_group(myadios_group, "ModelReader")
  call open_file_adios_read_and_init_method(myadios_file, &
       myadios_group, trim(dir_model_old)//'/model_gll_perturbed.bp')

  if (myrank == 0) print *, "old:", trim(dir_model_old)

  do iker=1,NPARAMS
     call read_adios_array(myadios_file, myadios_group, &
          myrank, NSPEC, trim(model_name(iker)), &
          model_old(:, :, :, :, iker))
  enddo
  call close_file_adios_read_and_finalize_method(myadios_file)

#else
  call init_adios_group(myadios_group, "ModelReader")
  call open_file_adios_read_and_init_method(myadios_file, &
       myadios_group, trim(dir_model_old)//'/model_gll.bp')

  if (myrank == 0) print *, "new:", trim(dir_model_old)

  do iker=1,NPARAMS
     call read_adios_array(myadios_file, myadios_group, &
          myrank, NSPEC, trim(model_name(iker)), &
          model_old(:, :, :, :, iker))
  enddo
  call close_file_adios_read_and_finalize_method(myadios_file)
#endif


  ! Read perturbation model
  call open_file_adios_read_and_init_method(myadios_file, &
       myadios_group, trim(dir_model_pert)//'/model_pert.bp')

  if (myrank == 0) print *, "old:", trim(dir_model_pert)

  do iker=1,NPARAMS
     call read_adios_array(myadios_file, myadios_group, &
          myrank, NSPEC, trim(model_name_pert(iker)), &
          model_pert(:, :, :, :, iker))
  enddo
  call close_file_adios_read_and_finalize_method(myadios_file)

  ! Calculate Perturbed Model
  do iker=1,NPARAMS
#if defined(USE_QMU)
     ! model: Qmu , model_pert: Qmu_new^{-1} - Qmu_old^{-1}
     if (iker == 7) then
         model_new(:,:,:,:,iker) = 1.0/ ((1.0/model_old(:,:,:,:,iker)) + 1.0/(model_old(:,:,:,:,iker)))
         if (myrank == 0) print*, "Updated", model_name(iker)
     else 
          model_new(:,:,:,:,iker) = model_old(:,:,:,:,iker)
     endif
#else
     if (iker==7) then
        model_new(:,:,:,:,iker) = model_old(:,:,:,:,iker)
        
     else
        model_new(:,:,:,:,iker) = model_old(:,:,:,:,iker) * exp(model_pert(:,:,:,:,iker))
        if (myrank == 0) print*, "Updated", model_name(iker)
     endif
#endif
  enddo

  ! Open output (perturbation) file
  call init_adios_group(myadios_val_group, "Writer")
#if defined(USE_AZIMUTHAL_ANISOTROPY)
  call open_file_adios_write(myadios_val_file, myadios_val_group, &
       trim(dir_model_new) // '/model_gll_perturbed_azi.bp', "Writer")
#elif defined(USE_QMU)
  call open_file_adios_write(myadios_val_file, myadios_val_group, &
       trim(dir_model_new) // '/model_gll_perturbed_qmu.bp', "Writer")
#else
  call open_file_adios_write(myadios_val_file, myadios_val_group, &
       trim(dir_model_new) // '/model_gll_perturbed.bp', "Writer")
#endif

  ! define adios groups
  group_size_inc = 0
  call define_adios_scalar(myadios_val_group, group_size_inc, '', &
       "NSPEC", NSPEC)
  call define_adios_scalar(myadios_val_group, group_size_inc, '', &
       "reg1/nspec", NSPEC)

  do iker = 1,NPARAMS
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC
    if (myrank == 0) print*, local_dim
    call define_adios_global_array1D(myadios_val_group, &
         group_size_inc, local_dim, '', &
         trim(model_name(iker)), model_new(:,:,:,:,iker))
  enddo

  call set_adios_group_size(myadios_val_file, group_size_inc)

  ! write values
  call write_adios_scalar(myadios_val_file, myadios_val_group, &
       "NSPEC", NSPEC)
  call write_adios_scalar(myadios_val_file,myadios_val_group, &
       "reg1/nspec", NSPEC)

  do iker=1,NPARAMS
     if (myrank == 0) print*,iker
     call write_adios_array_gll(myadios_val_file, myadios_val_group, &
          myrank,sizeprocs_adios,NSPEC, &
          trim(model_name(iker)), model_new(:,:,:,:,iker))
  enddo
  call write_adios_perform(myadios_val_file)

  call close_file_adios(myadios_val_file)
  call finalize_adios()
  call finalize_mpi()
end program model_pert_ref
