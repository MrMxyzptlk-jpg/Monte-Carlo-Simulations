program electron_scatter
    use precision
    use formats
    use materials
    use mzranmod
    use electron_module
    implicit none

    integer                     :: Num_electrons, Num_trajectories
    integer                     :: i, nBSE, step
    type(Electron), allocatable :: sample(:)
    real(pr)                    :: R1, R2, R3, covering_thickness, ray_tilt, BSE_tol
    character(len=30)           :: traj_file, electron_id
    character(len=15)           :: particle_type, bulk_material, surface_material
    integer                     :: traj_unit, unitnum
    logical                     :: covering, track_trajectories
    real(pr)                    :: E_min, E_max


!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /physical/ bulk_material, covering, surface_material, ray_tilt, covering_thickness, E_min, E_max
    namelist /calculation/ Num_electrons, Num_trajectories, BSE_tol
    namelist /tasks/ particle_type, track_trajectories

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
        bulk_material    = "Au"
        covering         = .false.
        surface_material = "Au"
        ray_tilt         = 45._pr
        covering_thickness = 2e-6_pr
        E_min              = 0.5_pr ! keV
        E_max              = 30._pr ! keV

        !Calculation settings
        Num_electrons = 250
        Num_trajectories = 10
        BSE_tol = 1e-6

        ! Tasks
        particle_type = "electrons"
        track_trajectories = .false.

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
        read(unitnum, nml=tasks)
    close(unitnum)

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    call init_materials(bulk_material, surface_material, covering, covering_thickness)
    call init_angles(ray_tilt)
    call init_energy_range(E_max, E_min)
    call mzran_init()

    allocate(sample(Num_electrons))
    nBSE = 0
    select case(track_trajectories)
        case(.true.)
            do i = 1, Num_trajectories
                ! Open trajectory file
                write(electron_id, '(I6.6".dat")') i
                electron_id = adjustl(trim(electron_id))
                traj_file = "datos/trajectory_id_"//electron_id
                open(newunit=traj_unit, file=traj_file, status="replace", action="write")
                write(traj_unit, format_style_header) "## Step", "|", "x","|","y","|", "z","|","Energy"

                ! Initialize particle
                call sample(i)%init(0.0_pr, 0.0_pr, 0.0_pr)
                call MC_step(sample, i)
                step = 0

                ! Walk randomly
                do while (sample(i)%energy > E_min)
                    call MC_step(sample, i)

                    if (sample(i)%z < -BSE_tol ) then
                        nBSE = nBSE + 1
                        exit
                    end if

                    step = step + 1
                    call track_trajectory(traj_unit, step, sample(i)%x, sample(i)%y, sample(i)%z, sample(i)%energy)
                end do
            end do
            do i = Num_trajectories + 1, Num_electrons

                ! Initialize particle
                call sample(i)%init(0.0_pr, 0.0_pr, 0.0_pr)
                call MC_step(sample, i)
                step = 0

                ! Walk randomly
                do while (sample(i)%energy > E_min)
                    call MC_step(sample, i)
                    if (sample(i)%z < -BSE_tol ) then
                        nBSE = nBSE + 1
                        exit
                    end if
                end do

            end do

        case(.false.)
            do i = 1, Num_electrons

                ! Initialize particle
                call sample(i)%init(0.0_pr, 0.0_pr, 0.0_pr)
                call MC_step(sample, i)
                step = 0

                ! Walk randomly
                do while (sample(i)%energy > E_min)
                    call MC_step(sample, i)
                    if (sample(i)%z < -BSE_tol ) then
                        nBSE = nBSE + 1
                        exit
                    end if
                end do

            end do
    end select

    close(traj_unit)

    print *, "Backscattered electrons:", real(nBSE)/real(Num_electrons)

    deallocate(sample)

end program electron_scatter
