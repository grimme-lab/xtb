module xtb_solv_cpx
    use cpx, only: calculation_type, parameter_type, atomicmass, density,&
    &initialize_param,load_solvent,read_cosmo
    use xtb_solv_cosmo, only: TCosmo
    use xtb_type_environment, only: TEnvironment
    use iso_fortran_env, only: input_unit, output_unit
    use mctc_env, only: error_type, fatal_error, wp
    implicit none
    private
    public :: TCpcmx

    character(len=*), parameter :: source = 'xtb_solv_cpx'

    !> CPCM-X calculation type
    type, extends(calculation_type) :: TCpcmx
        contains

        procedure :: setup => setup_cpcmx
        procedure :: calc_solv => calculate_cpcmx
        procedure :: print => print_cpcmx
    end type Tcpcmx

    contains

    subroutine setup_cpcmx(self, env, solvent)
        implicit none
        class(Tcpcmx), intent(inout) :: self
        type(TEnvironment), intent(inout) :: env
        !> Solvent name for internal CPCM-X database
        character(len=*), intent(in) :: solvent
        !> Error handling
        type(error_type), allocatable :: error

        !> Load xTB parameters from internal CPCM-X database
        call initialize_param('xtb',solvent,self,error)
        if (allocated(error)) Call env%error(error%message, source)
        call load_solvent(solvent,self%solvent,error)
        if (allocated(error)) Call env%error(error%message, source)
        !> TODO: SETUP solute directly from TCosmo instead of reading COSMO file
        Call read_cosmo('xtb.cosmo',self%solute,'NONE',error)
        if (allocated(error)) Call env%error(error%message, source)
    end subroutine setup_cpcmx

    subroutine calculate_cpcmx(self,env,solvent,probe,T,max_cycle,conv_crit)
        implicit none
        class(Tcpcmx), intent(inout) :: self
        type(TEnvironment), intent(inout) :: env
        !> Solvent name for solvent properties (density, atomic mass) from internal CPCM-X database
        character(len=*), intent(in) :: solvent
        !> Error handling
        type(error_type), allocatable :: error
        !> Probe radius
        real(wp), intent(in) :: probe
        !> Temperature
        real(wp), intent(in) :: T
        !> Convergence criteria for restoring free energy
        integer, intent(in) :: max_cycle
        real(wp), intent(in) :: conv_crit

        !> Charge averaging
        call self%average_charge(error)
        if (allocated(error)) Call env%error(error%message, source)
        call self%init_bonding()
        !> Calculate restoring free energy
        call self%solv('crs',error,T,max_cycle,conv_crit)
        if (allocated(error)) Call env%error(error%message, source)
        !> Calculate SMD contribution and state correction
        call self%state_correction(density(solvent),atomicmass(self%solvent%element),T)
        call self%cds(probe,solvent)
    end subroutine calculate_cpcmx

    !> Print CPCM-X results
    subroutine print_cpcmx(self,verbose)
        implicit none
        class(Tcpcmx), intent(in) :: self
        !> Output level
        logical, intent(in), optional :: verbose

        real(wp), parameter :: autokcal=627.509_wp
        logical :: pr

        if (present(verbose)) then
            pr=verbose
        else
            pr=.false.
        end if

        write(output_unit,*) ""
        write(output_unit,'(5x,a,t55,a,t66,a)') &
            "Free Energy contributions:", "[Eh]", " [kcal/mol]"
        write(output_unit,'(4x,a)') repeat('-',73)
        if (pr) then
            write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
                "Ideal State (dG_IS):", self%dG_is, self%dG_is*autokcal, &
                "Averaging correction (dG_av):", self%dG_cc, self%dG_cc*autokcal, &
                "restoring free energy (dG_res):", self%dG_res, self%dG_res*autokcal, &
                "SMD Contribution (dG_CDS):", self%dG_smd, self%dG_smd*autokcal, &
                "Standard state correction (dG_corr):", self%dG_ss, self%dG_ss*autokcal, &
                "Systematic empirical shift (dG_shift)", self%dG_shift, self%dG_shift*autokcal
            write(output_unit,'(4x,a)') repeat('-',73)
        else
        end if
        write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
        "solvation free energy (dG_solv): ", self%dG(), self%dG()*autokcal, &
        "gas phase energy (E)", self%solute%energy_gas
        write(output_unit,'(4x,a)') repeat('-',73)
        write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
        "total free energy (dG)", self%dG()+self%solute%energy_gas
        write(output_unit,*) ""

    end subroutine print_cpcmx

end module xtb_solv_cpx