! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Two-step SCF of the PTB method

module xtb_ptb_scf
   use mctc_io, only: structure_type
   use mctc_env, only: wp, error_type, fatal_error

   use tblite_basis_type, only: basis_type, get_cutoff
   use tblite_blas, only: gemv
   use tblite_context, only: context_type
   use tblite_scf_solver, only: solver_type
   use tblite_scf_iterator, only: get_qat_from_qsh
   use tblite_adjlist, only: adjacency_list, new_adjacency_list
   use tblite_cutoff, only: get_lattice_points
   use tblite_wavefunction, only: wavefunction_type, get_alpha_beta_occupation
   use tblite_wavefunction_fermi, only: get_fermi_filling
   use tblite_wavefunction_type, only: get_density_matrix
   use tblite_wavefunction_mulliken, only: get_mulliken_atomic_multipoles, &
      & get_mulliken_shell_charges, get_mayer_bond_orders, get_molecular_quadrupole_moment, &
      & get_molecular_dipole_moment
   use tblite_scf_potential, only: potential_type, new_potential, add_pot_to_h1
   use tblite_integral_type, only: new_integral, integral_type
   use tblite_container, only: container_type, container_cache
   use tblite_external_field, only: electric_field

   use multicharge_model, only: mchrg_model_type

   use xtb_ptb_vdzp, only: add_vDZP_basis, nshell, max_shell
   use xtb_ptb_param, only: ptbGlobals, rf
   use xtb_ptb_integrals, only: get_integrals
   use xtb_ptb_mmlpopanalysis, only: get_mml_overlaps, get_mml_shell_charges
   use xtb_ptb_ncoord, only: ncoord_erf
   use xtb_ptb_corebasis, only: get_Vecp
   use xtb_ptb_data, only: TPTBData
   use xtb_ptb_hamiltonian, only: get_hamiltonian, get_selfenergy, get_occupation
   use xtb_ptb_guess, only: guess_qsh, get_psh_from_qat, get_psh_from_qsh
   use xtb_ptb_paulixc, only: calc_Vxc_pauli
   use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral
   use xtb_ptb_coulomb, only: coulomb_potential
   use xtb_ptb_plusu, only: plusu_potential_type

   use xtb_readin, only: bool2string, bool2int

   implicit none
   private

   public :: twostepscf, get_density

   character(len=*), private, parameter :: outfmt = &
                                           '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'
   character(len=*), parameter :: intfmt = &
                                  '(10x,":",2x,a,i18,      10x,":")'
   character(len=*), parameter :: chrfmt = &
                                  '(10x,":",2x,a,a18,      10x,":")'
   character(len=*), parameter :: scifmt = &
                                  '(10x,":",2x,a,e22.7,1x,a,1x,":")'
   character(len=*), parameter :: dblfmt = &
                                  '(10x,":",2x,a,f18.7,5x,a,1x,":")'
   character(len=*), parameter :: source = "twostepscf"

   real(wp), parameter :: default_cutoff = 25.0_wp
   !> Conversion factor from temperature to energy
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains

   !> Final Hamiltonian to solve (2nd iteration) on ints%hamiltonian
   subroutine twostepscf(ctx, wfn, data, mol, bas, cbas, ints, auxints, &
         & eeqmodel, dipole, vecp, list, levels, v_es_sh, cn_star, wbo, efield)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(inout) :: wfn
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> PTB parameterization data
      type(TPTBData), intent(in) :: data
      !> Basis set and core-valence basis set data
      type(basis_type), intent(in) :: bas, cbas
      !> Integral type
      type(integral_type), intent(out) :: ints
      !> Auxiliary integral type
      type(aux_integral_type), intent(out) :: auxints
      !> Initialized EEQ model
      type(mchrg_model_type), intent(in) :: eeqmodel
      !> Molecular dipole moment
      real(wp), intent(out) :: dipole(3)
      real(wp) :: quadrupole(6)
      !> Wiberg bond orders
      real(wp), allocatable, intent(out) :: wbo(:, :, :)
      !> Effective core potential
      real(wp), allocatable, intent(out) :: vecp(:, :)
      !> Adjacency list
      type(adjacency_list), intent(out) :: list
      !> Effective self-energies
      real(wp), allocatable, intent(out) :: levels(:)
      !> Electrostatic potential in second iteration
      real(wp), allocatable, intent(out) :: v_es_sh(:)
      !> (optional) Electric field
      real(wp), intent(in), optional :: efield(:)
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Error type
      type(error_type), allocatable :: error
      !> Coulomb potential
      type(coulomb_potential) :: coulomb
      !> +U potential
      type(plusu_potential_type) :: plusu
      !> Potential type
      type(potential_type) :: pot
      !> H0 basis in second iteration
      type(basis_type), allocatable :: bas_h0
      !> Restart data for interaction containers
      type(container_cache) :: icache
      !> Electric field object
      type(electric_field) :: efield_object
      real(wp), allocatable :: expscal_h0_2nditer(:, :)
      !> Loop variables
      integer :: i, j, isp, izp, iat, ish, iid, is
      !> Coordination numbers
      real(wp), intent(out) :: cn_star(mol%nat)
      real(wp) :: cn(mol%nat), cn_eeq(mol%nat)
      real(wp) :: radii(mol%nid)
      !> Lattice points
      real(wp), allocatable :: lattr(:, :)
      !> Cutoff for lattice points
      real(wp) :: cutoff
      !> Number of electrons
      real(wp) :: nel
      !> Pauli XC potential
      real(wp), allocatable :: psh(:, :)
      !> Electronic entropy
      real(wp) :: ts
      !> Tmp variable for dipole moment
      real(wp) :: tmpdip(3)
      real(wp), allocatable :: mulliken_qsh(:, :), mulliken_qat(:, :)

      !> Solver for the effective Hamiltonian
      call ctx%new_solver(solver, bas%nao)

      if (present(efield)) then
         efield_object = electric_field(efield)
      end if

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)

      !            _____                    _____                _____                    _____                    _____
      !           /\    \                  /\    \              /\    \                  /\    \                  /\    \
      !          /::\    \                /::\    \            /::\    \                /::\____\                /::\    \
      !         /::::\    \              /::::\    \           \:::\    \              /:::/    /               /::::\    \
      !        /::::::\    \            /::::::\    \           \:::\    \            /:::/    /               /::::::\    \
      !       /:::/\:::\    \          /:::/\:::\    \           \:::\    \          /:::/    /               /:::/\:::\    \
      !      /:::/__\:::\    \        /:::/__\:::\    \           \:::\    \        /:::/    /               /:::/__\:::\    \
      !      \:::\   \:::\    \      /::::\   \:::\    \          /::::\    \      /:::/    /               /::::\   \:::\    \
      !    ___\:::\   \:::\    \    /::::::\   \:::\    \        /::::::\    \    /:::/    /      _____    /::::::\   \:::\    \
      !   /\   \:::\   \:::\    \  /:::/\:::\   \:::\    \      /:::/\:::\    \  /:::/____/      /\    \  /:::/\:::\   \:::\____\
      !  /::\   \:::\   \:::\____\/:::/__\:::\   \:::\____\    /:::/  \:::\____\|:::|    /      /::\____\/:::/  \:::\   \:::|    |
      !  \:::\   \:::\   \::/    /\:::\   \:::\   \::/    /   /:::/    \::/    /|:::|____\     /:::/    /\::/    \:::\  /:::|____|
      !   \:::\   \:::\   \/____/  \:::\   \:::\   \/____/   /:::/    / \/____/  \:::\    \   /:::/    /  \/_____/\:::\/:::/    /
      !    \:::\   \:::\    \       \:::\   \:::\    \      /:::/    /            \:::\    \ /:::/    /            \::::::/    /
      !     \:::\   \:::\____\       \:::\   \:::\____\    /:::/    /              \:::\    /:::/    /              \::::/    /
      !      \:::\  /:::/    /        \:::\   \::/    /    \::/    /                \:::\__/:::/    /                \::/____/
      !       \:::\/:::/    /          \:::\   \/____/      \/____/                  \::::::::/    /                  ~~
      !        \::::::/    /            \:::\    \                                    \::::::/    /
      !         \::::/    /              \:::\____\                                    \::::/    /
      !          \::/    /                \::/    /                                     \::/____/
      !           \/____/                  \/____/                                       ~~

      !    _____           _                                   _
      !   |_   _|         | |                                 | |
      !     | |    _ __   | |_    ___    __ _   _ __    __ _  | |  ___
      !     | |   | '_ \  | __|  / _ \  / _` | | '__|  / _` | | | / __|
      !    _| |_  | | | | | |_  |  __/ | (_| | | |    | (_| | | | \__ \
      !   |_____| |_| |_|  \__|  \___|  \__, | |_|     \__,_| |_| |___/
      !                                  __/ |
      !                                 |___/

      call new_integral(ints, bas%nao)
      call new_aux_integral(auxints, bas%nao)
      call get_integrals(mol, bas, lattr, list, ints%overlap, ints%dipole, ints%quadrupole, norm=auxints%norm)
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f8.5)', advance="no") ints%overlap(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      ! write (*, *) "Dipole moment integrals:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(2i3,3f8.4)') i, j, ints%dipole(:, i, j)
      !    end do
      ! end do
      ! write (*, *) "Quadrupole moment integrals:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(2i3,6f8.4)') i, j, ints%quadrupole(:, i, j)
      !    end do
      ! end do
      !#####################
      call get_integrals(mol, lattr, list, auxints%overlap_h0_1, alpha_scal=id_to_atom(mol, data%hamiltonian%kalphah0l))
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      call get_integrals(mol, lattr, list, auxints%overlap_xc, alpha_scal=id_to_atom(mol, data%pauli%klalphaxc))
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_xc(i, j)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################
      call get_mml_overlaps(bas, ints%overlap, ptbGlobals%mlmix, auxints%overlap_to_x, &
      & auxints%overlap_to_1_x)

      !     _____   _   _
      !    / ____| | \ | |
      !   | |      |  \| |
      !   | |      | . ` |
      !   | |____  | |\  |
      !    \_____| |_| \_|

      !> Get first coordination number (" CN' ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn_star)
      !##### DEV WRITE #####
      ! write (*, *) "CN star:"
      ! do i = 1, mol%nat
      !    write (*, *) "Atom ", i, ":", cn_star(i)
      ! end do
      !#####################
      !> Get radii from PTB parameters for second coordination number (" CN ")
      do isp = 1, mol%nid
         izp = mol%num(isp)
         radii(isp) = rf(izp)
      end do
      !> Get second coordination number (" CN ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn, covrad=radii)
      !##### DEV WRITE #####
      ! write (*, *) "CN:"
      ! do i = 1, mol%nat
      !    write (*, *) "Atom ", i, ":", cn(i)
      ! end do
      !#####################
      !> Get third coordination number for EEQ model (" CN-EEQ ")
      call ncoord_erf(mol, ptbGlobals%kerfcn_eeq, default_cutoff, cn_eeq)

      !##### DEV WRITE #####
      ! write (*, *) "CN-EEQ:"
      ! do i = 1, mol%nat
      !    write (*, *) "Atom ", i, ":", cn_eeq(i)
      ! end do
      !#####################

      !    _    _    ___
      !   | |  | |  / _ \
      !   | |__| | | | | |
      !   |  __  | | | | |
      !   | |  | | | |_| |
      !   |_|  |_|  \___/

      !> V_ECP via PTB core basis
      call get_Vecp(mol, data%corepotential, bas, cbas, auxints%norm, vecp)
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") vecp(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      !> Get the effective self-energies
      call get_selfenergy(mol, bas, data%hamiltonian, cn, cn_star, levels)

      !    ______   _                 _
      !   |  ____| | |               | |
      !   | |__    | |   ___    ___  | |_   _ __    ___
      !   |  __|   | |  / _ \  / __| | __| | '__|  / _ \
      !   | |____  | | |  __/ | (__  | |_  | |    | (_) |
      !   |______| |_|  \___|  \___|  \__| |_|     \___/

      !> EEQ call
      call eeqmodel%solve(mol, cn_eeq, qvec=wfn%qat(:, 1))
      !##### DEV WRITE #####
      ! write (*, *) "EEQ charges:"
      ! do i = 1, mol%nat
      !    write (*, '(a,i0,a,f12.6)') "Atom ", i, ":", wfn%qat(i, 1)
      ! end do
      !#####################
      !> Project reference occupation on wavefunction and use EEQ charges as guess
      call get_occupation(mol, bas, data%hamiltonian%refocc, wfn%nocc, wfn%n0at, wfn%n0sh)
      call guess_qsh(wfn, bas)
      !##### DEV WRITE #####
      ! do i = 1, bas%nsh
      !    write(*,*) wfn%qsh(i, 1)
      ! enddo
      !#####################
      call coulomb%init(mol, bas, wfn%qat(:, 1), data%coulomb%shellHardnessFirstIter, &
         & data%coulomb%kQHubbard, data%coulomb%kOK1, data%coulomb%kTO)
      call coulomb%update(mol, bas)

      if (ctx%verbosity > 1) then
         write (ctx%unit, '(/,10x,51("."))')
         write (ctx%unit, '(10x,":",22x,a,22x,":")') "SETUP"
         write (ctx%unit, '(10x,":",49("."),":")')
         write (ctx%unit, intfmt) "# atomic orbitals  ", bas%nao
         write (ctx%unit, intfmt) "# shells           ", bas%nsh
         write (ctx%unit, intfmt) "# electrons        ", nint(wfn%nocc)
         write (ctx%unit, intfmt) "# open shells      ", mol%uhf
         write (ctx%unit, intfmt) "max. iterations    ", 2
         write (ctx%unit, chrfmt) "Hamiltonian        ", "PTB"
         write (ctx%unit, chrfmt) "PC potential       ", bool2string(.false.)
         ! if (lpcem) then
         !    write (env%unit, intfmt) "-> # point charges ", pcem%n
         !    write (env%unit, dblfmt) "-> sum of PC       ", sum(pcem%q), "e   "
         ! end if
         write (ctx%unit, dblfmt) "electronic temp.   ", wfn%kt / kt, "K   "
         ! write (env%unit, dblfmt) "accuracy           ", acc, "    "
         write (ctx%unit, scifmt) "-> integral cutoff ", bas%intcut, "    "
         ! write (env%unit, scifmt) "-> integral neglect", neglect, "    "
         write (ctx%unit, intfmt) "verbosity level    ", ctx%verbosity
         write (ctx%unit, '(10x,51("."))')
      end if

      if (ctx%verbosity > 1) then
         write (ctx%unit, '(/,10x,a)') "--- Calculation progress ---"
         write (ctx%unit, '(14x,a)') "1st iteration..."
      end if

      !           _____                    _____                    _____
      !         /\    \                  /\    \                  /\    \
      !        /::\    \                /::\    \                /::\    \
      !       /::::\    \              /::::\    \              /::::\    \
      !      /::::::\    \            /::::::\    \            /::::::\    \
      !     /:::/\:::\    \          /:::/\:::\    \          /:::/\:::\    \
      !    /:::/__\:::\    \        /:::/  \:::\    \        /:::/__\:::\    \
      !    \:::\   \:::\    \      /:::/    \:::\    \      /::::\   \:::\    \
      !  ___\:::\   \:::\    \    /:::/    / \:::\    \    /::::::\   \:::\    \
      ! /\   \:::\   \:::\    \  /:::/    /   \:::\    \  /:::/\:::\   \:::\    \
      !/::\   \:::\   \:::\____\/:::/____/     \:::\____\/:::/  \:::\   \:::\____\
      !\:::\   \:::\   \::/    /\:::\    \      \::/    /\::/    \:::\   \::/    /
      ! \:::\   \:::\   \/____/  \:::\    \      \/____/  \/____/ \:::\   \/____/
      !  \:::\   \:::\    \       \:::\    \                       \:::\    \
      !   \:::\   \:::\____\       \:::\    \                       \:::\____\
      !    \:::\  /:::/    /        \:::\    \                       \::/    /
      !     \:::\/:::/    /          \:::\    \                       \/____/
      !      \::::::/    /            \:::\    \
      !       \::::/    /              \:::\____\
      !        \::/    /                \::/    /
      !         \/____/                  \/____/

      !   _     _     _ _                 _   _
      !  / |___| |_  (_) |_ ___ _ __ __ _| |_(_) ___  _ __
      !  | / __| __| | | __/ _ \ '__/ _` | __| |/ _ \| '_ \
      !  | \__ \ |_  | | ||  __/ | | (_| | |_| | (_) | | | |
      !  |_|___/\__| |_|\__\___|_|  \__,_|\__|_|\___/|_| |_|

      !> Set up the effective Hamiltonian in the first iteration
      call new_potential(pot, mol, bas, wfn%nspin)
      !    _    _  ___
      !   | |  | |/ _ \
      !   | |__| | | | |
      !   |  __  | | | |
      !   | |  | | |_| |
      !   |_|  |_|\___/
      !
      !>  --------- Get H0 (wavefunction-independent (but iteration-dependent)) ----------

      ints%hamiltonian = 0.0_wp
      call get_hamiltonian(mol, list, bas, data%hamiltonian, data%hamiltonian%kla, auxints%overlap_h0_1, &
      & levels, ints%hamiltonian, ptbGlobals%kpol, ptbGlobals%kitr, ptbGlobals%kitocod)
      ints%hamiltonian = ints%hamiltonian + vecp
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f8.4)', advance="no") ints%hamiltonian(i, j)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      !    _____      _             _   _       _
      !   |  __ \    | |           | | (_)     | |
      !   | |__) |__ | |_ ___ _ __ | |_ _  __ _| |
      !   |  ___/ _ \| __/ _ \ '_ \| __| |/ _` | |
      !   | |  | (_) | ||  __/ | | | |_| | (_| | |
      !   |_|   \___/ \__\___|_| |_|\__|_|\__,_|_|
      !
      !>  --------- Get potential (wavefunction-dependent) -----

      call pot%reset()
      !> Coulomb potential
      call coulomb%get_potential(wfn, pot)
      call add_pot_to_h1(bas, ints, pot, wfn%coeff)

      !> Pauli XC potential
      allocate (psh(bas%nsh, wfn%nspin), source=0.0_wp)
      psh = get_psh_from_qat(wfn, bas)
      call calc_Vxc_pauli(mol, bas, psh(:, 1), auxints%overlap_xc, levels, data%pauli%kxc1, wfn%coeff(:, :, 1))

      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f11.7)', advance="no") wfn%coeff(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      !   __       _
      !  / _\ ___ | |_   _____
      !  \ \ / _ \| \ \ / / _ \
      !  _\ \ (_) | |\ V /  __/
      !  \__/\___/|_| \_/ \___|

      !> Project occupations on alpha and beta orbitals
      nel = sum(wfn%n0at) - mol%charge
      if (mod(mol%uhf, 2) == mod(nint(nel), 2)) then
         wfn%nuhf = mol%uhf
      else
         call fatal_error(error, "Number of electrons and spin multiplicity do not match.")
      end if
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2))
      call get_density(wfn, solver, ints, ts, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if

      call get_mml_shell_charges(bas, auxints%overlap_to_x, auxints%overlap_to_1_x, &
         & wfn%density, wfn%n0sh, wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      psh = get_psh_from_qsh(wfn, bas)
      !##### DEV WRITE #####
      ! write (*, *) "Shell charges after 1st iteration ..."
      ! do i = 1, bas%nsh
      !    write (*, '(f8.4)') wfn%qsh(i, 1)
      ! end do
      ! write (*, *) "Atomic charges after 1st iteration ..."
      ! do i = 1, mol%nat
      !    write (*, '(f8.4)') wfn%qat(i, 1)
      ! end do
      !#####################

      if (ctx%verbosity > 1) then
         write (ctx%unit, '(14x,a)') "2nd iteration..."
      end if

      !   ____            _   _ _                 _   _
      !  |___ \ _ __   __| | (_) |_ ___ _ __ __ _| |_(_) ___  _ __
      !    __) | '_ \ / _` | | | __/ _ \ '__/ _` | __| |/ _ \| '_ \
      !   / __/| | | | (_| | | | ||  __/ | | (_| | |_| | (_) | | | |
      !  |_____|_| |_|\__,_| |_|\__\___|_|  \__,_|\__|_|\___/|_| |_|

      !    _    _  ___
      !   | |  | |/ _ \
      !   | |__| | | | |
      !   |  __  | | | |
      !   | |  | | |_| |
      !   |_|  |_|\___/

      !>  --------- Get H0 (wavefunction-independent (but iteration-dependent)) ----------
      !> Allocate temporary basis set and charge-dependent scaling factors
      !> Consequently, the basis set is not equal for same atom ids but different
      !> for each symmetry-unique atom
      allocate (bas_h0)
      allocate (expscal_h0_2nditer(max_shell, mol%nat), source=1.0_wp)
      do iat = 1, mol%nat
         iid = mol%id(iat)
         is = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iat)
            expscal_h0_2nditer(ish, iat) = data%hamiltonian%kalphah0l(ish, iid) * &
               & (1.0_wp + data%hamiltonian%kits0(iid) * wfn%qsh(is + ish, 1))
            !##### DEV WRITE #####
            ! write(*,*) "qsh: ", wfn%qsh(is + ish, 1)
            ! write (*, *) "Atom ", iat, ":", expscal_h0_2nditer(ish, iat)
            !#####################
         end do
      end do
      !> Get temporary basis set with charge-dependent scaled exponents
      call add_vDZP_basis(mol, expscal_h0_2nditer, bas_h0)
      !> Get integrals with temporary basis set
      call get_integrals(mol, bas_h0, lattr, list, auxints%overlap_h0_2)
      !> Deallocate temporary basis set and charge-dependent scaling factors
      deallocate (bas_h0, expscal_h0_2nditer)
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      ints%hamiltonian = 0.0_wp
      call get_hamiltonian(mol, list, bas, data%hamiltonian, data%hamiltonian%kla, auxints%overlap_h0_2, &
      & levels, ints%hamiltonian, ptbGlobals%kpol)
      ints%hamiltonian = ints%hamiltonian + vecp
      !##### DEV WRITE #####
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f11.7)', advance="no") ints%hamiltonian(i, j)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      !    _____      _             _   _       _
      !   |  __ \    | |           | | (_)     | |
      !   | |__) |__ | |_ ___ _ __ | |_ _  __ _| |
      !   |  ___/ _ \| __/ _ \ '_ \| __| |/ _` | |
      !   | |  | (_) | ||  __/ | | | |_| | (_| | |
      !   |_|   \___/ \__\___|_| |_|\__|_|\__,_|_|
      !
      !>  --------- Get potential (wavefunction-dependent) -----
      call coulomb%init(mol, bas, wfn%qat(:, 1), data%coulomb%shellHardnessSecondIter, &
         & 0.0_wp, data%coulomb%kOK2, data%coulomb%kTO)
      call plusu%init(data%plusU, mol, bas, wfn%qat(:, 1), cn_star)

      call pot%reset()
      call coulomb%update(mol, bas)
      call coulomb%get_potential(wfn, pot)
      allocate (v_es_sh(bas%nsh), source=0.0_wp)
      v_es_sh(:) = pot%vsh(:, 1)
      if (present(efield)) then
         call efield_object%update(mol, icache)
         call efield_object%get_potential(mol, icache, wfn, pot)
      end if
      call add_pot_to_h1(bas, ints, pot, wfn%coeff)
      call calc_Vxc_pauli(mol, bas, psh(:, 1), auxints%overlap_xc, levels, data%pauli%kxc2l, wfn%coeff(:, :, 1))
      call plusu%get_potential(mol, bas, wfn%density(:, :, 1), wfn%coeff(:, :, 1))

      ints%hamiltonian = wfn%coeff(:, :, 1)
      call get_density(wfn, solver, ints, ts, error, ptbGlobals%geps, ptbGlobals%geps0)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if

      call get_mml_shell_charges(bas, auxints%overlap_to_x, auxints%overlap_to_1_x, &
         & wfn%density, wfn%n0sh, wfn%qsh)
      psh = get_psh_from_qsh(wfn, bas)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
         & wfn%dpat)
      call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
         & wfn%qpat)
      !> This step is only required because the dipole integrals in tblite are not
      !> centered to a fixed point (e.g., origin) but are atomic dipole moments
      !> Calculation of dipole moments requires Mulliken atomic charges.
      allocate (mulliken_qsh(bas%nsh, wfn%nspin), mulliken_qat(mol%nat, wfn%nspin))
      call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & mulliken_qsh)
      call get_qat_from_qsh(bas, mulliken_qsh, mulliken_qat)
      !> Get the molecular dipole moment
      call gemv(mol%xyz, mulliken_qat(:, 1), tmpdip)
      dipole(:) = tmpdip + sum(wfn%dpat(:, :, 1), 2)

      !##### DEV TEST #####
      !> Alternative way to get the molecular dipole moment
      ! call get_molecular_dipole_moment(mol, mulliken_qat(:, 1), wfn%dpat(:, :, 1), &
      !      & dipole)
      !#####################
      call get_molecular_quadrupole_moment(mol, mulliken_qat(:, 1), wfn%dpat(:, :, 1), &
         & wfn%qpat(:, :, 1), quadrupole)

      write (*, *) "Quadrupole moments:"
      do i = 1, 6
         write (*, '(f10.6)') quadrupole(i)
      end do

      !> Get the WBOs
      allocate (wbo(mol%nat, mol%nat, wfn%nspin))
      call get_mayer_bond_orders(bas, ints%overlap, wfn%density, wbo)

      !> Save some memory by deallocating the second-iteration-specific HO overlap integrals
      deallocate (auxints%overlap_h0_2)

      if (ctx%verbosity > 1) then
         write (ctx%unit, '(10x,a,/)') "--- Two-step SCF done. ---"
      end if

   end subroutine twostepscf

   !> TAKEN OVER FROM TBLITE - modified in order to use the PTB parameters
   subroutine get_density(wfn, solver, ints, ts, error, keps_param, keps0_param)
      !> Tight-binding wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Solver for the general eigenvalue problem
      class(solver_type), intent(inout) :: solver
      !> Integral container
      type(integral_type), intent(in) :: ints
      !> Electronic entropy
      real(wp), intent(out) :: ts
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Linear regression cofactors for the orbital energies
      real(wp), intent(in), optional :: keps_param, keps0_param
      real(wp) :: keps, keps0

      real(wp) :: e_fermi, stmp(2)
      real(wp), allocatable :: focc(:)
      integer :: spin, i, j

      if (present(keps_param)) then
         keps = keps_param
      else
         keps = 0.0_wp
      end if
      if (present(keps0_param)) then
         keps0 = keps0_param
      else
         keps0 = 0.0_wp
      end if

      select case (wfn%nspin)
      case default
         ! write (*, *) "Matrix to solve ..."
         ! do i = 1, size(wfn%coeff, 1)
         !    do j = 1, size(wfn%coeff, 2)
         !       write (*, '(f11.7)', advance="no") wfn%coeff(i, j, 1)
         !    end do
         !    write (*, '(/)', advance="no")
         ! end do
         call solver%solve(wfn%coeff(:, :, 1), ints%overlap, wfn%emo(:, 1), error)
         ! write (*, *) "Coefficients after solving ..."
         ! do i = 1, size(wfn%coeff, 1)
         !    do j = 1, size(wfn%coeff, 2)
         !       write (*, '(f11.7)', advance="no") wfn%coeff(i, j, 1)
         !    end do
         !    write (*, '(/)', advance="no")
         ! end do
         if (allocated(error)) return
         wfn%emo(:, 1) = wfn%emo(:, 1) * (1.0_wp + keps) + keps0

         allocate (focc(size(wfn%focc, 1)))
         wfn%focc(:, :) = 0.0_wp
         do spin = 1, 2
            call get_fermi_filling(wfn%nel(spin), wfn%kt, wfn%emo(:, 1), &
               & wfn%homo(spin), focc, e_fermi)
            call get_electronic_entropy(focc, wfn%kt, stmp(spin))
            wfn%focc(:, 1) = wfn%focc(:, 1) + focc
         end do
         ts = sum(stmp)

         call get_density_matrix(wfn%focc(:, 1), wfn%coeff(:, :, 1), wfn%density(:, :, 1))
      case (2)
         wfn%coeff = 2 * wfn%coeff
         do spin = 1, 2
            call solver%solve(wfn%coeff(:, :, spin), ints%overlap, wfn%emo(:, spin), error)
            if (allocated(error)) return
            wfn%emo(:, 1) = wfn%emo(:, 1) * (1.0_wp + keps) + keps0

            call get_fermi_filling(wfn%nel(spin), wfn%kt, wfn%emo(:, spin), &
               & wfn%homo(spin), wfn%focc(:, spin), e_fermi)
            call get_electronic_entropy(wfn%focc(:, spin), wfn%kt, stmp(spin))
            call get_density_matrix(wfn%focc(:, spin), wfn%coeff(:, :, spin), &
               & wfn%density(:, :, spin))
         end do
         ts = sum(stmp)
      end select
      ! write (*, *) "Density matrix after solving ..."
      ! do i = 1, size(wfn%density, 1)
      !    do j = 1, size(wfn%density, 2)
      !       write (*, '(f11.7)', advance="no") wfn%density(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
   end subroutine get_density

   pure function id_to_atom(mol, idparam) result(atomparam)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> PTB parameterization data for each species
      real(wp), intent(in) :: idparam(:, :)
      !> PTB parameterization data for each atom
      real(wp), allocatable :: atomparam(:, :)

      integer :: iat, iid

      allocate (atomparam(size(idparam, 1), mol%nat), source=0.0_wp)
      do iat = 1, mol%nat
         iid = mol%id(iat)
         atomparam(:, iat) = idparam(:, iid)
      end do

   end function id_to_atom

   pure subroutine get_electronic_entropy(occ, kt, s)
      real(wp), intent(in) :: occ(:)
      real(wp), intent(in) :: kt
      real(wp), intent(out) :: s

      s = sum(log(occ**occ * (1 - occ)**(1 - occ))) * kt
   end subroutine get_electronic_entropy

end module xtb_ptb_scf
