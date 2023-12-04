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
   use tblite_context, only: context_type
   use tblite_scf_solver, only: solver_type
   use tblite_scf_iterator, only: get_qat_from_qsh
   use tblite_adjlist, only: adjacency_list, new_adjacency_list
   use tblite_cutoff, only: get_lattice_points
   use tblite_wavefunction, only: wavefunction_type, get_alpha_beta_occupation
   use tblite_wavefunction_fermi, only: get_fermi_filling
   use tblite_wavefunction_type, only: get_density_matrix
   use tblite_scf_potential, only: potential_type, new_potential, add_pot_to_h1
   use tblite_integral_type, only: new_integral, integral_type

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

   implicit none
   private

   public :: twostepscf

   real(wp), parameter :: default_cutoff = 25.0_wp

contains

   subroutine twostepscf(ctx, wfn, data, mol, bas, cbas, eeqmodel)
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
      !> Initialized EEQ model
      type(mchrg_model_type), intent(in) :: eeqmodel
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Error type
      type(error_type), allocatable :: error
      !> Coulomb potential
      type(coulomb_potential) :: coulomb
      !> +U potential
      type(plusu_potential_type) :: plusu
      !> Adjacency list
      type(adjacency_list) :: list
      !> Integral type
      type(integral_type) :: ints
      !> Auxiliary integral type
      type(aux_integral_type) :: auxints
      !> Potential type
      type(potential_type) :: pot
      !> H0 basis in second iteration
      type(basis_type), allocatable :: bas_h0
      real(wp), allocatable :: expscal_h0_2nditer(:, :)
      !> Loop variables
      integer :: i, j, isp, izp, iat, ish, iid, is
      !> Coordination numbers
      real(wp) :: cn_star(mol%nat), cn(mol%nat), cn_eeq(mol%nat)
      real(wp) :: radii(mol%nid)
      !> Effective core potential
      real(wp), allocatable :: vecp(:, :)
      !> Effective self-energies
      real(wp), allocatable :: levels(:)
      !> Lattice points
      real(wp), allocatable :: lattr(:, :)
      !> Cutoff for lattice points
      real(wp) :: cutoff
      !> Number of electrons
      real(wp) :: nel
      !> Pauli XC potential
      real(wp), allocatable :: Vxc(:, :), psh(:, :)
      !> Electronic entropy
      real(wp) :: ts

      !> Solver for the effective Hamiltonian
      call ctx%new_solver(solver, bas%nao)

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
      call get_integrals(mol, bas, lattr, list, ints%overlap, ints%dipole, norm=auxints%norm)
      !##### DEV WRITE #####
      write (*, *) "Standard overlap ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f8.5)', advance="no") ints%overlap(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      write (*, *) "Dipole:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f12.6)', advance="no") dipole(1, i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      call get_integrals(mol, lattr, list, auxints%overlap_h0, alpha_scal=id_to_atom(mol, data%hamiltonian%kalphah0l))
      !##### DEV WRITE #####
      write (*, *) "Overlap H0 scaled (SS) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      call get_integrals(mol, lattr, list, auxints%overlap_xc, alpha_scal=id_to_atom(mol, data%pauli%klalphaxc))
      !##### DEV WRITE #####
      write (*, *) "Overlap XC scaled (SS) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_xc(i, j)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################
      call get_mml_overlaps(bas, ints%overlap, ptbGlobals%mlmix, auxints%overlap_to_x, &
      & auxints%overlap_to_1_x)
      !##### DEV WRITE #####
      write (*, *) "Overlap S(1-x) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_soneminusx(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      write (*, *) "Overlap S(x) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_sx(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !     _____   _   _
      !    / ____| | \ | |
      !   | |      |  \| |
      !   | |      | . ` |
      !   | |____  | |\  |
      !    \_____| |_| \_|

      !> Get first coordination number (" CN' ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn_star)
      !##### DEV WRITE #####
      write (*, *) "CN star:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn_star(i)
      end do
      !#####################
      !> Get radii from PTB parameters for second coordination number (" CN ")
      do isp = 1, mol%nid
         izp = mol%num(isp)
         radii(isp) = rf(izp)
      end do
      !> Get second coordination number (" CN ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn, covrad=radii)
      !##### DEV WRITE #####
      write (*, *) "CN:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn(i)
      end do
      !#####################
      !> Get third coordination number for EEQ model (" CN-EEQ ")
      call ncoord_erf(mol, ptbGlobals%kerfcn_eeq, default_cutoff, cn_eeq)

      !##### DEV WRITE #####
      write (*, *) "CN-EEQ:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn_eeq(i)
      end do
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
      write (*, *) "V_ECP ..."
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
      write (*, *) "EEQ charges:"
      do i = 1, mol%nat
         write (*, '(a,i0,a,f12.6)') "Atom ", i, ":", wfn%qat(i, 1)
      end do
      !#####################
      !> Project reference occupation on wavefunction and use EEQ charges as guess
      call get_occupation(mol, bas, data%hamiltonian%refocc, wfn%nocc, wfn%n0at, wfn%n0sh)
      call guess_qsh(wfn, bas)
      !##### DEV WRITE #####
      write (*, *) "Shell populations ..."
      ! do i = 1, bas%nsh
      !    write(*,*) wfn%qsh(i, 1)
      ! enddo
      !#####################
      call coulomb%init(mol, bas, wfn%qat(:, 1), data%coulomb%shellHardnessFirstIter, &
         & data%coulomb%kQHubbard, data%coulomb%kOK1, data%coulomb%kTO)
      call coulomb%update(mol, bas)
      !##### DEV WRITE #####
      write (*, *) "Coulomb matrix ..."
      ! do i = 1, bas%nsh
      !    do j = 1, bas%nsh
      !       write (*, '(f10.6)', advance="no") coulomb%cmat(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

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
      call get_hamiltonian(mol, list, bas, data%hamiltonian, auxints%overlap_h0, &
      & levels, ints%hamiltonian, ptbGlobals%kpol, ptbGlobals%kitr, ptbGlobals%kitocod)
      ints%hamiltonian = ints%hamiltonian + vecp
      !##### DEV WRITE #####
      write (*, *) "H0 ..."
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
      call coulomb%get_potential(wfn, pot)
      call add_pot_to_h1(bas, ints, pot, wfn%coeff)
      !##### DEV WRITE #####
      write (*, *) "V_Coulomb ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f8.4)', advance="no") wfn%coeff(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      allocate (psh(bas%nsh, wfn%nspin), source=0.0_wp)
      psh = get_psh_from_qat(wfn, bas)
      call calc_Vxc_pauli(mol, bas, psh(:, 1), auxints%overlap_xc, levels, data%pauli%kxc1, wfn%coeff(:, :, 1))

      !##### DEV WRITE #####
      write (*, *) "Hamiltonian matrix to solve ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f8.4)', advance="no") wfn%coeff(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      !> Get additional potentials
      !> TODO

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

      !##### DEV WRITE #####
      write (*, *) "Coefficients after 1st iteration ..."
      ! do i = 1, size(wfn%coeff, 1)
      !    do j = 1, size(wfn%coeff, 2)
      !       write (*, '(f8.4)', advance="no") wfn%coeff(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      write (*, *) "Density matrix after 1st iteration ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.5)', advance="no") wfn%density(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      call get_mml_shell_charges(bas, auxints%overlap_to_x, auxints%overlap_to_1_x, &
         & wfn%density, wfn%n0sh, wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      !##### DEV WRITE #####
      write (*, *) "Shell charges after 1st iteration ..."
      do i = 1, bas%nsh
         write (*, '(f8.4)') wfn%qsh(i, 1)
      end do
      write (*, *) "Atom charges after 1st iteration ..."
      do i = 1, mol%nat
         write (*, '(f8.4)') wfn%qat(i, 1)
      end do
      !#####################

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
      call get_integrals(mol, bas_h0, lattr, list, auxints%overlap_h0)
      !> Deallocate temporary basis set and charge-dependent scaling factors
      deallocate (bas_h0, expscal_h0_2nditer)
      !##### DEV WRITE #####
      write (*, *) "Overlap H0 scaled (SS) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") auxints%overlap_h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      ints%hamiltonian = 0.0_wp
      call get_hamiltonian(mol, list, bas, data%hamiltonian, auxints%overlap_h0, &
      & levels, ints%hamiltonian, ptbGlobals%kpol)
      ints%hamiltonian = ints%hamiltonian + vecp
      !##### DEV WRITE #####
      write (*, *) "H0 ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.5)', advance="no") ints%hamiltonian(i, j)
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
      call add_pot_to_h1(bas, ints, pot, wfn%coeff)
      psh = get_psh_from_qsh(wfn, bas)
      call calc_Vxc_pauli(mol, bas, psh(:, 1), auxints%overlap_xc, levels, data%pauli%kxc2l, wfn%coeff(:, :, 1))
      call plusu%get_potential(mol, bas, wfn%density(:, :, 1), wfn%coeff(:, :, 1))

      !##### DEV WRITE #####
      write (*, *) "Hamiltonian matrix to solve ..."
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f9.5)', advance="no") wfn%coeff(i, j, 1)
         end do
         write (*, '(/)', advance="no")
      end do
      !#####################

      call get_density(wfn, solver, ints, ts, error, ptbGlobals%geps, ptbGlobals%geps0)

      !##### DEV WRITE #####
      write (*, *) "Coefficients after 1st iteration ..."
      ! do i = 1, size(wfn%coeff, 1)
      !    do j = 1, size(wfn%coeff, 2)
      !       write (*, '(f9.5)', advance="no") wfn%coeff(i, j, 1)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      write (*, *) "Density matrix after 2nd iteration ..."
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f9.5)', advance="no") wfn%density(i, j, 1)
         end do
         write (*, '(/)', advance="no")
      end do
      ! write (*, *) "Diagonal elements of density ..."
      ! do i = 1, bas%nao
      !    write (*, '(f9.5)') wfn%density(i, i, 1)
      ! end do
      !#####################

      call get_mml_shell_charges(bas, auxints%overlap_to_x, auxints%overlap_to_1_x, &
         & wfn%density, wfn%n0sh, wfn%qsh)
      psh = get_psh_from_qsh(wfn, bas)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      !##### DEV WRITE #####
      write (*, *) "Shell charges after 2nd iteration ..."
      do i = 1, bas%nsh
         write (*, '(f8.4)') wfn%qsh(i, 1)
      end do
      write (*, *) "Shell populations after 2nd iteration ..."
      do i = 1, bas%nsh
         write (*, '(f8.4)') psh(i, 1)
      end do
      write (*, *) "Atom charges after 2nd iteration ..."
      do i = 1, mol%nat
         write (*, '(f8.4)') wfn%qat(i, 1)
      end do
      !#####################

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
      real(wp), optional :: keps_param, keps0_param
      real(wp) :: keps, keps0

      real(wp) :: e_fermi, stmp(2)
      real(wp), allocatable :: focc(:)
      integer :: spin

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
         call solver%solve(wfn%coeff(:, :, 1), ints%overlap, wfn%emo(:, 1), error)
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
