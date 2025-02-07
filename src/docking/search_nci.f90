! This file is part of xtb.
!
! Copyright (C) 2022 Stefan Grimme, Christoph Plett
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

module xtb_docking_search_nci
   use xtb_mctc_accuracy, only: wp
   use xtb_type_environment, only: TEnvironment
   use xtb_docking_param
   use xtb_iff_iffenergy, only: iff_e, intermole_probe, alignmol
   use xtb_scanparam
   use xtb_mctc_symbols, only: toSymbol
   use xtb_type_molecule, only: TMolecule
   use xtb_type_restart, only: TRestart
   use xtb_setparam
   use xtb_gfnff_calculator, only: TGFFCalculator,newGFFCalculator
   use xtb_xtb_calculator, only: TxTBCalculator
   use xtb_type_calculator, only: TCalculator
   use xtb_main_setup, only: newCalculator
   use xtb_geoopt
   use xtb_main_defaults, only: initDefaults
   use xtb_solv_state, only: solutionState
   use xtb_dynamic, only: xyzsort2
   use xtb_gfnff_setup, only: gfnff_setup, gfnff_input
   use xtb_disp_dftd4, only: newD3Model
   use xtb_single, only: singlepoint
   use xtb_type_data, only: scc_results
   use xtb_readin, only: xfind
   use xtb_mctc_timings
   use xtb_gfnff_topology, only: TGFFTopology
   use xtb_type_topology, only: TTopology
   use xtb_gfnff_ini, only: gfnff_ini
   use xtb_topology, only: makeBondTopology
   use xtb_sphereparam, only : number_walls, cavity_egrad, maxwalls, wpot
   use xtb_gfnff_param, only: ini, gfnff_set_param, gfnff_param_alloc,&
   & gff_print, gfnff_param_dealloc
   use xtb_constrain_param, only: read_userdata
   use xtb_fixparam
   use xtb_disp_ncoord, only: ncoord_gfn, ncoord_erf, ncoord_d3
   use xtb_scc_core, only: iniqshell
   use xtb_eeq, only: goedecker_chrgeq
   use xtb_basis, only: newBasisset
   use xtb_gfnff_neighbor, only: TNeigh
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_filetypes, only : generateFileName
   use xtb_iniq, only: iniqcn
   implicit none

   private
   public :: docking_search

contains

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! find best structure globally by genetic algorithm
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   subroutine docking_search(env,molA,molB,n,n1,n2,at1,at2,neigh,xyz1,xyz2,q1,q2,c6ab,&
                       &        z1, z2, nl1, nl2, l1, l2, cl1, cl2, qdr1, qdr2,&
                       &        cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                       &        den1, den2, gab1, gab2, molA_e, molB_e, cprob, emin, coord, comb)

      INTERFACE
         RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices)
            REAL, DIMENSION(:), INTENT(INOUT) :: Item       ! array of values
            INTEGER, INTENT(IN)    :: First, Last
            INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
         END SUBROUTINE Quicksort
      END INTERFACE

      !> Source of errors in the main program unit
      character(len=*), parameter :: source = "docking_search"
      type(TEnvironment), intent(inout) :: env
      integer, intent(in) :: n, n1, n2
      integer, intent(in) :: nl1, nl2
      integer, intent(in) :: neigh(0:n, n)
      real(wp), intent(inout) ::molA_e, molB_e
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)
      real(wp), intent(in) :: cl1(4, n1*10)
      real(wp), intent(in) :: cl2(4, n2*10)
      real(wp), intent(in) :: q1(n1)
      real(wp), intent(in) :: q2(n2)
      real(wp), intent(in) :: c6ab(n, n)
      real(wp), intent(in) :: alpab(n2, n1)
      real(wp), intent(in) :: cprob(n1)
      real(wp), intent(in) :: emin
      real(wp), intent(in) :: coord(6)
      real(wp), intent(in) :: qdr1(n1)
      real(wp), intent(in) :: qdr2(n2)
      real(wp), intent(in) :: z1(n1)
      real(wp), intent(in) :: cn1(n1)
      real(wp), intent(in) :: z2(n2)
      real(wp), intent(in) :: cn2(n2)
      real(wp), intent(in) :: alp1(n1)
      real(wp), intent(in) :: alp2(n2)
      real(wp), intent(in) :: qct1(n1, 2), qct2(n2, 2)
      real(wp), intent(in) :: den1(2,4,n1),den2(2,4,n2),gab1(n1,n1),gab2(n2,n2)
      integer, intent(in) :: at1(n1)
      integer, intent(in) :: at2(n2)
      integer, intent(in) :: l1(n1*10)
      integer, intent(in) :: l2(n2*10)
      type(TMolecule), intent(inout) :: molA, molB
      type(TMolecule), intent(inout) :: comb

      real(wp) :: E_int
      real(wp), allocatable :: found(:, :), found2(:, :)
      real(wp), allocatable :: xyzprobc(:, :)
      real(wp), allocatable :: xyzprob3(:, :), xyzprob4(:, :)
      real(wp), allocatable :: xyztmp(:, :), rprobc(:, :)
      real(wp), allocatable :: xyzprob0(:), xyzprob1(:), xyzprob2(:), rprob(:)
      integer, allocatable :: ind0(:), ind1(:), ind2(:), rind(:), attmp(:)
      integer :: counter
      real(wp) :: found_tmp(6), displ(6), c2(3, n2), f(6), bestsofar
      real(wp) :: r,e,dx,dy,dz,xx1,yy1,zz1,bord,bord_rg,stepr3(3),stepr4,eshift
      real(wp) :: dx_rg, dy_rg, dz_rg, xx1_rg, yy1_rg, zz1_rg
      real(wp) :: r1, r2, r3, a1, a2, a3, av, sig, epq, tmpx(15), t0, t1, w0, w1
      real(wp) :: dum3(3, 3), dum2(3), dum(3), pmass, rmin, mvec(3, 14)
      real(wp) :: rmsd(n_opt*(n_opt + 1)/2)
      real(wp) :: dx2, dy2, dz2
      real*4   :: r4
      integer :: i, j, k, nang, nr, ii, jj, ns, icycle, maxs, iseed(2)
      integer :: kk = 0
      integer :: m1, m2, m3, l, nn, ndim, ll, kkk, kkkk, mfrag, maxlow, nclust, ierr
      integer :: mm1, mm2, mm3, maxpock, icssym
      integer, allocatable :: atfrag(:), molvec(:), iatnum(:)
      real(wp), allocatable :: xyzfrag(:, :)
      logical :: tooclose
      integer :: maxstruc
      character(len=100) :: warning
      real(wp) :: xyz_tmp(3, molB%n)
      integer, external :: ncore
      !> GFN-FF Calculator
      class(TCalculator), allocatable :: calc
      type(TGFFCalculator) :: gff_calc
      type(TGFFTopology) :: topo_backup, topo_xTB
      type(TNeigh) :: neigh_backup
      type(TTopology) :: bonds 
      !> Parameterfile
      character(len=:), allocatable :: pfile
      character(len=:), allocatable :: fnv
      character(len=*), parameter :: p_fname_param_gfn0 = 'param_gfn0-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn1 = 'param_gfn1-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn2 = 'param_gfn2-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfnff = '.param_gfnff-xtb'
      !> Restart
      type(TRestart) :: chk
      !> SCC results
      type(scc_results) :: results
      !> Opt variable
      real(wp) :: egap, et, etot, grad(3, n), sigma(3, 3)
!      integer :: maxiter, maxcycle_in,olev
      logical :: pr, initial_sp, fail
      logical :: restart = .false.
      !> Dummy
      integer, allocatable ::  tmp_unit
      integer :: itemp = 64, iopt = 96, iprob = 30, ifinal = 31
      integer :: itopo = 32, ipocket = 53, ichr = 33, iscreen = 34
      !> Pocket optimization
      real(wp), allocatable :: xyz_pocket(:, :, :)
      real(wp), allocatable :: pocket_e(:)
      !> Structure after Optimization
      real(wp), allocatable :: xyz_opt(:, :, :)
      !> Final energies
      real(wp), allocatable :: final_e(:)
      real(wp) :: tmp_e
      !> Minimumposition
      integer :: minpos
      !> For outprint
      character(len=:),allocatable :: extension, fin_name

      !> Parameter
      real(wp), parameter :: pi2 = 2*3.14159265358979_wp
      real(wp), parameter :: pi = 3.14159265358979_wp
      real(wp), parameter :: au = 627.509541_wp
      real(wp), parameter ::autoang = 0.52917726_wp
      character*80 atmp
      data mvec/1., 0., 0., 0., 1., 0., 0., 0., 1.,&
                 &-1., 0., 0., 0, -1., 0., 0., 0., -1.,&
                 & 1., 1., 1., -1., -1., -1.,&
                 &-1., 1., 1., 1., -1., 1., 1., 1., -1.,&
                 & 1., -1., -1., -1., 1., -1., -1., -1., 1./

      call start_timing(4)

      write (env%unit, *)
      write (env%unit, *) '=============================================='
      write (env%unit, *) '|         Starting Energy Screening          |'
      write (env%unit, *) '=============================================='
      write (env%unit, *)
      if (.not. fulle) write (env%unit, *) ' No ATM term employed (recommended)'
      if (.not. fulle) write (env%unit, *) ' If ATM term should be included, use -atm option.'
      if (.not. fulle) write (env%unit, *)

!     just output
      call axis_docking(n1, at1, xyz1, dum, pmass, dum2, dum3)
      r1 = sqrt(dum2(1) + dum2(2) + dum2(3))
      call axis_docking(n2, at2, xyz2, dum, pmass, dum2, dum3)
      r2 = sqrt(dum2(1) + dum2(2) + dum2(3))
      call rcma(n1, xyz1, at1, n2, xyz2, at2, r, rmin)
      write (env%unit, '('' Method for final opts.    :'',1x,a  )') optlvl
      write (env%unit, '('' # of genetic optimizations:'',1x,i0  )') maxgen
      write (env%unit, '('' # of parents              :'',1x,i0  )') maxparent
      write (env%unit, '('' # of final geo. opts.     :'',1x,i0  )') n_opt
      write (env%unit, '('' Rare gas grid step size   :'',F8.2)') stepr
      write (env%unit, '('' ang step size /deg        :'',F8.2)') stepa
      write (env%unit, '('' # angular grid points     :'',1x,i0  )') (360/int(stepa))**3
      write(*,*)
      if (pocket_grid) write (env%unit, '('' Performing pocket search'')')
      if (stack_grid) write (env%unit, '('' Performing stack search'')')
      if (angular_grid)  write (env%unit, '('' Performing angular search'')')
      write(*,*)

      !Calculator prep
      call combine_mol(comb, molA, molB) !molB is shifted far away from molA

      !Setting up z
      do i = 1, comb%n
         comb%z(i) = comb%at(i) - ncore(comb%at(i))
         ! lanthanides without f are treated as La
         if (comb%at(i) .gt. 57 .and. comb%at(i) .lt. 72) comb%z(i) = 3
      end do

      call open_file(itemp, 'opt_tmp', 'w')
      tmp_unit = env%unit
      env%unit = itemp
      if (optlvl == 'gfnff') fnv = xfind(p_fname_param_gfnff)
      if (optlvl == 'gfn2') fnv = xfind(p_fname_param_gfn2)
      if (optlvl == 'gfn0') fnv = xfind(p_fname_param_gfn0)
      if (optlvl == 'gfn1') fnv = xfind(p_fname_param_gfn1)
      call newCalculator(env, comb, calc, fnv, restart, acc)
      call initDefaults(env, calc, comb, gsolvstate_iff)
      write(*,*) 'initialization done'

      select type (calc)
      type is (TGFFCalculator)
         topo_backup = calc%topo !topo is from molA and molB that are far away
         neigh_backup = calc%neigh
      end select
      pr = .false.
      initial_sp = .true.
      !>Silent GFN calculations (No charge file or anything is written)
      set%silent = .true.
      gff_print = .false.
      set%pr_finalstruct = .false.

      !> SP molA
      set%ichrg = nint(molA%chrg)
      select type (calc)
      type is (TGFFCalculator)
         !With p_modh_old the next step craches sometimes for really large systmes
         call restart_gff(env, molA, calc)
         call calc%singlepoint(env,molA,chk,1,.false.,etot,grad,sigma,egap,results)
         do i = 1, molA%n
            topo_backup%qa(i) = calc%topo%qa(i) !Done to ensure right charge assignment to molecules
         end do
         molA_e = etot
      type is (TxTBCalculator)
         molA_e = pre_e_A
      end select

      !> SP molB
      set%ichrg = nint(molB%chrg)
      select type (calc)
      type is (TGFFCalculator)
         call restart_gff(env, molB, calc)
         call calc%singlepoint(env,molB,chk,1,.false.,etot,grad,sigma,egap,results)
         do i = 1, molB%n
            topo_backup%qa(i + molA%n) = calc%topo%qa(i)!Charges from molB for right partial charges
         end do
         molB_e = etot

         !> Topo of both single molecules together
      type is (TxTBCalculator)
         molB_e = pre_e_B
      end select

      etot = 0.0_wp; grad = 0.0_wp; sigma = 0.0_wp; egap = 0.0_wp

      env%unit = tmp_unit
      call close_file(itemp)
      call remove_file(itemp)

      write (env%unit, '(''  Total '', a,'' energy molecule 1: '', F16.10)') optlvl, molA_e
      write (env%unit, '(''  Total '', a,'' energy molecule 2: '', F16.10)') optlvl, molB_e

      !> Setting charge back to combined ones
      set%ichrg = nint(comb%chrg)

      write (env%unit, *)
      write (env%unit, *) '-----------------------------'
      write (env%unit, *) ' Grid based energy screening '
      write (env%unit, *) '-----------------------------'
      write (env%unit, *)

      !>Getting borders
      xx1 = maxval(xyz1(1, 1:n1)) - minval(xyz1(1, 1:n1))
      yy1 = maxval(xyz1(2, 1:n1)) - minval(xyz1(1, 1:n1))
      zz1 = maxval(xyz1(3, 1:n1)) - minval(xyz1(1, 1:n1))
      dum(1) = xx1
      dum(2) = yy1
      dum(3) = zz1
      dx2 = maxval(xyz2(1, 1:n2)) - minval(xyz2(1, 1:n2)) !On every side, molB has to have space
      dy2 = maxval(xyz2(2, 1:n2)) - minval(xyz2(2, 1:n2))
      dz2 = maxval(xyz2(3, 1:n2)) - minval(xyz2(3, 1:n2))
      dum2(1) = dx2
      dum2(2) = dy2
      dum2(3) = dz2
      bord = maxval(abs(dum)) + maxval(dum2) + 2 !Max bond distance of 2 considered
      bord_rg = 30 !Additional border increase for RG scan

      !> Angular search grid
      dx = maxval(xyz1(1, 1:n1)) - minval(xyz1(1, 1:n1)) + bord
      dy = maxval(xyz1(2, 1:n1)) - minval(xyz1(2, 1:n1)) + bord
      dz = maxval(xyz1(3, 1:n1)) - minval(xyz1(3, 1:n1)) + bord
      xx1 = minval(xyz1(1, 1:n1)) - bord*0.5
      yy1 = minval(xyz1(2, 1:n1)) - bord*0.5
      zz1 = minval(xyz1(3, 1:n1)) - bord*0.5
      stepr3(1) = dx/3
      stepr3(2) = dy/3
      stepr3(3) = dz/3
      mm1 = 4
      mm2 = 4
      mm3 = 4 !Grid with 256 points

      !> RG Grid (larger to have more points for pocket search, cheap)
      dx_rg = dx + bord_rg
      dy_rg = dy + bord_rg
      dz_rg = dz + bord_rg
      xx1_rg = xx1 - bord_rg*0.5
      yy1_rg = yy1 - bord_rg*0.5
      zz1_rg = zz1 - bord_rg*0.5
     !> Increase Gridsize for RG-Scan additionally to have enough points for pocket search
      m1 = 1 + int((dx_rg)/stepr)
      m2 = 1 + int((dy_rg)/stepr)
      m3 = 1 + int((dz_rg)/stepr)

      icssym = 0
      if (cssym) then
         dx2 = dx - bord
         dy2 = dy - bord
         dz2 = dz - bord
         if (dx2 .lt. dy2 .and. dx2 .lt. dz2) then
            m1 = m1/2
            mm1 = mm1/2
            xx1 = xx1 + bord*0.5
            icssym = 1
         end if
         if (dy2 .lt. dx2 .and. dy2 .lt. dz2) then
            m2 = m2/2
            mm2 = mm2/2
            yy1 = yy1 + bord*0.5
            icssym = 2
         end if
         if (dz2 .lt. dx2 .and. dz2 .lt. dy2) then
            m3 = m3/2
            mm3 = mm3/2
            zz1 = zz1 + bord*0.5
            icssym = 3
         end if
         write (*, '(''  # Cs symmetry red. :'',i1  )') icssym
      end if

      nn = m1*m2*m3
      if (maxparent < 28) then
         maxparent = 28 !Otherwise the computation crashes as more structures than parents are present
         call env%warning('Too few maxparent, taking 28', source)
      end if
      maxpock = 10 ! max number of pockets
      ndim = maxparent**2 !size of gene pool
      nclust = 0

      allocate (found(7, ndim), found2(7, ndim))  ! main arrays for gen ensemble
      found = 0
      found2 = 0

      write (*, '(''  # probe RG points   :'',i0  )') nn
      allocate (xyzprob0(nn), xyzprob1(nn),&
              &xyzprob2(nn), xyzprobc(3, nn),&
              &ind0(nn), ind1(nn), ind2(nn))
      xyzprob0 = 0
      xyzprob1 = 0
      xyzprob2 = 0

!ccccccccccccc
! R grid for RG atom
!ccccccccccccc
      dx_rg = xx1_rg
      l = 0
      do i = 1, m1
         dum(1) = dx_rg
         dy_rg = yy1_rg
         do j = 1, m2
            dum(2) = dy_rg
            dz_rg = zz1_rg
            do k = 1, m3
               dum(3) = dz_rg
               l = l + 1
               xyzprobc(1:3, l) = dum(1:3)
               dz_rg = dz_rg + stepr
            end do
            dy_rg = dy_rg + stepr
         end do
         dx_rg = dx_rg + stepr
      end do

      if (debug) then
        call wr_grid('RG_grid.xyz', molA%n, nn, molA%at, 36, molA%xyz, xyzprobc)
      end if

!> Raregas atom check around MolA
!$omp parallel do default(none) &
!$omp shared(xyzprob0,xyzprob1,xyzprob2,n1,nl1,at1,q1,xyz1,cl1,l1,cn1,cprob,nn,xyzprobc) &
!$omp private(l,e,epq,dum)
      do l = 1, nn
         dum(1:3) = xyzprobc(1:3, l)
         call intermole_probe(n1, nl1, at1, q1, xyz1, cl1,&
                             &dum, l1, cn1, cprob,&
                             &e, epq)
         xyzprob0(l) = e
         xyzprob1(l) = e + epq
         xyzprob2(l) = e - epq
      end do
!$omp end parallel do

      do i = 1, nn
         ind0(i) = i
         ind1(i) = i
         ind2(i) = i
      end do
      call Qsort(xyzprob0, 1, nn, ind0) ! sort for E neutral
      call Qsort(xyzprob1, 1, nn, ind1) ! sort for probe + 0.1 charge
      call Qsort(xyzprob2, 1, nn, ind2) ! sort for probe - 0.1 charge

      if (debug) then
        call wr_grid('best_RG_position.xyz', molA%n, 1, molA%at, 36, molA%xyz, xyzprobc(1:3, ind0(1)))
      end if

      write (*, '(''  Best rare gas probe energy/kcal   :'',F8.2)') xyzprob0(1)*au
      write (*, '(''  +0.1 charged probe energy/kcal:'',F8.2)') xyzprob1(1)*au
      write (*, '(''  -0.1 charged probe energy/kcal:'',F8.2)') xyzprob2(1)*au

!cccccccccccccccccccccccccccccccccccccccccc
! clustering based on RG filled pockets
!cccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------- start pocket search
      if (pocket_grid) then
         write(*,*)
         write (*,*) '  Starting pocket search'
         write(*,*)

         if (nn/4 < mxcent_clust) mxcent_clust = nn/4 !Default mxcent_clust=500

         allocate(atfrag(mxcent_clust),molvec(mxcent_clust),iatnum(mxcent_clust),xyzfrag(3,mxcent_clust))

         do i = 1, mxcent_clust
            xyzfrag(1:3, i) = xyzprobc(1:3, ind0(i))
            atfrag(i) = probe_atom_type
         end do
         molvec = 0
         call mrec(mfrag, xyzfrag, mxcent_clust, atfrag, molvec) ! recursive mol splitting into non-bonded fragments
         ! mfrag is the number of fragments that the best 500 gridpoints xyzfrag were split into
         ! molvec tells to which nci fragment the respective gridpoint belongs to
         ! if all scanned gridpoints are close, only one fragment results and every value of molvec will be 1
         iatnum = 0
         do i = 1, mxcent_clust
            k = molvec(i)
            iatnum(k) = iatnum(k) + 1
         end do

         ! now generate fragment geom
         kk = 0
         do i = 1, mfrag
            k = iatnum(i)
            allocate (xyztmp(3, k), attmp(k))
            xyztmp = 0.0d0
            attmp = 0
            l = 0
            do j = 1, mxcent_clust
               if (molvec(j) .eq. i) then !gridpoint j belongs to fragment number i
                  l = l + 1
                  xyztmp(1:3, l) = xyzfrag(1:3, j)
                  attmp(l) = atfrag(j)
               end if
            end do
            if (debug) then
               call wr_grid('pocket_grid.xyz', molA%n, k, molA%at, 36, xyz1, xyztmp)
            end if

            ! xyztmp contains now each gridpoint belonging to fragment i
            ! same for attmp (anyway always probe_atom_type)
            ! determine size (=R) of gridpoint cluster belonging to fragment i
            call axis_docking(k, attmp, xyztmp, dum, pmass, dum2, dum3)
            r = sqrt(dum2(1) + dum2(2) + dum2(3))

            !> check if moleculeB is small enough to fit in gridpoint cluster of fragmen i
            if (r/r2 .gt. 0.5 .and. k .gt. 10) then  
               displ = 0
               !> Align xyz2 into gridpoint cluster number i by placing it in the middle and
               !  performing random mutations on the xyz coords and the rotation angles
               !  Each position is scanned by dispersion energy and the lowest one is given back as displ
               call alignmol(k, n2, attmp, at2, xyztmp, xyz2, displ)
               kk = kk + 1
               found2(1:6, kk) = displ(1:6)
            end if
            deallocate (xyztmp, attmp)
         end do
         write (*, '(''# pocket optimizations  :'',2i5  )') kk
         allocate (xyz_pocket(3, comb%n, kk))
         allocate (pocket_e(kk))

         !> Prepare optimization
         set%opt_logfile = 'pocket_opt.xyz'
         call open_file(ipocket, 'pocket.xyz', 'w')
         call open_file(itemp, 'opt_tmp', 'w')
         tmp_unit = env%unit
         env%unit = itemp

         do i = 1, kk
            displ(1:6) = found2(1:6, i)
            found_tmp = found2(1:6, i)
            call move2(molB%n, molB%xyz, xyz_tmp, found_tmp) !return xyz_tmp that is with found_tmp transformed molB%xyz
            counter = 0
            do j = 1, molA%n
               comb%xyz(1:3, j) = molA%xyz(1:3, j) !comb overwritten with A, as it is changed upon geo_opt
            end do
            do j = molA%n + 1, molA%n + molB%n
               counter = counter + 1
               comb%xyz(1:3, j) = xyz_tmp(1:3, counter) !combined molA and shifted molB
            end do
            !> Select which kind of optimization is done
            select type (calc)
            type is (TGFFCalculator)
               call restart_gff(env, comb, calc)
               !Keeping Fragments and charges
               calc%neigh%nbond = neigh_backup%nbond
               calc%neigh%nb = neigh_backup%nb
               calc%topo%qfrag = topo_backup%qfrag
               calc%topo%qa = topo_backup%qa
               calc%topo%fraglist = topo_backup%fraglist
               calc%topo%nfrag = topo_backup%nfrag
            type is (TxTBCalculator)
               call restart_xTB(env, comb, chk, calc)
            end select
            call geometry_optimization&
            &     (env, comb, chk, calc, egap, set%etemp, set%maxscciter, set%optset%maxoptcycle,&
            &      etot, grad, sigma, set%optset%optlev, pr, initial_sp, fail)
            grad = 0.0_wp; sigma = 0.0_wp; egap = 0.0_wp
            do k = 1, 3
               do j = 1, comb%n
                  xyz_pocket(k, j, kk) = comb%xyz(k, j)
               end do
            end do
            pocket_e(i) = etot
            write (*, '(''pocket '',i3,'' E /au :'',f9.2)') i, etot
            write (ipocket, '(i0)') comb%n
            write (ipocket, '(f20.14)') etot
            do j = 1, comb%n
               write (ipocket, '(a4,2x,3f20.14)') comb%sym(j), comb%xyz(:, j)*autoang
            end do
         end do

         call close_file(ipocket)
         call remove_file(itemp)
         call delete_file(set%opt_logfile)
         env%unit = tmp_unit
         if (.not. stack_grid .AND. .not. angular_grid) then
            call stop_timing(4)
            return !Than nothing has to be done anymore
         end if
      end if
!---------------------------------------------------------- end pocket search

!---------------------------------------------------------- start stack grid
      if (stack_grid) then
         write(*,*)
         write (*,*) '  Starting stack search'

!        move along axis seperately. This solves the identical molecule/stack case
!        mxcma= 20 !Setting for nice pictures
         ii = mxcma   ! SETING
         nn = 4*14*ii ! all directions incl. diagonal axis, none+ three rotsations of molB on each position
         if (cssym) nn = 4*9*ii
         write (*, '(''   Grid points: '',i0  )') nn
         allocate (rprobc(6, nn), rprob(nn), source=0.0_wp)
         allocate (rind(nn), source=0)
         stepr4 = 0.2*float(n1 + n2)**0.10d0 ! increase step for larger systems
!        stepr4=2 !Setting for nice pictures

         l = 0
         do j = 1, 14
            if (cssym .and. mvec(icssym, j) .lt. 0) cycle ! exclude sym. equiv.
            r = stepr4
            if (j .gt. 6) r = stepr4/sqrt(3.)
            dum2(1:3) = mvec(1:3, j)*r !mvec just vector in every direction in 3D
            do i = 1, ii !Moving stepwise along one vector of mvec (one direction in 3D)
               displ = 0
               displ(1:3) = float(ii + 1 - i)*dum2(1:3) !Point in distance ii-i in every direction
               l = l + 1
               rprobc(1:6, l) = displ(1:6)
               do k = 4, 6
                  displ(k) = pi !Rotationangles = Pi
                  l = l + 1
                  rprobc(1:6, l) = displ(1:6)
               end do
            end do
         end do

         !$omp parallel do default(none) reduction(+:rprob) &
         !$omp shared(env,nn,n,n1,n2,at1,at2,neigh,xyz1,xyz2,q1,q2,c6ab,&
         !$omp z1,z2,nl1,nl2,l1,l2,cl1,cl2,qdr1,qdr2,cn1,cn2,alp1,&
         !$omp alp2,alpab,qct1,qct2,den1,den2,gab1,gab2,rprobc) &
         !$omp private(l,e,displ)
         do l = 1, nn
            displ(1:6) = rprobc(1:6, l)
            call iff_e(env, n, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab,&
                          &z1, z2, nl1, nl2, l1, l2, cl1, cl2,&
                          &qdr1, qdr2,&
                          &cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                          &den1, den2, gab1, gab2,&
                          &.false., -1, e, displ)
            rprob(l) = e
         end do
         !$omp end parallel do

         if (debug) then
            call wr_grid('stack_grid.xyz', molA%n, nn, molA%at, &
                    &    36, xyz1, rprobc(1:3, 1:nn))
         end if

         do i = 1, nn
            rind(i) = i
         end do
         call Qsort(rprob, 1, nn, rind) ! sort for E neutral
         write (*, '(''   lowest found /kcal :'',F8.2)') rprob(1)*au
         do i = 1, nn/100 ! take some lowest
            if (rprob(i) .lt. 0.1) then
               nclust = nclust + 1
               found(1:6, nclust) = rprobc(1:6, rind(i))
               found(7, nclust) = rprob(rind(i))
            end if
         end do
         ll = nclust

         if(debug) then
            if (rprob(1)*au .lt. bestsofar) then
               call wrc3('best_after_angular_stack.xyz', n1, n2, at1, at2, xyz1, xyz2, 0, found)
               bestsofar = rprob(1)*au
            end if

            deallocate (rprobc, rprob, rind)
         endif

      end if
!--------------------------------------------------------------- end stack grid

!--------------------------------------------------------------- start angular rot
      if (angular_grid) then
         write(*,*)
         write (*, *) '  Starting angular search'

         m1 = 360/int(stepa)       ! angular grid with stepa deg steps
         maxlow = mm1*mm2*mm3 + 2    ! crude R grid dim
         nn = maxlow*m1*m1*m1
         write (*, '(''   Grid points:'',i0  )') nn
         allocate (xyzprob3(3, maxlow), xyzprob4(7, nn), source=0.0_wp)

!        include the best +/- polar points for the angular grid
         xyzprob3(1:3, 1) = xyzprobc(1:3, ind1(1)) ! take only lowest
         xyzprob3(1:3, 2) = xyzprobc(1:3, ind2(1)) !  "    "     "

         deallocate (xyzprobc, xyzprob0, xyzprob1, xyzprob2, ind0, ind1, ind2)

!        just a very crude grid over entire molecule 1 for filling up points
         dx = xx1
         l = 2 !The first two are best +/- polarized points
         do i = 1, mm1
            dum(1) = dx
            dy = yy1
            do j = 1, mm2
               dum(2) = dy
               dz = zz1
               do k = 1, mm3
                  dum(3) = dz
                  l = l + 1
                  xyzprob3(1:3, l) = dum(1:3)
                  dz = dz + stepr3(3) !stepr3=10*stepr
               end do
               dy = dy + stepr3(2)
            end do
            dx = dx + stepr3(1)
         end do
         stepa = pi2/float(m1 - 1)

         if (debug) then
            call wr_grid('angular_grid.xyz', molA%n, maxlow, molA%at, 36, xyz1, xyzprob3)
         end if

!        loop over R points
         !$omp parallel do default(none) &
         !$omp shared(env,xyzprob4,maxlow,m1,stepa,xyzprob3,n,n1,n2,at1,at2,neigh,&
         !$omp xyz1,xyz2,q1,q2,c6ab,z1,z2,nl1,nl2,l1,l2,cl1,cl2,&
         !$omp qdr1,qdr2,cn1,cn2,alp1,alp2,alpab,qct1,qct2,den1,den2,gab1,gab2)&
         !$omp private(icycle,i,j,k,e,displ,l)
         do icycle = 1, maxlow !Maxlow is dimension of R grid with stepr3
            l = 0
            !> Rotating at each crude gridpoint moleculeB stepwise
            do i = 1, m1
               displ(4) = (i - 1)*stepa
               do j = 1, m1
                  displ(5) = (j - 1)*stepa
                  do k = 1, m1
                     displ(6) = (k - 1)*stepa
                     displ(1:3) = xyzprob3(1:3, icycle)
                     call iff_e(env, n, n1, n2, at1, at2, neigh,&
                               &xyz1, xyz2, q1, q2, c6ab, z1, z2,&
                               &nl1, nl2, l1, l2, cl1, cl2,&
                               &qdr1, qdr2,&
                               &cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                               &den1, den2, gab1, gab2,&
                               &.false., -1, e, displ)
                     l = k + m1*(j - 1) + m1**2*(i - 1) + (icycle - 1)*maxlow
                     xyzprob4(1:6, l) = displ(1:6)
                     xyzprob4(7, l) = e
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         call sort6(nn, xyzprob4)
!--------------------------------------------------------------------- end angular rot

         !> Adding stack and angular search
         ll = 5*maxparent + nclust

         !> Performing an energy computation of stack and angular search structures
         !> to rank them energetically
         !$omp parallel do default(none) &
         !$omp shared(env,found,ll,xyzprob4,nclust,n,n1,n2,at1,at2,neigh,xyz1,&
         !$omp xyz2,q1,q2,c6ab,z1,z2,nl1,nl2,l1,l2,cl1,cl2,qdr1,qdr2,cn1,cn2,&
         !$omp alp1,alp2,alpab,qct1,qct2,den1,den2,gab1,gab2) &
         !$omp private(i,displ,e)
         do i = 1, ll                     ! approx nparent*2 because doubles appear
            if (i .le. nclust) then
               displ(1:6) = found(1:6, i) !Lowest of stack search
            else
               displ(1:6) = xyzprob4(1:6, i - nclust) !Lowest of angular search
            end if
            call iff_e(env, n, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab,&
                          &z1, z2, nl1, nl2, l1, l2, cl1, cl2,&
                          &qdr1, qdr2,&
                          &cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                          &den1, den2, gab1, gab2,&
                          &.false., -1, e, displ)
            found(1:6, i) = displ(1:6)
            found(7, i) = e*au
         end do
         !$omp end parallel do

         deallocate (xyzprob3, xyzprob4)
      end if

      call doubles(ll, ndim, 3.0d0, 0.1d0, found)    ! Sorting out doubles
      call sort6(ll, found)
      if(debug) call wrc3('best_before_gen.xyz', n1, n2, at1, at2, xyz1, xyz2, 1, found)

      write (*, *)
      write (*, *) '  Interaction energy of lowest structures so far in kcal/mol:'
      do i = 1, 10
         write (*, '(   F14.2)') found(7, i)
      end do
      write(*,*)

!      write (*, *) '(Rx,Ry,Rz,alp,bet,gam,Eint in kcal/mol):'
!      do i = 1, 10
!         write (*, '(6F8.3,5x,F14.2)') found(1:7, i)
!      end do

      if (debug) then
         call wrc2('genstart.xyz', 1, n1, n2, at1, at2, xyz1, xyz2,&
                   &maxparent, found)
      end if

!-------------------------------------------------------------- start genetic optimization

      write (*, *) '------------------------------'
      write (*, *) 'genetic optimization algorithm'
      write (*, *) '------------------------------'
      write (*, *) '  cycle  Eint/kcal/mol  average Eint'
      do icycle = 1, maxgen

         !$omp parallel do default(none) &
         !$omp shared(env,found2,maxparent,found,n,n1,n2,at1,at2,neigh,&
         !$omp xyz1,xyz2,q1,q2,c6ab,z1,z2,nl1,nl2,l1,l2,cl1,cl2,qdr1,qdr2,cn1,cn2,&
         !$omp alp1,alp2,alpab,qct1,qct2,den1,den2,gab1,gab2) &
         !$omp private(i,j,ii,displ,e,f)
!     LOOP HEAD -------------------------------------------------------
         do i = 1, maxparent
            do j = 1, maxparent
               ii = (i - 1)*maxparent + j
               call crossover(0.00d0, f)
               displ(1:6) = found(1:6, i)*f(1:6) + found(1:6, j)*(1.0d0 - f(1:6))
               if (i .ne. j) call rand6(0.5d0, 1.0d0, 1.0d0, displ)   ! mutation only on childs
               call iff_e(env, n, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab,&
                               &z1, z2, nl1, nl2, l1, l2, cl1, cl2,&
                               &qdr1, qdr2,&
                               &cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                               &den1, den2, gab1, gab2,&
                               &.false., -1, e, displ)            ! only e returned
               found2(1:6, ii) = displ(1:6)
               found2(7, ii) = e*au
            end do
         end do
         !$omp end parallel do

         ii = maxparent**2
         call doubles(ii, ndim, 3.0d0, 0.10d0, found2)          ! SETING
         call sort6(ii, found2)

         found = found2
         if (debug) call wrc2('structures_after_gen.xyz', 0, n1, n2, at1, at2, xyz1,&
         & xyz2, n_opt, found)
         call wrc3('best_after_gen.xyz', n1, n2, at1, at2, xyz1, xyz2, 1, found)

         av = sum(found(7, 1:maxparent))/float(maxparent)
         sig = 0
         do i = 1, n_opt
            sig = sig + (found(7, i) - av)**2
         end do
         sig = sqrt(sig/float(n_opt - 1))
         write(*,'(4x,i0,7x,F7.1,5x,F8.1,5x)')&
         &icycle, found(7, 1), av
         if (sig .lt. 0.3d0) exit
!     LOOP END ------
      end do
      call stop_timing(4)
!-------------------------------------------------------------- end genetic optimization

      !> Check if ensemble runtype is requested
      if (docking_ens) then
         n_opt = 0
         do i = 1, ii
            !> Increase number of optimizations for each structure with an Eint of lower than -0.1 kcal/mol
            if (found(7, i) < -0.1) then
               n_opt = n_opt + 1
            else
               exit
            end if
         end do
      end if

!-------------------------------------------------------------- final optimization

      call open_file(iopt, 'optimized_structures.xyz', 'w')
      call open_file(itemp, 'opt_tmp', 'w')
      tmp_unit = env%unit
      env%unit = itemp
      set%opt_logfile = 'final_opt.xyz'

      if (pocket_grid) then !include pocket search
         allocate (xyz_opt(3, n, n_opt + kk), source=0.0_wp)
         allocate (final_e(n_opt + kk), source=0.0_wp)
      else
         allocate (xyz_opt(3, n, n_opt), source=0.0_wp)
         allocate (final_e(n_opt), source=0.0_wp)
      end if

      write (*, *)
      write (*, '(''Optimizing '',i0,'' best structures with '',a )') n_opt, optlvl
      do icycle = 1, n_opt
         found_tmp = found2(1:6, icycle)
         call move2(molB%n, molB%xyz, xyz_tmp, found_tmp) !return xyz_tmp that is with found_tmp transformed molB%xyz
         do j = 1, molA%n
            comb%xyz(1:3, j) = molA%xyz(1:3, j) !comb overwritten with A, as it is changed upon geo_opt
         end do

         counter = 0
         do j = molA%n + 1, molA%n + molB%n
            counter = counter + 1
            comb%xyz(1:3, j) = xyz_tmp(1:3, counter) !combined molA and shifted molB
         end do

         select type (calc)
         type is (TGFFCalculator)
            call restart_gff(env, comb, calc)
            calc%neigh%nbond = neigh_backup%nbond
            calc%neigh%nb = neigh_backup%nb
            calc%topo%qfrag = topo_backup%qfrag
            calc%topo%qa = topo_backup%qa
            calc%topo%fraglist = topo_backup%fraglist
            calc%topo%nfrag = topo_backup%nfrag
         type is (TxTBCalculator)
            call restart_xTB(env, comb, chk, calc)
         end select
         write (*, *) icycle
         call start_timing(5)
         call geometry_optimization &
         &     (env, comb, chk, calc, egap, set%etemp, set%maxscciter, set%optset%maxoptcycle,&
         &      etot, grad, sigma, set%optset%optlev, pr, initial_sp, fail)
         grad = 0.0_wp; sigma = 0.0_wp; egap = 0.0_wp
         call stop_timing(5)
         do i = 1, 3
            do j = 1, comb%n
               xyz_opt(i, j, icycle) = comb%xyz(i, j)
            end do
         end do
         final_e(icycle) = etot

         write (iopt, '(i0)') comb%n
         write (iopt, '(f20.14)') etot
         do j = 1, comb%n
            write (iopt, '(a4,2x,3f20.14)') comb%sym(j), comb%xyz(:, j)*autoang
         end do
         found(7, i) = etot
      end do

      env%unit = tmp_unit
      call remove_file(itemp)
      call close_file(iopt)

      !> Include pocket structure if present
      if (pocket_grid) then
         do k = 1, kk
            final_e(n_opt + k) = pocket_e(k)
            do j = 1, comb%n
               do i = 1, 3
                  xyz_opt(i, j, n_opt + k) = xyz_pocket(i, j, k)
               end do
            end do
         end do
      end if

      !>Sorting strucutures
      call sortxyz_e(n, n_opt+kk, final_e, xyz_opt)
      call open_file(ifinal, 'final_structures.xyz', 'w')
      do i = 1, n_opt + kk !kk=0 if no pocket search
         write (ifinal, '(i0)') comb%n
         write (ifinal, '(f20.14)') final_e(i)
         do j = 1, comb%n
            write (ifinal, '(a4,2x,3f20.14)') comb%sym(j), xyz_opt(1, j, i)*autoang, &
                    &                         xyz_opt(2, j, i)*autoang, xyz_opt(3, j, i)*autoang
         end do
      end do
      call close_file(ifinal)

      call remove_file(iopt)

      ! Write best structure in format of the largest input molecule
      comb%xyz(:,:) = xyz_opt(:,:,1)
      if(molA%n >= molB%n) then
         comb%ftype=molA%ftype
      else
         comb%ftype=molB%ftype
      end if
      call generateFileName(fin_name, 'best', '', comb%ftype)
      call open_file(ifinal, fin_name, 'w')
      call writeMolecule(comb, ifinal, energy=final_e(1))
      !If not xyz then best.xyz is written to not have api break
      if(comb%ftype /= 1)then
         call open_file(ifinal, 'best.xyz', 'w')
         write (ifinal, '(i0)') comb%n
         write (ifinal, '(f20.14)') final_e(1)
         do j = 1, comb%n
            write (ifinal, '(a4,2x,3f20.14)') comb%sym(j), xyz_opt(1, j, 1)*autoang, &
            &                                 xyz_opt(2, j, 1)*autoang, xyz_opt(3, j, 1)*autoang
         end do
      end if
      call close_file(ifinal)


      call delete_file(set%opt_logfile)

      !> Printout Interaction Energy
      write (env%unit, *)
      write (env%unit, *) '  ---------------------------'
      write (env%unit, *) '     Interaction energies'
      write (env%unit, *) '  ---------------------------'
      write (env%unit, *) 'Attention: monomers are not optimized'
      write (env%unit, *) 'Interaction energies are not physical'
      write(*,*) 
      write (env%unit, '(2x,''Lowest Interaction Energy:'',1x, F10.2, 1x, ''kcal/mol'' )') &
              & (final_e(1) - molA_e - molB_e)*au

      write (env%unit, *) ' #   E_int (kcal/mol)'
      do i = 1, n_opt + kk
         E_int = (final_e(i) - molA_e - molB_e)*au
         write (*, '(2x,i0,3x,F10.2)') i, E_int
      end do

      if (pocket_grid) then
         deallocate (xyz_pocket, pocket_e)
      end if
      deallocate (xyz_opt, final_e)

   end subroutine docking_search

!Sorts xyz according to an energy (minimum first)
!Necessary as the xyzsort and xyzsort2 does somehow not work
   subroutine sortxyz_e(n, nall, e, xyz)

      INTERFACE
         RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices)
            REAL, DIMENSION(:), INTENT(INOUT) :: Item       ! array of values
            INTEGER, INTENT(IN)    :: First, Last
            INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
         END SUBROUTINE Quicksort
      END INTERFACE

      integer, intent(in) :: n, nall
      real(wp), intent(inout) :: e(nall), xyz(3, n, nall)

      integer :: ind(nall), i, j, k
      real(wp) :: e_tmp(nall), xyz_tmp(3, n, nall)

      do i = 1, nall
         ind(i) = i
      end do

      call Qsort(e, 1, nall, ind)

      xyz_tmp = xyz
      do i = 1, nall
         do k = 1, n
            do j = 1, 3
               xyz(j, k, i) = xyz_tmp(j, k, ind(i))
            end do
         end do
      end do

   end subroutine sortxyz_e

   subroutine restart_gff(env, mol, calc)

      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      type(TMolecule), intent(inout) :: mol
      type(TGFFCalculator), intent(inout)::calc

      integer :: itopo = 32

      !> Read the constrain again with new xyz only if necessary
      if (constraint_xyz) then
         call read_userdata(xcontrol, env, mol)
         call constrain_xTB_gff(env, mol)
      end if

      !if(auto_walls) call number_walls=0; read_userdata(xcontrol,env,mol) !everytime new wall pot

      call open_file(itopo, 'gfnff_topo', 'r')
      call remove_file(itopo)
      call open_file(itopo, 'charges', 'r')
      call remove_file(itopo)
      call calc%topo%zero
      calc%update = .true.
      call gfnff_param_dealloc(calc%topo)
      call newD3Model(calc%topo%dispm, mol%n, mol%at)
      call gfnff_set_param(mol%n, calc%gen, calc%param)
      if (allocated(calc%neigh%nb)) deallocate(calc%neigh%nb)
      allocate (calc%neigh%nb(calc%neigh%numnb, mol%n, calc%neigh%numctr), source=0)
      if (allocated(calc%topo%qfrag)) deallocate(calc%topo%qfrag)
      allocate( calc%topo%qfrag(mol%n), source = 0.0d0 )
      if (allocated(calc%topo%fraglist)) deallocate(calc%topo%fraglist)
      allocate( calc%topo%fraglist(mol%n), source = 0 )
      if (allocated(calc%neigh%iTrUsed)) deallocate(calc%neigh%iTrUsed)
      if (allocated(calc%neigh%bpair)) deallocate(calc%neigh%bpair)
      if (allocated(calc%neigh%blist)) deallocate(calc%neigh%blist)
      if (allocated(calc%neigh%vbond)) deallocate(calc%neigh%vbond)
      if (allocated(calc%neigh%nr_hb)) deallocate(calc%neigh%nr_hb)
      calc%topo%qfrag(1) = set%ichrg
      calc%topo%qfrag(2:mol%n) = 0.0_wp
      call gfnff_ini(env, .false., ini, mol, calc%gen,&
           &         calc%param, calc%topo, calc%neigh, set%efield, calc%accuracy)

   end subroutine restart_gff

   subroutine restart_xTB(env, mol, chk, calc, basisset)

      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      type(TMolecule), intent(inout) :: mol
      type(TRestart), intent(inout) :: chk
      type(TxTBCalculator), intent(inout)::calc
      logical, optional, intent(in) :: basisset

      real(wp), allocatable :: cn(:), dcn(:, :, :), dq(:, :, :), g(:, :)
      real(wp) :: er
      logical :: okbas
      integer :: i

      allocate (cn(mol%n), g(3, mol%n))

      !> Read the constrain again with new xyz only if necessary
      if (constraint_xyz) then
         number_walls = 0 !Reset wall potentials
         do i=1, maxwalls
!            deallocate(wpot(i)%list)
            if(allocated(wpot(i)%list)) deallocate(wpot(i)%list)
         end do
         call read_userdata(xcontrol, env, mol)
         call constrain_xTB_gff(env, mol)
      end if

      if (present(basisset)) then
         if (basisset) then
            deallocate (calc%basis)
            allocate (calc%basis)
            call newBasisset(calc%xtbData, mol%n, mol%at, calc%basis, okbas)
         end if
      end if

      calc%etemp = set%etemp
      calc%maxiter = set%maxscciter

      call chk%wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)
      ! Make sure number of electrons is initialized an multiplicity is consistent
      chk%wfn%nel = nint(sum(mol%z) - mol%chrg)
      chk%wfn%nopen = mol%uhf
      if (chk%wfn%nopen == 0 .and. mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1

      !> EN charges and CN
      if (set%gfn_method < 2) then
         call ncoord_d3(mol%n, mol%at, mol%xyz, cn)
      else
         call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
      end if
      if (mol%npbc > 0) then
         chk%wfn%q = real(set%ichrg, wp)/real(mol%n, wp)
      else
            if (set%guess_charges == p_guess_gasteiger) then
               call iniqcn(mol%n, mol%at, mol%z, mol%xyz, set%ichrg, 1.0_wp, chk%wfn%q, cn, set%gfn_method, .true.)
            else if (set%guess_charges == p_guess_goedecker) then
               call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
               call goedecker_chrgeq(mol%n, mol%at, mol%xyz, real(set%ichrg, wp), cn, dcn, chk%wfn%q, dq, er, g, &
                                     .false., .false., .false.)
            else
               call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
               chk%wfn%q = real(set%ichrg, wp) / real(mol%n, wp)
            end if
      end if
      !> initialize shell charges from gasteiger charges
      call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell,chk%wfn%q,chk%wfn%qsh,set%gfn_method)
      call delete_file('.sccnotconverged')
      call env%checkpoint("Setup for calculation failed")

   end subroutine restart_xTB

   subroutine constrain_xTB_gff(env, mol)

      type(TEnvironment), intent(inout) :: env

      type(TMolecule) :: mol
      real(wp), save :: fc
      logical, save :: called = .false.

      !> restraining potential
      if (allocated(potset%xyz)) then
      if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, potset%xyz)
      if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, potset%xyz)
      if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, potset%xyz)
         call setup_constrain_pot(mol%n, mol%at, potset%xyz)
      else
         if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, mol%xyz)
         if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, mol%xyz)
         if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, mol%xyz)
         call setup_constrain_pot(mol%n, mol%at, mol%xyz)
      end if
      if (.not. called) then
         called = .true.
         fc = potset%pos%fc
      else
         potset%pos%fc = fc
      end if

   end subroutine constrain_xTB_gff

   subroutine doubles(nlow, ndim, rthr, ethr, found)
      integer, intent(in) :: nlow, ndim
      real(wp), intent(in) :: rthr
      real(wp), intent(in) :: ethr
      real(wp), intent(inout) :: found(7, ndim)

      real(wp) :: scr(7, ndim)
      integer :: i, j, k
      real(wp) :: r, de, dem

      scr = found

      do i = 1, nlow - 1
         do j = i + 1, nlow
            call rnorm(i, j, found, r)
            de = abs(found(7, i) - found(7, j))
            dem = abs(found(7, i) + found(7, j))*0.5
            if (r .lt. rthr .and. de/dem .lt. ethr) then
               scr(7, j) = 1.d+42
            end if
            if (abs(found(7, i) - found(7, j)) .lt. 1.d-3) scr(7, j) = 1.d+42 ! chiral or symmetry case
         end do
      end do

      found = scr

   end subroutine doubles

end module xtb_docking_search_nci
