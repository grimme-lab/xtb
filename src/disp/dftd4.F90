! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2020, NVIDIA CORPORATION. All rights reserved.
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

module xtb_disp_dftd4
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_param, only: gam => chemical_hardness
   use xtb_mctc_blas, only : mctc_gemv
   !! ========================================================================
   !  mix in the covalent coordination number from the ncoord module
   !  also get the CN-Parameters to inline the CN-derivative in the gradient
   use xtb_disp_dftd4param, only: zeff, thopi, ootpi, p_mbd_none, p_mbd_rpalike, &
      &                   p_mbd_exact_atm, p_mbd_approx_atm, p_refq_goedecker, &
      &                   p_refq_gfn2xtb
   use xtb_param_sqrtzr4r2, only : sqrtZr4r2
   use xtb_type_dispersionmodel
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourList
   use xtb_type_param, only : dftd_parameter
   implicit none


   interface d4_gradient
      module procedure :: d4_full_gradient_neigh
      module procedure :: d4_full_gradient_latp
      module procedure :: d4_gradient_neigh
      module procedure :: d4_gradient_latp
   end interface d4_gradient


   interface d4_atm_gradient
      module procedure :: d4_atm_gradient_neigh
      module procedure :: d4_atm_gradient_latp
   end interface d4_atm_gradient


contains

subroutine newD3Model(dispm,nat,at)
   use xtb_disp_dftd4param
   implicit none
   type(TDispersionModel), intent(out) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)

   integer  :: i,ia,is,icn,j,ja,ii,jj,ij
   integer  :: cncount(0:18)
   real(wp) :: alpha(23),c6

   intrinsic :: nint

   call init(dispm, maxElem=maxval(at))

   dispm%atoms = 0
   dispm%nref = 0.0_wp

   !> set up ncount und alpha, also obtain the dimension of the dispmat
   do i = 1, nat
      ia = at(i)
      if (dispm%atoms(ia).eq.0) then
         dispm%nref(ia) = refn(ia)
         do j = 1, dispm%nref(ia)
            is = refsys(j,ia)
            alpha = sscale(is)*secaiw(:,is)
            dispm%cn(j,ia) = refcn(j,ia)
            dispm%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia) &
               &                - hcount(j,ia)*alpha), 0.0_wp)
         enddo
      endif
      dispm%atoms(ia) = dispm%atoms(ia)+1
   enddo

   ! integrate C6 coefficients
   do i = 1, maxval(at)
      do j = 1, i
         if (dispm%atoms(i) > 0 .and. dispm%atoms(j) > 0) then
            do ii = 1, dispm%nref(i)
               do jj = 1, dispm%nref(j)
                  alpha = dispm%alpha(:,ii,i)*dispm%alpha(:,jj,j)
                  c6 = thopi * trapzd(alpha)
                  dispm%c6(ii,jj,i,j) = c6
                  dispm%c6(jj,ii,j,i) = c6
               enddo
            enddo
         endif
      enddo
   enddo

end subroutine newD3Model

subroutine newD4Model(dispm,g_a,g_c,mode)
   use xtb_disp_dftd4param
   type(TDispersionModel), intent(out) :: dispm
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode

   integer  :: i,ia,is,icn,j,ii,jj
   integer  :: cncount(0:18)
   real(wp) :: sec_al(23),iz,c6,alpha(23)
   real(wp) :: tmp_hq(7,118)

   intrinsic :: nint

   call init(dispm)

   secq = 0.0_wp
   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
!     print'(1x,''* using PBE0/def2-TZVP Hirshfeld charges'')'
      refq = dftq
      refh = dfth
      secq = dfts
!  case(2)
!     refq = pbcq
!     refh = pbch
!     secq = pbcs
   case(p_refq_gasteiger)
!     print'(1x,''* using classical Gasteiger charges'')'
      refq = gffq
      refh = gffh
      secq = gffs
   case(p_refq_goedecker)
      refq = clsq
      refh = clsh
      secq = clss
   case(p_refq_gfn2xtb_gbsa_h2o)
!     print'(1x,''* using GFN2-xTB//GBSA(H2O) charges'')'
      refq = solq
      refh = solh
      secq = sols
   end select

   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
      dispm%q = dftq
      tmp_hq = dfth
   case(p_refq_gasteiger)
      dispm%q = gffq
      tmp_hq = gffh
   case(p_refq_goedecker)
      dispm%q = clsq
      tmp_hq = clsh
   case(p_refq_gfn2xtb_gbsa_h2o)
      dispm%q = solq
      tmp_hq = solh
   case default
      dispm%q = refq
      tmp_hq = refh
   end select

   dispm%atoms = 0
   dispm%nref = 0

   do ia = 1, 118
      cncount = 0
      cncount(0) = 1
      dispm%nref(ia) = refn(ia)
      do j = 1, refn(ia)
         is = refsys(j,ia)
         iz = zeff(is)
         sec_al = sscale(is)*secaiw(:,is) &
            &  * zeta(g_a,gam(is)*g_c,secq(is)+iz,tmp_hq(j,ia)+iz)
         dispm%cn(j,ia) = refcovcn(j,ia)
         icn =nint(refcn(j,ia))
         cncount(icn) = cncount(icn) + 1
         dispm%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia)-hcount(j,ia)*sec_al),0.0_wp)
      enddo
      do j = 1, refn(ia)
         icn = cncount(nint(refcn(j,ia)))
         dispm%ncount(j,ia) = icn*(icn+1)/2
      enddo
   end do

   ! integrate C6 coefficients
   do i = 1, 118
      do j = 1, i
         do ii = 1, dispm%nref(i)
            do jj = 1, dispm%nref(j)
               alpha = dispm%alpha(:,ii,i)*dispm%alpha(:,jj,j)
               c6 = thopi * trapzd(alpha)
               dispm%c6(jj,ii,j,i) = c6
               dispm%c6(ii,jj,i,j) = c6
            enddo
         enddo
      enddo
   enddo

end subroutine newD4Model

subroutine d4dim(dispm,nat,at,ndim)
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   integer, intent(out) :: ndim

   integer :: i

   ndim = 0

   do i = 1, nat
      ndim = ndim + dispm%nref(at(i))
   enddo

end subroutine d4dim

subroutine prmolc6(molc6,molc8,molpol,nat,at,  &
      & cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_io, only : stdout
   use xtb_mctc_symbols, only : toSymbol
   real(wp),intent(in)  :: molc6
   real(wp),intent(in)  :: molc8
   real(wp),intent(in)  :: molpol
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in),optional :: cn(nat)
   real(wp),intent(in),optional :: covcn(nat)
   real(wp),intent(in),optional :: q(nat)
   real(wp),intent(in),optional :: qlmom(3,nat)
   real(wp),intent(in),optional :: c6ab(nat,nat)
   real(wp),intent(in),optional :: alpha(nat)
   real(wp),intent(in),optional :: rvdw(nat)
   real(wp),intent(in),optional :: hvol(nat)
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(stdout,'(a)')
   write(stdout,'(''   #   Z   '')',advance='no')
   if(present(cn))   write(stdout,'(''        CN'')',advance='no')
   if(present(covcn))write(stdout,'(''     covCN'')',advance='no')
   if(present(q))    write(stdout,'(''         q'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(stdout,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(stdout,'(''      C6AA'')',advance='no')
   if(present(alpha))write(stdout,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(stdout,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(stdout,'(''    relVol'')',advance='no')
   write(*,'(a)')
   do i=1,nat
      write(*,'(i4,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),toSymbol(at(i))
      if(present(cn))   write(stdout,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(stdout,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(stdout,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(stdout, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(stdout,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(stdout,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(stdout,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(stdout,'(f10.3)',advance='no')hvol(i)
      write(*,'(a)')
   enddo
   endif
   write(stdout,'(/,1x,''Mol. C6AA /au·bohr⁶  :'',f18.6,'// &
   &         '/,1x,''Mol. C8AA /au·bohr⁸  :'',f18.6,'// &
   &         '/,1x,''Mol. α(0) /au        :'',f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

subroutine mdisp(dispm,nat,ndim,at,q,xyz,g_a,g_c, &
      &     gw,c6abns,molc6,molc8,molpol,aout,cout,rout,vout)
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: molc6
   real(wp),intent(out) :: molc8
   real(wp),intent(out) :: molpol
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: cout(nat,nat)
   real(wp),intent(out),optional :: rout(nat)
   real(wp),intent(out),optional :: vout(nat)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: qmod,oth,iz
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: rvdw(:)
   real(wp),allocatable :: phv(:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   parameter (oth=1._wp/3._wp)
   
   allocate( zetvec(ndim),rvdw(nat),phv(nat),c6ab(nat,nat),aw(23,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   molc6  = 0._wp
   molc8  = 0._wp
   molpol = 0._wp

   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = k+1
         itbl(ii,i) = k
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz)
         aw(:,i) = aw(:,i) + zetvec(k) * dispm%alpha(:,ii,ia)
      enddo
!     van-der-Waals radius, alpha = 4/3 pi r**3 <=> r = (3/(4pi) alpha)**(1/3)
      rvdw(i) = (0.25_wp*thopi*aw(1,i))**oth
!     pseudo-hirshfeld volume
      phv(i) = aw(1,i)/dispm%alpha(1,1,ia)
      c6ab(i,i) = thopi * trapzd(aw(:,i)**2)
      molpol = molpol + aw(1,i)
      molc6  = molc6  + c6ab(i,i)
      molc8 = molc8 + 3*sqrtZr4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*sqrtZr4r2(ia)*sqrtZr4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp

pure elemental function zeta(a,c,qref,qmod)
   !$acc routine seq
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: zeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1._wp - exp( c * ( 1._wp - qref/qmod ) ) ) )
   endif

end function zeta

pure elemental function dzeta(a,c,qref,qmod)
   !$acc routine seq
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: dzeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      dzeta = 0._wp
   else
      dzeta = - a * c * exp( c * ( 1._wp - qref/qmod ) ) &
      &           * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

pure function trapzd(pol)
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp)            :: tmp1, tmp2
   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)
!  just precalculate all weights and get the job done
   real(wp),parameter :: weights(23) = 0.5_wp * (/ &
&  ( freq (2) - freq (1) ),  &
&  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
&  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
&  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
&  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
&  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
&  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
&  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
&  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
&  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
&  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
&  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
&  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
&  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
&  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
&  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
&  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
&  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
&  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
&  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
&  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
&  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
&  ( freq(23) - freq(22) ) /)

!!  do average between trap(1)-trap(22) .and. trap(2)-trap(23)
!   tmp1 = 0.5_wp * ( &
!&  ( freq (2) - freq (1) ) * ( pol (2) + pol (1) )+ &
!&  ( freq (4) - freq (3) ) * ( pol (4) + pol (3) )+ &
!&  ( freq (6) - freq (5) ) * ( pol (6) + pol (5) )+ &
!&  ( freq (8) - freq (7) ) * ( pol (8) + pol (7) )+ &
!&  ( freq(10) - freq (9) ) * ( pol(10) + pol (9) )+ &
!&  ( freq(12) - freq(11) ) * ( pol(12) + pol(11) )+ &
!&  ( freq(14) - freq(13) ) * ( pol(14) + pol(13) )+ &
!&  ( freq(16) - freq(15) ) * ( pol(16) + pol(15) )+ &
!&  ( freq(18) - freq(17) ) * ( pol(18) + pol(17) )+ &
!&  ( freq(20) - freq(19) ) * ( pol(20) + pol(19) )+ &
!&  ( freq(22) - freq(21) ) * ( pol(22) + pol(21) ))
!   tmp2 = 0.5_wp * ( &
!&  ( freq (3) - freq (2) ) * ( pol (3) + pol (2) )+ &
!&  ( freq (5) - freq (4) ) * ( pol (5) + pol (4) )+ &
!&  ( freq (7) - freq (6) ) * ( pol (7) + pol (6) )+ &
!&  ( freq (9) - freq (8) ) * ( pol (9) + pol (8) )+ &
!&  ( freq(11) - freq(10) ) * ( pol(11) + pol(10) )+ &
!&  ( freq(13) - freq(12) ) * ( pol(13) + pol(12) )+ &
!&  ( freq(15) - freq(14) ) * ( pol(15) + pol(14) )+ &
!&  ( freq(17) - freq(16) ) * ( pol(17) + pol(16) )+ &
!&  ( freq(19) - freq(18) ) * ( pol(19) + pol(18) )+ &
!&  ( freq(21) - freq(20) ) * ( pol(21) + pol(20) )+ &
!&  ( freq(23) - freq(22) ) * ( pol(23) + pol(22) ))

   trapzd = sum(pol*weights)

end function trapzd

pure elemental function cngw(wf,cn,cnref)
   !$acc routine seq
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight
   real(wp)            :: val

   intrinsic :: exp

   val = -wf * ( cn - cnref )**2
   if (val < -200.0_wp) then ! technically, exp(-200) -> 1.383897e-87
     cngw = 0.0_wp
   else
     cngw = exp ( val )
   end if

end function cngw

pure elemental function dcngw(wf,cn,cnref)
   real(wp),intent(in) :: wf,cn,cnref
   real(wp) :: dcngw

   dcngw = 2*wf*(cnref-cn)*cngw(wf,cn,cnref)

end function dcngw

!* BJ damping function ala DFT-D3(BJ)
!  f(n,rab) = sn*rab**n/(rab**n + R0**n)  w/ R0 = a1*sqrt(C6/C8)+a2
!  see: https://doi.org/10.1002/jcc.21759
pure elemental function fdmpr_bj(n,r,c) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: fdmp
   fdmp = 1.0_wp / ( r**n + c**n )
end function fdmpr_bj
pure elemental function fdmprdr_bj(n,r,c) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: dfdmp
   dfdmp = -n*r**(n-1) * fdmpr_bj(n,r,c)**2
end function fdmprdr_bj

!* original DFT-D3(0) damping
!  f(n,rab) = sn/(1+6*(4/3*R0/rab)**alp)  w/ R0 of unknown origin
pure elemental function fdmpr_zero(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (c/r)**(n+alp)))
end function fdmpr_zero
pure elemental function fdmprdr_zero(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(c/r)**(alp)) &
             - alp*r**n/c*six*(c/r)**(alp-1) ) &
           * fdmpr_zero(n,r,c,alp)**2
!  fdmp = 1.0_wp / (r**n*(1 + 6.0_wp * (c/r)**(n+alp)))
end function fdmprdr_zero

!* fermi damping function from TS and MBD methods
!  f(n,rab) = sn/(1+exp[-alp*(rab/R0-1)]) w/ R0 as experimenal vdW-Radii
pure elemental function fdmpr_fermi(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1.0_wp+exp(-alp*(r/c - 1.0))))
end function fdmpr_fermi
pure elemental function fdmprdr_fermi(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = -(-alp/c*r**n*exp(-alp*(r/c - 1.0)) &
             + n*r**(n-1)*(1.0_wp+exp(-alp*(r/c - 1.0)))) &
             * fdmpr_fermi(n,r,c,alp)**2
end function fdmprdr_fermi

!* optimized power zero damping (M. Head-Gordon)
!  f(n,rab) = sn*rab**(n+alp)/(rab**(n+alp) + R0**(n+alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.7b00176
pure elemental function fdmpr_op(n,r,c,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmpr_op
pure elemental function fdmprdr_op(n,r,c,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = (alp*r*(alp-1) - (n+alp)*r**alp*r**(n+alp-1)) &
           * fdmpr_op(n,r,c,alp)**2
!  fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmprdr_op

!* Sherrill's M-zero damping function
!  f(n,rab) = sn/(1+6*(4/3*R0/rab+a2*R0)**(-alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.6b00780
pure elemental function fdmpr_zerom(n,r,c,rsn,alp) result(fdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (r/c+rsn*c)**(-alp)))
end function fdmpr_zerom
pure elemental function fdmprdr_zerom(n,r,c,rsn,alp) result(dfdmp)
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(r/c+rsn*c)**(-alp)) &
              - alp*r**n/c*six*(r/c+rsn*c)**(-alp-1) ) &
           * fdmpr_zerom(n,r,c,rsn,alp)**2
end function fdmprdr_zerom


subroutine d4(dispm,nat,ndim,at,wf,g_a,g_c,covcn,gw,c6abns)
   use xtb_mctc_accuracy, only : wp
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: covcn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: c6abns(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   c6abns = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, dispm%nref(ia)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),dispm%cn(ii,ia))
         enddo
      enddo
      if (norm > 1e-80_wp) then
         norm = 1._wp / norm
      else
         norm = 0.0_wp
      end if
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),dispm%cn(ii,ia)) * norm
         enddo
         if (norm == 0.0_wp) then
            if (abs(maxval(dispm%cn(:dispm%nref(ia),ia)) &
               & - dispm%cn(ii,ia)) < 1e-12_wp) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         do j = 1, i-1
            ja = at(j)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6abns(l,k) = dispm%c6(ii,jj,ia,ja)
               c6abns(k,l) = c6abns(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine d4


subroutine build_wdispmat(dispm,nat,ndim,at,itbl,xyz,par,c6abns,gw,wdispmat)
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(:)
   integer, intent(in)  :: itbl(:)
   real(wp),intent(in)  :: xyz(:, :)
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: c6abns(:, :)
   real(wp),intent(in)  :: gw(:)
   real(wp),intent(out) :: wdispmat(:, :)

   integer :: i,ii,ia,j,jj,ja,k,l
   real(wp) :: c8abns,c10abns,r2,cutoff,oor6,oor8,oor10,r,gwgw,r4r2ij
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient
   real(wp), parameter :: gwcut = 1.0e-7_wp

   ! acc enter data create(wdispmat) copyin(at, xyz, itbl, dispm, dispm%nref, &
   ! acc& c6abns, gw, par)
 
   ! acc kernels default(present)
   wdispmat = 0.0_wp
   ! acc end kernels

!#ifdef XTB_GPU
   ! acc parallel default(present)
   ! acc loop gang collapse(2) private(r4r2ij, cutoff, r2, oor6, oor8, oor10)
!#else
   !$omp parallel do shared(wdispmat) &
   !$omp shared(nat, at, xyz, itbl, par, dispm, c6abns, gw) &
   !$omp private(ia, k, j, ja, l, r4r2ij, cutoff, r2, oor6, oor8, oor10, &
   !$omp& ii, jj, gwgw, c8abns, c10abns)
!#endif
   do i = 1, nat
      do j = 1, nat
         if (j >= i) cycle
         ia = at(i)
         ja = at(j)
         k = itbl(i)
         l = itbl(j)
         r4r2ij = 3.0_wp*sqrtZr4r2(ia)*sqrtZr4r2(ja)
         cutoff = par%a1*sqrt(r4r2ij)+par%a2
         r2 = sum( (xyz(:,j)-xyz(:,i))**2 )
         if (r2.gt.rthr*rthr) cycle
         oor6  = 1.0_wp/(r2**3 + cutoff**6 )
         oor8  = 1.0_wp/(r2**4 + cutoff**8 )
         oor10 = 1.0_wp/(r2**5 + cutoff**10)
         ! acc loop seq
         do ii = 1, dispm%nref(ia)
            ! acc loop seq private(gwgw, c8abns, c10abns)
            do jj = 1, dispm%nref(ja)
               gwgw = gw(k+ii)*gw(l+jj)
               if (gwgw.lt.gwcut) cycle
               c8abns  = r4r2ij * c6abns(k+ii,l+jj)
               c10abns = 49.0_wp/40.0_wp * r4r2ij**2 * c6abns(k+ii,l+jj)
               wdispmat(k+ii,l+jj) = gw(k+ii)*gw(l+jj) * ( &
               &  - par%s6  * ( c6abns(k+ii,l+jj)  * oor6 ) &
               &  - par%s8  * ( c8abns       * oor8 ) &
               &  - par%s10 * ( c10abns      * oor10) )
               wdispmat(l+jj,k+ii) = wdispmat(k+ii,l+jj)
            enddo
         enddo
      enddo
   enddo
!#ifdef XTB_GPU
   ! acc end parallel

   ! acc exit data copyout(wdispmat) delete(itbl, sqrtZr4r2, dispm, dispm%nref, &
   ! acc& c6abns, gw)
!#endif

end subroutine build_wdispmat


subroutine disppot(dispm,nat,ndim,at,itbl,q,g_a,g_c,wdispmat,gw,hdisp)
   use xtb_mctc_blas, only : mctc_symv
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(:)
   integer, intent(in)  :: itbl(:)
   real(wp),intent(in)  :: q(:)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(:,:)
   real(wp),intent(in)  :: gw(:)
   real(wp),intent(inout) :: hdisp(:)

   integer  :: iat,ii,k,ati
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),zerovec(ndim),dumvec(ndim), source = 0._wp )

   !$acc enter data create(zetavec, zerovec, dumvec) copyin(hdisp, itbl, at, zeff, &
   !$acc& dispm, dispm%nref, dispm%q, q, gw, wdispmat, g_a, g_c)

   !$acc kernels default(present)
   zetavec = 0.0_wp
   zerovec = 0.0_wp
   dumvec  = 0.0_wp
   !$acc end kernels

#ifdef XTB_GPU
   !$acc parallel default(present)
   !$acc loop gang private(k, ati, iz)
#else
   !$omp parallel do shared(itbl, at, dispm, gw, q, g_a, g_c) &
   !$omp private(k, ati, iz, ii)
#endif
   do iat = 1, nat
      k = itbl(iat)
      ati = at(iat)
      iz = zeff(ati)
      !$acc loop vector
      do ii = 1, dispm%nref(ati)
         if (gw(k+ii).lt.gw_cut) cycle
         zerovec(k+ii) = dzeta(g_a,gam(ati)*g_c,dispm%q(ii,ati)+iz,q(iat)+iz)
         zetavec(k+ii) =  zeta(g_a,gam(ati)*g_c,dispm%q(ii,ati)+iz,q(iat)+iz)
      enddo
   enddo
#ifdef XTB_GPU
   !$acc end parallel
   !$acc exit data copyout(zetavec)
#endif

!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call mctc_symv(wdispmat,zetavec,dumvec)

!  get atomic reference contributions
#ifdef XTB_GPU
   !$acc parallel default(present)
   !$acc loop gang private(k, ati)
#else
   !$omp parallel do reduction(+:hdisp) shared(itbl, at, dumvec, zerovec) &
   !$omp private(ati, ii, k)
#endif
   do iat = 1, nat
      k = itbl(iat)
      ati = at(iat)
      !$acc loop vector
      do ii = 1, dispm%nref(ati)
         !$acc atomic
         hdisp(iat) = hdisp(iat) + dumvec(k+ii)*zerovec(k+ii)
      end do
   enddo
#ifdef XTB_GPU
   !$acc end parallel

   !$acc exit data copyout(hdisp) delete(zerovec, dumvec, itbl, at, zeff, &
   !$acc& dispm, dispm%nref, dispm%q, q, g_a, g_c, wdispmat, gw)
#endif

end subroutine disppot


function edisp_scc(dispm,nat,ndim,at,itbl,q,g_a,g_c,wdispmat,gw) result(ed)
   use xtb_mctc_blas, only : mctc_symv, mctc_dot
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(:)
   integer, intent(in)  :: itbl(:)
   real(wp),intent(in)  :: q(:)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(:,:)
   real(wp),intent(in)  :: gw(:)
   real(wp) :: ed

   integer  :: iat,ii,k,ati
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),dumvec(ndim))

   zetavec = 0.0_wp
   dumvec = 0.0_wp
   ed = 0.0_wp

   !$omp parallel do shared(zetavec, dispm, nat, itbl, at, q, g_c, g_a) &
   !$omp private(k, ati, iz, ii)
   do iat = 1, nat
       k = itbl(iat)
       ati = at(iat)
       iz = zeff(ati)
       do ii = 1, dispm%nref(ati)
          if (gw(k+ii).lt.gw_cut) cycle
          zetavec(k+ii) = zeta(g_a,gam(ati)*g_c,dispm%q(ii,ati)+iz,q(iat)+iz)
      enddo
   enddo

!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim)
   call mctc_symv(wdispmat,zetavec,dumvec,alpha=0.5_wp)
   ed = mctc_dot(dumvec,zetavec)

end function edisp_scc


! --- PBC
subroutine pbc_d4(dispm,nat,ndim,at,wf,g_a,g_c,covcn,gw,refc6)
   type(TDispersionModel), intent(in) :: dispm
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: covcn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: refc6(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   refc6 = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, dispm%nref(ia)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),dispm%cn(ii,ia))
         enddo
      enddo
      if (norm > 1e-80_wp) then
         norm = 1._wp / norm
      else
         norm = 0.0_wp
      end if
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),dispm%cn(ii,ia)) * norm
         enddo
         if (norm == 0.0_wp) then
            if (abs(maxval(dispm%cn(:dispm%nref(ia),ia)) &
               & - dispm%cn(ii,ia)) < 1e-12_wp) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         ! diagonal terms
         do jj = 1, ii
            l = itbl(jj,i)
            refc6(l,k) = dispm%c6(ii,jj,ia,ia)
            refc6(k,l) = refc6(l,k)
         enddo
         ! offdiagonal terms
         do j = 1, i-1
            ja = at(j)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               refc6(l,k) = dispm%c6(ii,jj,ia,ja)
               refc6(k,l) = refc6(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine pbc_d4


!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(dispm, nat, atoms, g_a, g_c, wf, q, cn, zeff, gam, &
      &                      zetavec, zerovec, zetadcn, zerodcn, zetadq)
   type(TDispersionModel), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Charge scaling height.
   real(wp), intent(in) :: g_a
   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)
   real(wp), intent(in) :: zeff(:)
   real(wp), intent(in) :: gam(:)
   !> weighting and scaling function for the atomic reference systems
   real(wp), intent(out) :: zetaVec(:, :)
   !> weighting and scaling function for the atomic reference systems for q=0
   real(wp), intent(out) :: zeroVec(:, :)
   !> derivative of the weight'n'scale function w.r.t. the partial charges
   real(wp), intent(out) :: zetadq(:, :)
   !> derivative of the weight'n'scale function w.r.t. the coordination number
   real(wp), intent(out) :: zetadcn(:, :)
   !> derivative of the weight'n'scale function w.r.t. the CN for q=0
   real(wp), intent(out) :: zerodcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, twf, gw, expw, expd, gwk, dgwk
   real(wp) :: gi, zi

   ! acc enter data create(zetavec, zerovec, zetadq, zetadq, zetadcn, zerodcn) &
   ! acc& copyin(dispm, dispm%nref, dispm%ncount, dispm%cn, dispm%q, atoms, &
   ! acc& cn, q, zeff, gam, g_a, g_c, wf)

   ! acc kernels default(present)
   zetavec = 0.0_wp
   zerovec = 0.0_wp
   zetadcn = 0.0_wp
   zerodcn = 0.0_wp
   zetadq  = 0.0_wp
   ! acc end kernels

!#ifdef XTB_GPU
   ! acc parallel default(present)
   ! acc loop gang private(zi, gi, norm, dnorm)
!#else
   !$omp parallel do shared(zetavec, zetadcn, zetadq, zerodcn) &
   !$omp shared(nat, atoms, dispm, cn, q, g_a, g_c, wf, zerovec) &
   !$omp private(iat, ati, zi, gi, norm, dnorm, iref, icount, twf, gw, expw, &
   !$omp& expd, gwk, dgwk)
!#endif
   do iat = 1, nat
      ati = atoms(iat)

      zi = zeff(ati)
      gi = g_c * gam(ati)

      norm = 0.0_wp
      dnorm = 0.0_wp
      ! acc loop vector
      do iref = 1, dispm%nref(ati)
         ! acc loop seq private(twf, gw)
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            norm = norm + gw
            dnorm = dnorm + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo
      end do
      if (norm > 1e-80_wp) then
         norm = 1.0_wp / norm
      else
         norm = 0.0_wp
      end if
      ! acc loop vector private(expw, expd)
      do iref = 1, dispm%nref(ati)
         expw = 0.0_wp
         expd = 0.0_wp
         ! acc loop seq private(twf, gw)
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            expw = expw + gw
            expd = expd + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo

         gwk = expw * norm
         if (norm == 0.0_wp) then
            if (abs(maxval(dispm%cn(:dispm%nref(ati), ati)) &
               & - dispm%cn(iref, ati)) < 1e-12_wp) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zetavec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerovec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         zetadcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * dgwk
         zetadq(iref, iat) = dzeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerodcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * dgwk

      end do
   end do
!#ifdef XTB_GPU
   ! acc end parallel

   ! acc exit data copyout(zetavec, zerovec, zetadq, zetadq, zetadcn, zerodcn) &
   ! acc& delete(dispm, dispm%nref, dispm%ncount, dispm%cn, dispm%q, atoms, &
   ! acc& cn, q, zeff, gam, g_a, g_c, wf)
!#endif

end subroutine weight_references

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(dispm, nat, atoms, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)
   type(TDispersionModel), intent(in) :: dispm
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting and scaling function for the atomic reference systems
   real(wp), intent(in) :: zetaVec(:, :)
   !> derivative of the weight'n'scale function w.r.t. the partial charges
   real(wp), intent(in) :: zetadq(:, :)
   !> derivative of the weight'n'scale function w.r.t. the coordination number
   real(wp), intent(in) :: zetadcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)
   !> derivative of the C6 w.r.t. the partial charge
   real(wp), intent(out) :: dc6dq(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj

   !$acc enter data create(c6, dc6dcn, dc6dq) copyin(atoms, dispm, dispm%nref, dispm%c6, &
   !$acc& zetavec, zetadcn, zetadq)

   !$acc kernels default(present)
   c6 = 0.0_wp
   dc6dcn = 0.0_wp
   dc6dq = 0.0_wp
   !$acc end kernels

#ifdef XTB_GPU
   !$acc parallel default(present)
   !$acc loop gang collapse(2)
#else
   !$omp parallel do default(none) shared(c6, dc6dcn, dc6dq) &
   !$omp shared(nat, atoms, dispm, zetavec, zetadcn, zetadq) &
   !$omp private(iat, ati, jat, atj, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj, &
   !$omp& iref, jref, refc6)
#endif
   do iat = 1, nat
      do jat = 1, nat
         if (jat > iat) cycle
         ati = atoms(iat)
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         dc6dqi = 0.0_wp
         dc6dqj = 0.0_wp
         !$acc loop vector collapse(2)
         do iref = 1, dispm%nref(ati)
            do jref = 1, dispm%nref(atj)
               refc6 = dispm%c6(iref, jref, ati, atj)
               dc6 = dc6 + zetavec(iref, iat) * zetavec(jref, jat) * refc6
               dc6dcni = dc6dcni + zetadcn(iref, iat) * zetavec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + zetavec(iref, iat) * zetadcn(jref, jat) * refc6
               dc6dqi = dc6dqi + zetadq(iref, iat) * zetavec(jref, jat) * refc6
               dc6dqj = dc6dqj + zetavec(iref, iat) * zetadq(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
         dc6dq(iat, jat) = dc6dqi
         dc6dq(jat, iat) = dc6dqj
      end do
   end do
#ifdef XTB_GPU
   !$acc end parallel

   !$acc exit data copyout(c6, dc6dcn, dc6dq) delete(atoms, dispm, dispm%nref, dispm%c6, &
   !$acc& zetavec, zetadcn, zetadq)
#endif

end subroutine get_atomic_c6


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_full_gradient_neigh &
      & (mol, dispm, neighs, neighs3, neighlist, par, g_a, g_c, wf, &
      &  cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma, e2, e3)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Static neighbourlist.
   type(TNeighbourList), intent(in) :: neighlist

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)

   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs3(:)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)

   !> Derivative of partial charges w.r.t. atomic coordinates.
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   real(wp), intent(out), optional :: e2
   real(wp), intent(out), optional :: e3

   integer :: nat, max_ref

   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:), dEdq(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), energies3(nat), dEdcn(nat), dEdq(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)

   call disp_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, &
      & c6, dc6dcn, dc6dq, energies, gradient, sigma, dEdcn, dEdq)

   if (present(e2)) e2 = sum(energies)
   if (par%s9 /= 0.0_wp) then
      call get_atomic_c6(dispm, nat, mol%at, zerovec, zerodcn, zerodq, &
         & c6, dc6dcn, dc6dq)
      call atm_gradient_neigh(mol, neighs3, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
         & energies3, gradient, sigma, dEdcn)
   end if
   if (present(e3)) e3 = sum(energies3)

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (present(dqdr)) then
      call mctc_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
   end if
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   if (present(dqdL)) then
      call mctc_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   energy = energy + sum(energies) + sum(energies3)

end subroutine d4_full_gradient_neigh


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_gradient_neigh &
      & (mol, dispm, neighs, neighlist, par, g_a, g_c, wf, &
      &  cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Static neighbourlist.
   type(TNeighbourList), intent(in) :: neighlist

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)

   !> Derivative of partial charges w.r.t. atomic coordinates.
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref

   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:), dEdq(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), dEdcn(nat), dEdq(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)

   call disp_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, &
      & c6, dc6dcn, dc6dq, energies, gradient, sigma, dEdcn, dEdq)

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (present(dqdr)) then
      call mctc_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
   end if
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   if (present(dqdL)) then
      call mctc_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   energy = energy + sum(energies)

end subroutine d4_gradient_neigh


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine disp_gradient_neigh &
      & (mol, neighs, neighlist, par, r4r2, c6, dc6dcn, dc6dq, &
      &  energies, gradient, sigma, dEdcn, dEdq)

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Static neighbourlist.
   type(TNeighbourList), intent(in) :: neighlist

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: r4r2(:)

   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)

   real(wp), intent(in) :: c6(:, :)

   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(in) :: dc6dq(:, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energies(:)

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   real(wp), intent(inout) :: dEdcn(:)

   real(wp), intent(inout) :: dEdq(:)

   integer :: iat, jat, ati, atj, ij, img

   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn, dEdq) &
   !$omp shared(mol, neighs, neighlist, par, r4r2, c6, dc6dcn, dc6dq) &
   !$omp private(ij, img, jat, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dist2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)

         r4r2ij = 3*r4r2(ati)*r4r2(atj)
         r0 = par%a1*sqrt(r4r2ij) + par%a2

         t6 = 1._wp/(r2**3+r0**6)
         t8 = 1._wp/(r2**4+r0**8)
         t10 = 1._wp/(r2**5+r0**10)

         d6 = -6*r2**2*t6**2
         d8 = -8*r2**3*t8**2
         d10 = -10*r2**4*t10**2

         disp = par%s6*t6 + par%s8*r4r2ij*t8 &
            &  + par%s10*49.0_wp/40.0_wp*r4r2ij**2*t10
         ddisp= par%s6*d6 + par%s8*r4r2ij*d8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*d10

         dE = -c6(iat, jat)*disp * 0.5_wp
         dG = -c6(iat, jat)*ddisp*rij
         dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

         energies(iat) = energies(iat) + dE
         dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
         dEdq(iat) = dEdq(iat) - dc6dq(iat, jat) * disp
         sigma = sigma + dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + dE
            dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
            dEdq(jat) = dEdq(jat) - dc6dq(jat, iat) * disp
            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG
            sigma = sigma + dS
         endif

      enddo
   enddo
   !$omp end parallel do

end subroutine disp_gradient_neigh


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_atm_gradient_neigh &
      & (mol, dispm, neighs, neighlist, par, g_a, g_c, wf, cn, dcndr, dcndL, &
      &  energy, gradient, sigma)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Static neighbourlist.
   type(TNeighbourList), intent(in) :: neighlist

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref

   real(wp), allocatable :: q(:)
   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     q(nat), c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zerovec, zerodcn, zerodq, &
      & c6, dc6dcn, dc6dq)

   call atm_gradient_neigh &
      & (mol, neighs, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d4_atm_gradient_neigh


subroutine atm_gradient_neigh &
      & (mol, neighs, neighlist, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, kat, ati, atj, atk, jtr, ktr, ij, jk, ik
   real(wp) :: rij(3), rjk(3), rik(3), r2ij, r2jk, r2ik
   real(wp) :: c6ij, c6jk, c6ik, cij, cjk, cik, scale
   real(wp) :: dE, dG(3, 3), dS(3, 3), dCN(3)
   real(wp), parameter :: sr = 4.0_wp/3.0_wp

   !$omp parallel do default(none) reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(mol, neighs, neighlist, par, r4r2, c6, dc6dcn) &
   !$omp private(iat, ati, ij, jtr, r2ij, rij, jat, atj, ik, ktr, kat, atk, rik, &
   !$omp& r2ik, rjk, r2jk, c6ij, cij, c6ik, c6jk, cik, cjk, scale, dE, dG, dS, dCN)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         jtr = neighlist%ineigh(ij, iat)
         r2ij = neighlist%dist2(ij, iat)
         rij = neighlist%coords(:, jtr) - neighlist%coords(:, iat)
         jat = neighlist%image(jtr)
         atj = mol%at(jat)

         c6ij = c6(jat,iat)
         cij = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atj))+par%a2

         do ik = 1, ij-1
            ktr = neighlist%ineigh(ik, iat)
            rik = neighlist%coords(:, ktr) - neighlist%coords(:, iat)
            r2ik = neighlist%dist2(ik, iat)
            rjk = neighlist%coords(:, ktr) - neighlist%coords(:, jtr)
            r2jk = sum(rjk**2)
            kat = neighlist%image(ktr)
            atk = mol%at(kat)

            c6ik = c6(kat,iat)
            c6jk = c6(kat,jat)

            cik = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atk))+par%a2
            cjk = par%a1*sqrt(3.0_wp*r4r2(atj)*r4r2(atk))+par%a2

            call deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
               & r2ij, r2jk, r2ik, dc6dcn(iat,jat), dc6dcn(jat,iat), &
               & dc6dcn(jat,kat), dc6dcn(kat,jat), dc6dcn(iat,kat), &
               & dc6dcn(kat,iat), rij, rjk, rik, par%alp, dE, dG, dS, dCN)

            scale = par%s9 * triple_scale(iat, jat, kat)
            energies(iat) = energies(iat) + dE * scale/3
            energies(jat) = energies(jat) + dE * scale/3
            energies(kat) = energies(kat) + dE * scale/3
            gradient(:, iat) = gradient(:, iat) + dG(:, 1) * scale
            gradient(:, jat) = gradient(:, jat) + dG(:, 2) * scale
            gradient(:, kat) = gradient(:, kat) + dG(:, 3) * scale
            sigma(:, :) = sigma + dS * scale
            dEdcn(iat) = dEdcn(iat) + dCN(1) * scale
            dEdcn(jat) = dEdcn(jat) + dCN(2) * scale
            dEdcn(kat) = dEdcn(kat) + dCN(3) * scale

         end do
      end do
   end do
   !$omp end parallel do

end subroutine atm_gradient_neigh


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_full_gradient_latp &
      & (mol, dispm, trans, par, g_a, g_c, wf, cutoff, cutoff3, &
      &  cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma, e2, e3)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   !> Cutoff for threebody interactions
   real(wp), intent(in) :: cutoff3

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)

   !> Derivative of partial charges w.r.t. atomic coordinates.
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   real(wp), intent(out), optional :: e2
   real(wp), intent(out), optional :: e3

   integer :: nat, max_ref

   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:), dEdq(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), energies3(nat), dEdcn(nat), dEdq(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)

   call disp_gradient_latp(mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, dc6dq, &
      &  energies, gradient, sigma, dEdcn, dEdq)

   if (present(e2)) e2 = sum(energies)
   if (par%s9 /= 0.0_wp) then
      call get_atomic_c6(dispm, nat, mol%at, zerovec, zerodcn, zerodq, &
         & c6, dc6dcn, dc6dq)
#ifdef XTB_GPU
      call atm_gradient_latp_gpu(mol, trans, cutoff3, par, sqrtZr4r2, c6, dc6dcn, &
         & energies3, gradient, sigma, dEdcn)
#else
      call atm_gradient_latp(mol, trans, cutoff3, par, sqrtZr4r2, c6, dc6dcn, &
         & energies3, gradient, sigma, dEdcn)
#endif
   end if
   if (present(e3)) e3 = sum(energies3)

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (present(dqdr)) then
      call mctc_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
   end if
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   if (present(dqdL)) then
      call mctc_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   energy = energy + sum(energies) + sum(energies3)

end subroutine d4_full_gradient_latp


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_gradient_latp &
      & (mol, dispm, trans, par, g_a, g_c, wf, cutoff, &
      &  cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Partial charge of every atom.
   real(wp), intent(in) :: q(:)

   !> Derivative of partial charges w.r.t. atomic coordinates.
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref

   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:), dEdq(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), dEdcn(nat), dEdq(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zetavec, zetadcn, zetadq, &
      & c6, dc6dcn, dc6dq)

   call disp_gradient_latp(mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, dc6dq, &
      & energies, gradient, sigma, dEdcn, dEdq)

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (present(dqdr)) then
      call mctc_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
   end if
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   if (present(dqdL)) then
      call mctc_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   energy = energy + sum(energies)

end subroutine d4_gradient_latp


subroutine disp_gradient_latp &
      & (mol, trans, cutoff, par, r4r2, c6, dc6dcn, dc6dq, &
      &  energies, gradient, sigma, dEdcn, dEdq)

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: r4r2(:)

   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   real(wp), intent(in) :: c6(:, :)

   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(in) :: dc6dq(:, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energies(:)

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   real(wp), intent(inout) :: dEdcn(:)

   real(wp), intent(inout) :: dEdq(:)

   integer :: iat, jat, ati, atj, itr

   real(wp) :: cutoff2
   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   cutoff2 = cutoff**2
   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn, dEdq) &
   !$omp shared(mol, trans, cutoff2, par, r4r2, c6, dc6dcn, dc6dq) &
   !$omp private(iat, jat, itr, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do jat = 1, iat
         atj = mol%at(jat)

         r4r2ij = 3*r4r2(ati)*r4r2(atj)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-10_wp) cycle

            t6 = 1._wp/(r2**3+r0**6)
            t8 = 1._wp/(r2**4+r0**8)
            t10 = 1._wp/(r2**5+r0**10)

            d6 = -6*r2**2*t6**2
            d8 = -8*r2**3*t8**2
            d10 = -10*r2**4*t10**2

            disp = par%s6*t6 + par%s8*r4r2ij*t8 &
               &  + par%s10*49.0_wp/40.0_wp*r4r2ij**2*t10
            ddisp= par%s6*d6 + par%s8*r4r2ij*d8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*d10

            dE = -c6(iat, jat)*disp * 0.5_wp
            dG = -c6(iat, jat)*ddisp*rij
            dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

            energies(iat) = energies(iat) + dE
            dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
            dEdq(iat) = dEdq(iat) - dc6dq(iat, jat) * disp
            sigma = sigma + dS
            if (iat /= jat) then
               energies(jat) = energies(jat) + dE
               dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
               dEdq(jat) = dEdq(jat) - dc6dq(jat, iat) * disp
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma = sigma + dS
            end if

         end do
      end do
   end do
   !$omp end parallel do

end subroutine disp_gradient_latp


!> Evaluate gradient of DFT-D4, this routine can handle systems of arbitrary
!  periodicity due to the static neighbourlist.
subroutine d4_atm_gradient_latp &
      & (mol, dispm, trans, par, g_a, g_c, wf, cutoff, cn, dcndr, dcndL, &
      &  energy, gradient, sigma)
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   type(TDispersionModel), intent(in) :: dispm

   !> Molecular Structure information.
   type(TMolecule), intent(in) :: mol

   !> Becke--Johnson damping parameters.
   type(dftd_parameter), intent(in) :: par

   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Charge scaling height.
   real(wp), intent(in) :: g_a

   !> Charge scaling steepness.
   real(wp), intent(in) :: g_c

   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf

   !> Cutoff for pairwise interactions
   real(wp), intent(in) :: cutoff

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN w.r.t. atomic coordinates.
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of CN w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Dispersion energy.
   real(wp), intent(inout) :: energy

   !> Derivative of the dispersion energy w.r.t. atomic positions.
   real(wp), intent(inout) :: gradient(:, :)

   !> Stress tensor resulting from dispersion interactions.
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref

   real(wp), allocatable :: q(:)
   real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
   real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:), dEdq(:)

   nat = len(mol)
   max_ref = maxval(dispm%nref(mol%at))
   allocate(zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &     zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &     q(nat), c6(nat, nat), dc6dcn(nat, nat), dc6dq(nat, nat), &
      &     energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(dispm, nat, mol%at, g_a, g_c, wf, q, cn, zeff, gam, &
      & zetavec, zerovec, zetadcn, zerodcn, zetadq)

   call get_atomic_c6(dispm, nat, mol%at, zerovec, zerodcn, zerodq, &
      & c6, dc6dcn, dc6dq)

#ifdef XTB_GPU
   call atm_gradient_latp_gpu &
      & (mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)
#else
   call atm_gradient_latp &
      & (mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)
#endif

   call mctc_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   call mctc_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d4_atm_gradient_latp

subroutine atm_gradient_latp &
      & (mol, trans, cutoff, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, kat, ati, atj, atk, jtr, ktr
   real(wp) :: cutoff2
   real(wp) :: rij(3), rjk(3), rik(3), r2ij, r2jk, r2ik
   real(wp) :: c6ij, c6jk, c6ik, cij, cjk, cik, scale
   real(wp) :: dE, dG(3, 3), dS(3, 3), dCN(3)
   real(wp), parameter :: sr = 4.0_wp/3.0_wp

   cutoff2 = cutoff**2

   !$omp parallel do default(none) reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(mol, r4r2, par, trans, cutoff2, c6, dc6dcn) &
   !$omp private(iat, ati, jat, atj, kat, atk, c6ij, cij, c6ik, c6jk, cik, cjk, &
   !$omp& rij, r2ij, ktr, rik, r2ik, rjk, r2jk, scale, dE, dG, dS, dCN)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do jat = 1, iat
         atj = mol%at(jat)

         c6ij = c6(jat,iat)
         cij = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atj))+par%a2

         do kat = 1, jat
            atk = mol%at(kat)

            c6ik = c6(kat,iat)
            c6jk = c6(kat,jat)

            cik = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atk))+par%a2
            cjk = par%a1*sqrt(3.0_wp*r4r2(atj)*r4r2(atk))+par%a2

            do jtr = 1, size(trans, dim=2)
               rij = mol%xyz(:, jat) - mol%xyz(:, iat) + trans(:, jtr)
               r2ij = sum(rij**2)
               if (r2ij > cutoff2 .or. r2ij < 1.0e-14_wp) cycle
               do ktr = 1, size(trans, dim=2)
                  if (jat == kat .and. jtr == ktr) cycle
                  rik = mol%xyz(:, kat) - mol%xyz(:, iat) + trans(:, ktr)
                  r2ik = sum(rik**2)
                  if (r2ik > cutoff2 .or. r2ik < 1.0e-14_wp) cycle
                  rjk = mol%xyz(:, kat) - mol%xyz(:, jat) + trans(:, ktr) &
                     & - trans(:, jtr)
                  r2jk = sum(rjk**2)
                  if (r2jk > cutoff2 .or. r2jk < 1.0e-14_wp) cycle

                  call deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
                     & r2ij, r2jk, r2ik, dc6dcn(iat,jat), dc6dcn(jat,iat), &
                     & dc6dcn(jat,kat), dc6dcn(kat,jat), dc6dcn(iat,kat), &
                     & dc6dcn(kat,iat), rij, rjk, rik, par%alp, dE, dG, dS, dCN)

                  scale = par%s9 * triple_scale(iat, jat, kat)
                  energies(iat) = energies(iat) + dE * scale/3
                  energies(jat) = energies(jat) + dE * scale/3
                  energies(kat) = energies(kat) + dE * scale/3
                  gradient(:, iat) = gradient(:, iat) + dG(:, 1) * scale
                  gradient(:, jat) = gradient(:, jat) + dG(:, 2) * scale
                  gradient(:, kat) = gradient(:, kat) + dG(:, 3) * scale
                  sigma(:, :) = sigma + dS * scale
                  dEdcn(iat) = dEdcn(iat) + dCN(1) * scale
                  dEdcn(jat) = dEdcn(jat) + dCN(2) * scale
                  dEdcn(kat) = dEdcn(kat) + dCN(3) * scale

               end do
            end do

         end do
      end do
   end do
   !$omp end parallel do

end subroutine atm_gradient_latp

subroutine atm_gradient_latp_gpu &
      & (mol, trans, cutoff, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, kat, ati, atj, atk, jtr, ktr
   real(wp) :: cutoff2
   real(wp) :: rij(3), rjk(3), rik(3), r2ij, r2jk, r2ik
   real(wp) :: c6ij, c6jk, c6ik, cij, cjk, cik, scale
   real(wp) :: dE, dG(3, 3), dS(3, 3), dCN(3)
   real(wp), parameter :: sr = 4.0_wp/3.0_wp
   integer :: mlen, k, kk

   real(wp) :: c9, dc9, ccc1, rrr1, rrr2, rrr3, ang, dang, fdmp, dfdmp, dGr, cr

   cutoff2 = cutoff**2
   mlen = len(mol)

   !$acc enter data copyin(par,trans,r4r2,c6,dc6dcn,energies,gradient,sigma,dEdcn, &
   !$acc& mol,mol%at,mol%xyz)

   !$acc parallel default(present) private(rij,rjk,rik,dG,dS,dCN) 

   !$acc loop gang collapse(3)
   do iat = 1, mlen
      do jat = 1, mlen
         do kat = 1, mlen
            if (jat.gt.iat) cycle
            if (kat.gt.jat) cycle

            ati = mol%at(iat)
            atj = mol%at(jat)

            c6ij = c6(jat,iat)
            cij = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atj))+par%a2

            atk = mol%at(kat)

            c6ik = c6(kat,iat)
            c6jk = c6(kat,jat)

            cik = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atk))+par%a2
            cjk = par%a1*sqrt(3.0_wp*r4r2(atj)*r4r2(atk))+par%a2

            do jtr = 1, size(trans, dim=2)
               rij = mol%xyz(:, jat) - mol%xyz(:, iat) + trans(:, jtr)
               r2ij = sum(rij**2)
               if (r2ij > cutoff2 .or. r2ij < 1.0e-14_wp) cycle
               do ktr = 1, size(trans, dim=2)
                  if (jat == kat .and. jtr == ktr) cycle
                  rik = mol%xyz(:, kat) - mol%xyz(:, iat) + trans(:, ktr)
                  r2ik = sum(rik**2)
                  if (r2ik > cutoff2 .or. r2ik < 1.0e-14_wp) cycle
                  rjk = mol%xyz(:, kat) - mol%xyz(:, jat) + trans(:, ktr) &
                     & - trans(:, jtr)
                  r2jk = sum(rjk**2)
                  if (r2jk > cutoff2 .or. r2jk < 1.0e-14_wp) cycle

                  call deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
                     & r2ij, r2jk, r2ik, dc6dcn(iat,jat), dc6dcn(jat,iat), &
                     & dc6dcn(jat,kat), dc6dcn(kat,jat), dc6dcn(iat,kat), &
                     & dc6dcn(kat,iat), rij, rjk, rik, par%alp, dE, dG, dS, dCN)

   !c9 = -sqrt(c6ij*c6ik*c6jk)

   !ccc1 = cij*cjk*cik

   !rrr2 = r2ij*r2jk*r2ik
   !rrr1 = sqrt(rrr2)
   !rrr3 = rrr1*rrr2

   !ang = 0.375_wp * (r2ij+r2jk-r2ik)*(r2ij-r2jk+r2ik)*(-r2ij+r2jk+r2ik) &
   !   & / (rrr3*rrr2) + 1.0_wp/(rrr3)

   !cr = (ccc1/rrr1)**(1.0_wp/3.0_wp)
   !fdmp = 1.0_wp/(1.0_wp + 6.0_wp*cr**par%alp)
   !dfdmp = -(2.0_wp*par%alp*cr**par%alp) * fdmp**2

   !! Energy contribution
   !dE = -fdmp*ang*c9

   !! Derivative w.r.t. i-j distance
   !dang = -0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
   !   & +r2ij*(3.0_wp*r2jk**2+2.0_wp*r2jk*r2ik+3.0_wp*r2ik**2) &
   !   & -5.0_wp*(r2jk-r2ik)**2*(r2jk+r2ik)) / (rrr3*rrr2)
   !dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ij
   !dG(:, 1) = -dGr * rij
   !dG(:, 2) = +dGr * rij 
   !dS(:, :) = 0.5_wp * dGr * spread(rij, 1, 3) * spread(rij, 2, 3)

   !! Derivative w.r.t. i-k distance
   !dang = -0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
   !   & +r2ik*(3.0_wp*r2jk**2+2.0*r2jk*r2ij+3.0_wp*r2ij**2) &
   !   & -5.0_wp*(r2jk-r2ij)**2*(r2jk+r2ij)) / (rrr3*rrr2)
   !dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ik
   !dG(:, 1) = -dGr * rik + dG(:, 1)
   !dG(:, 3) = +dGr * rik 
   !dS(:, :) = 0.5_wp * dGr * spread(rik, 1, 3) * spread(rik, 2, 3) + dS

   !! Derivative w.r.t. j-k distance
   !dang=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
   !   & +r2jk*(3.0_wp*r2ik**2+2.0_wp*r2ik*r2ij+3.0_wp*r2ij**2) &
   !   & -5.0_wp*(r2ik-r2ij)**2*(r2ik+r2ij)) / (rrr3*rrr2)
   !dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2jk
   !dG(:, 2) = -dGr * rjk + dG(:, 2)
   !dG(:, 3) = +dGr * rjk + dG(:, 3)
   !dS(:, :) = 0.5_wp * dGr * spread(rjk, 1, 3) * spread(rjk, 2, 3) + dS

   !! CN derivative
   !dc9 = 0.5_wp*c9*(dc6dcn(iat,jat)/c6ij+dc6dcn(iat,kat)/c6ik)
   !dCN(1) = -ang*fdmp*dc9
   !dc9 = 0.5_wp*c9*(dc6dcn(jat,iat)/c6ij+dc6dcn(jat,kat)/c6jk)
   !dCN(2) = -ang*fdmp*dc9
   !dc9 = 0.5_wp*c9*(dc6dcn(kat,iat)/c6ik+dc6dcn(kat,jat)/c6jk)
   !dCN(3) = -ang*fdmp*dc9

                  scale = par%s9 * triple_scale(iat, jat, kat)
                  !$acc atomic
                  energies(iat) = energies(iat) + dE * scale/3
                  !$acc atomic
                  energies(jat) = energies(jat) + dE * scale/3
                  !$acc atomic
                  energies(kat) = energies(kat) + dE * scale/3
                  do k = 1,3
                    !$acc atomic
                    gradient(k, iat) = gradient(k, iat) + dG(k, 1) * scale
                    !$acc atomic
                    gradient(k, jat) = gradient(k, jat) + dG(k, 2) * scale
                    !$acc atomic
                    gradient(k, kat) = gradient(k, kat) + dG(k, 3) * scale
                  enddo
                  do k = 1,3
                    do kk = 1,3
                      !$acc atomic
                      sigma(kk, k) = sigma(kk, k) + dS(kk, k) * scale
                    enddo
                  enddo
                  !$acc atomic
                  dEdcn(iat) = dEdcn(iat) + dCN(1) * scale
                  !$acc atomic
                  dEdcn(jat) = dEdcn(jat) + dCN(2) * scale
                  !$acc atomic
                  dEdcn(kat) = dEdcn(kat) + dCN(3) * scale

               end do
            end do

         end do
      end do
   end do
   !$acc end parallel

   !$acc exit data copyout(energies,gradient,sigma,dEdcn)
   !$acc exit data delete(par,trans,r4r2,c6,dc6dcn,mol,mol%at,mol%xyz)


end subroutine atm_gradient_latp_gpu


pure subroutine deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
      & r2ij, r2jk, r2ik, dc6ij, dc6ji, dc6jk, dc6kj, dc6ik, dc6ki, &
      & rij, rjk, rik, alp, dE, dG, dS, dCN)

   !$acc routine vector
   real(wp), intent(in) :: c6ij, c6ik, c6jk
   real(wp), intent(in) :: cij, cjk, cik
   real(wp), intent(in) :: r2ij, r2jk, r2ik
   real(wp), intent(in) :: dc6ij, dc6ji, dc6jk, dc6kj, dc6ik, dc6ki
   real(wp), intent(in) :: rij(3), rjk(3), rik(3)
   integer, intent(in) :: alp
   real(wp), intent(out) :: dE, dG(3, 3), dS(3, 3), dCN(3)

   real(wp) :: c9, dc9, ccc1, rrr1, rrr2, rrr3, ang, dang, fdmp, dfdmp, dGr, cr

   c9 = -sqrt(c6ij*c6ik*c6jk)

   ccc1 = cij*cjk*cik

   rrr2 = r2ij*r2jk*r2ik
   rrr1 = sqrt(rrr2)
   rrr3 = rrr1*rrr2

   ang = 0.375_wp * (r2ij+r2jk-r2ik)*(r2ij-r2jk+r2ik)*(-r2ij+r2jk+r2ik) &
      & / (rrr3*rrr2) + 1.0_wp/(rrr3)

   cr = (ccc1/rrr1)**(1.0_wp/3.0_wp)
   fdmp = 1.0_wp/(1.0_wp + 6.0_wp*cr**alp)
   dfdmp = -(2.0_wp*alp*cr**alp) * fdmp**2

   ! Energy contribution
   dE = -fdmp*ang*c9

   ! Derivative w.r.t. i-j distance
   dang = -0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
      & +r2ij*(3.0_wp*r2jk**2+2.0_wp*r2jk*r2ik+3.0_wp*r2ik**2) &
      & -5.0_wp*(r2jk-r2ik)**2*(r2jk+r2ik)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ij
   dG(:, 1) = -dGr * rij
   dG(:, 2) = +dGr * rij 
   dS(:, :) = 0.5_wp * dGr * spread(rij, 1, 3) * spread(rij, 2, 3)

   ! Derivative w.r.t. i-k distance
   dang = -0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
      & +r2ik*(3.0_wp*r2jk**2+2.0*r2jk*r2ij+3.0_wp*r2ij**2) &
      & -5.0_wp*(r2jk-r2ij)**2*(r2jk+r2ij)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ik
   dG(:, 1) = -dGr * rik + dG(:, 1)
   dG(:, 3) = +dGr * rik 
   dS(:, :) = 0.5_wp * dGr * spread(rik, 1, 3) * spread(rik, 2, 3) + dS

   ! Derivative w.r.t. j-k distance
   dang=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
      & +r2jk*(3.0_wp*r2ik**2+2.0_wp*r2ik*r2ij+3.0_wp*r2ij**2) &
      & -5.0_wp*(r2ik-r2ij)**2*(r2ik+r2ij)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2jk
   dG(:, 2) = -dGr * rjk + dG(:, 2)
   dG(:, 3) = +dGr * rjk + dG(:, 3)
   dS(:, :) = 0.5_wp * dGr * spread(rjk, 1, 3) * spread(rjk, 2, 3) + dS

   ! CN derivative
   dc9 = 0.5_wp*c9*(dc6ij/c6ij+dc6ik/c6ik)
   dCN(1) = -ang*fdmp*dc9
   dc9 = 0.5_wp*c9*(dc6ji/c6ij+dc6jk/c6jk)
   dCN(2) = -ang*fdmp*dc9
   dc9 = 0.5_wp*c9*(dc6ki/c6ik+dc6kj/c6jk)
   dCN(3) = -ang*fdmp*dc9

end subroutine deriv_atm_triple


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(scale)
   !$acc routine seq

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: scale

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         scale = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         scale = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         scale = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         scale = 0.5_wp
      end if
   end if

end function triple_scale


end module xtb_disp_dftd4
