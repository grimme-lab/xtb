! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module tbmod_dftd4
   use iso_fortran_env, only : wp => real64
   use mctc_constants, only : pi
   use mctc_param, only: gam => chemical_hardness, &
      &                  en => pauling_en, &
      &                  rcov => covalent_radius_d3, &
      &                  r4r2 => sqrt_z_r4_over_r2
   !! ========================================================================
   !  mix in the covalent coordination number from the ncoord module
   !  also get the CN-Parameters to inline the CN-derivative in the gradient
   use ncoord, only : covncoord => ncoord_d4, kn,k1,k4,k5,k6
   use tbdef_param, only : dftd_parameter
   use tbpar_dftd4, only: zeff, thopi, ootpi, p_mbd_none, p_mbd_rpalike, &
      &                   p_mbd_exact_atm, p_mbd_approx_atm, p_refq_goedecker, &
      &                   p_refq_gfn2xtb
   use tbdef_dispersion_model
   implicit none

   type(tb_dispersion_model) :: dispm

contains

subroutine d4init(g_a,g_c,mode)
   use tbpar_dftd4
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode

   integer  :: i,ia,is,icn,j,ii,jj
   integer  :: cncount(0:15)
   real(wp) :: sec_al(23),iz,c6,alpha(23)
   real(wp) :: tmp_hq(7,118)

   intrinsic :: nint

   dispm = tb_dispersion_model()

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

end subroutine d4init

subroutine d4dim(nat,at,ndim)
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
           &       cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   use iso_fortran_env, only : id => output_unit
   use mctc_econv, only : autoaa
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
   character(len=2),   external :: asym
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(id,'(a)')
   write(id,'(''   #   Z   '')',advance='no')
   if(present(cn))   write(id,'(''        CN'')',advance='no')
   if(present(covcn))write(id,'(''     covCN'')',advance='no')
   if(present(q))    write(id,'(''         q'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(id,'(''      C6AA'')',advance='no')
   if(present(alpha))write(id,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(id,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(id,'(''    relVol'')',advance='no')
   write(*,'(a)')
   do i=1,nat
      write(*,'(i4,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),asym(at(i))
      if(present(cn))   write(id,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(id,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(id,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(id,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(id,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(id,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(id,'(f10.3)',advance='no')hvol(i)
      write(*,'(a)')
   enddo
   endif
   write(id,'(/,1x,''Mol. C6AA /au·bohr⁶  :'',f18.6,'// &
   &         '/,1x,''Mol. C8AA /au·bohr⁸  :'',f18.6,'// &
   &         '/,1x,''Mol. α(0) /au        :'',f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

subroutine mdisp(nat,ndim,at,q,xyz,g_a,g_c, &
           &     gw,c6abns,molc6,molc8,molpol,aout,cout,rout,vout)
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
      molc8 = molc8 + 3*r4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*r4r2(ia)*r4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp

pure elemental function zeta(a,c,qref,qmod)
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
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight

   intrinsic :: exp

   cngw = exp ( -wf * ( cn - cnref )**2 )

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

subroutine d4(nat,ndim,at,wf,g_a,g_c,covcn,gw,c6abns)
   use iso_fortran_env, only : wp => real64
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
      norm = 1._wp / norm
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),dispm%cn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical precision, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(dispm%cn(:dispm%nref(ia),ia)).eq.dispm%cn(ii,ia)) then
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

pure subroutine build_dispmat(nat,ndim,at,xyz,par,c6abns,dispmat)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: dispmat(ndim,ndim)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns,c10abns,r,r2,oor6,oor8,oor10,cutoff
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient

   allocate( itbl(7,nat), source = 0 )

   dispmat = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         cutoff = par%a1*sqrt(3._wp*r4r2(at(i))*r4r2(at(j)))+par%a2
         r = norm2(xyz(:,j)-xyz(:,i))
         if (r.gt.rthr) cycle
!        oor6  = 1.0_wp/(r**6  + cutoff**6 )
!        oor8  = 1.0_wp/(r**8  + cutoff**8 )
!        oor10 = 1.0_wp/(r**10 + cutoff**10)
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c8abns = 3.0_wp * r4r2(ia)*r4r2(ja) * c6abns(k,l)
               c10abns = 49.0_wp/40.0_wp * c8abns**2/c6abns(k,l)
               dispmat(k,l) = &
               &  - par%s6 * ( c6abns(k,l) * oor6 ) &
               &  - par%s8 * ( c8abns      * oor8 ) &
               &  - par%s8 * ( c10abns     * oor8 )
               dispmat(l,k) = dispmat(k,l)
            enddo
         enddo
      enddo
   enddo

end subroutine build_dispmat

subroutine build_wdispmat(nat,ndim,at,xyz,par,c6abns,gw,wdispmat)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: wdispmat(ndim,ndim)

   integer :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns,c10abns,r2,cutoff,oor6,oor8,oor10,r,gwgw,r4r2ij
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient
   real(wp), parameter :: gwcut = 1.0e-7_wp

   allocate( itbl(7,nat), source = 0 )
 
   wdispmat = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         r4r2ij = 3.0_wp*r4r2(ia)*r4r2(ja)
         cutoff = par%a1*sqrt(r4r2ij)+par%a2
!        r2 = sum( (xyz(:,j)-xyz(:,i))**2 )
!        oor6  = 1.0_wp/(r2**3 + cutoff**6 )
!        oor8  = 1.0_wp/(r2**4 + cutoff**8 )
!        oor10 = 1.0_wp/(r2**5 + cutoff**10)
         r = norm2(xyz(:,j)-xyz(:,i))
         if (r.gt.rthr) cycle
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               gwgw = gw(k)*gw(l)
               if (gwgw.lt.gwcut) cycle
               c8abns  = r4r2ij * c6abns(k,l)
               c10abns = 49.0_wp/40.0_wp * r4r2ij**2 * c6abns(k,l)
               wdispmat(k,l) = gw(k)*gw(l) * ( &
               &  - par%s6  * ( c6abns(k,l)  * oor6 ) &
               &  - par%s8  * ( c8abns       * oor8 ) &
               &  - par%s10 * ( c10abns      * oor10) )
               wdispmat(l,k) = wdispmat(k,l)
            enddo
         enddo
      enddo
   enddo

end subroutine build_wdispmat

subroutine disppot(nat,ndim,at,q,g_a,g_c,wdispmat,gw,hdisp)
   use mctc_la, only : symv
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: hdisp(nat)

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),zerovec(ndim),dumvec(ndim), source = 0._wp )

   zetavec = 0.0_wp
   zerovec = 0.0_wp
   dumvec  = 0.0_wp
   hdisp   = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, dispm%nref(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zerovec(k) = dzeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz)
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call symv('U',ndim,1._wp,wdispmat,ndim,zetavec,1,0._wp,dumvec,1)
!  call gemv('N',ndim,ndim,1._wp,wdispmat,ndim,zetavec, &
!  &     1,0._wp,dumvec,1)
!  get atomic reference contributions
   k = 0
   do i = 1, nat
      ia = at(i)
      hdisp(i) = sum(dumvec(k+1:k+dispm%nref(ia))*zerovec(k+1:k+dispm%nref(ia)))
      k = k + dispm%nref(ia)
   enddo

   deallocate(zetavec,zerovec,dumvec)

end subroutine disppot

function edisp_scc(nat,ndim,at,q,g_a,g_c,wdispmat,gw) result(ed)
   use mctc_la, only : symv
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp) :: ed

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),dumvec(ndim), source = 0._wp )

   ed = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, dispm%nref(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call symv('U',ndim,0.5_wp,wdispmat,ndim,zetavec,1,0.0_wp,dumvec,1)
!  call gemv('N',ndim,ndim,0.5_wp,wdispmat,ndim,zetavec, &
!  &           1,0.0_wp,dumvec,1)
   ed = dot_product(dumvec,zetavec)

   deallocate(zetavec,dumvec)

end function edisp_scc

subroutine edisp(nat,ndim,at,q,xyz,par,g_a,g_c, &
           &     gw,c6abns,mbd,E,aout,etwo,emany)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: E
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: emany

   integer  :: i,ii,ia,k,ij,l,j,jj,ja
   integer, allocatable :: itbl(:,:)
   real(wp) :: Embd,qmod,c6ij,c6ij_ns,oor6,oor8,oor10,r2,cutoff,iz,r
   real(wp),allocatable :: dispmat(:,:)
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)
   real(wp),allocatable :: c6ab(:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: oor6ab(:,:)

   intrinsic :: present,sqrt,sum
   
   allocate( zetvec(ndim),aw(23,nat),oor6ab(nat,nat), &
   &         zerovec(ndim),c6ab(nat*(nat+1)/2), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   e = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz)
         zerovec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,iz)
         aw(:,i) = aw(:,i) + zerovec(k) * dispm%alpha(:,ii,ia)
      enddo
   enddo

!$OMP parallel private(i,j,ia,ja,ij,k,l,r,oor6,oor8,oor10,cutoff,c6ij,c6ij_ns) &
!$omp&         shared(c6ab) reduction(+:E)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         r = norm2(xyz(:,i)-xyz(:,j))
!        r2 = sum( (xyz(:,i)-xyz(:,j))**2 )
         cutoff = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
!        oor6 = 1._wp/(r2**3 + cutoff**6)
         oor6 = fdmpr_bj( 6,r,cutoff)
         oor6ab(i,j) = oor6
         oor6ab(j,i) = oor6
!        oor8  = 1._wp/(r2**4 + cutoff**8)
!        oor10 = 1._wp/(r2**5 + cutoff**10)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         c6ij_ns = 0.0_wp
         c6ij = 0.0_wp
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij_ns = c6ij_ns + zerovec(k)*zerovec(l)*c6abns(k,l)
               c6ij = c6ij + zetvec(k)*zetvec(l)*c6abns(k,l)
            enddo
         enddo
         c6ab(ij) = c6ij_ns
         E = E - c6ij*(par%s6*oor6 + par%s8*3._wp*r4r2(ia)*r4r2(ja)*oor8 &
         &      + par%s10*49.0_wp/40._wp*(3.0_wp*r4r2(ia)*r4r2(ja))**2*oor10 )
      enddo
   enddo
!$omp enddo
!$omp end parallel

   if (present(Etwo)) Etwo = E

   select case(mbd)
   case(p_mbd_rpalike) ! full RPA-like MBD
!     print'(1x,''* MBD effects calculated by RPA like scheme'')'
      call dispmb(Embd,aw,xyz,oor6ab,nat)
      Embd = par%s9*Embd
      E = E + Embd
   case(p_mbd_exact_atm) ! Axilrod-Teller-Muto three-body term
!     print'(1x,''* MBD effects calculated by ATM formula'')'
      call dispabc(nat,at,xyz,aw,par,Embd)
      E = E + Embd
   case(p_mbd_approx_atm) ! D3-like approximated ATM term
!     print'(1x,''* MBD effects approximated by ATM formula'')'
      call apprabc(nat,at,xyz,c6ab,par,Embd)
      E = E + Embd
   case default
      Embd = 0.0_wp
   end select

   if (present(Emany)) Emany = Embd

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, dispm%nref(ia)
            aout(:,i) = aout(:,i) + zetvec(k) * dispm%alpha(:,ii,ia)
         enddo
      enddo
   endif

end subroutine edisp

!* compute D4 gradient
!  ∂E/∂rij = ∂/∂rij (W·D·W)
!          = ∂W/∂rij·D·W  + W·∂D/∂rij·W + W·D·∂W/∂rij
!  ∂W/∂rij = ∂(ζ·w)/∂rij = ζ·∂w/∂rij = ζ·∂w/∂CN·∂CN/∂rij
!  ∂ζ/∂rij = 0
subroutine dispgrad(nat,ndim,at,q,xyz, &
           &        par,wf,g_a,g_c, &
           &        c6abns,mbd, &
           &        g,eout,dq,aout)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
!  real(wp),intent(in)  :: cn(nat) ! calculate on-the-fly
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf,g_a,g_c
!  real(wp),intent(in)  :: gw(ndim) ! calculate on-the-fly
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout
   real(wp),intent(in), optional :: dq(3,nat,nat+1)
   real(wp),intent(out),optional :: aout(23,nat)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10
   real(wp) :: c8abns,disp,x1,x2,x3
   real(wp) :: c6ij,dic6ij,djc6ij,dizij,djzij
   real(wp) :: rcovij,expterm,den,dcndr
   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dr(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:)
   real(wp),allocatable :: dzvec(:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: dgw(:)
   real(wp),allocatable :: dc6dq(:)
   real(wp),allocatable :: dzdq(:)
   real(wp) :: cn_thr,r_thr,gw_thr
   parameter(cn_thr = 1600.0_wp)
   parameter( r_thr=5000._wp)
   parameter(gw_thr=0.000001_wp)
   real(wp),parameter :: sqrtpi = 1.77245385091_wp
   real(wp),parameter :: hlfosqrtpi = 0.5_wp/1.77245385091_wp
!  timing
!  real(wp) :: time0,time1
!  real(wp) :: wall0,wall1

   intrinsic :: present,sqrt,sum,maxval,exp,abs

!  print'(" * Allocating local memory")'
   allocate( dc6dr(nat*(nat+1)/2),dc6dcn(nat),  &
   &         r2ab(nat*(nat+1)/2),cn(nat),  &
   &         zvec(ndim),dzvec(ndim),  &
   &         gw(ndim),dgw(ndim),dc6dq(nat),dzdq(ndim),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   ed = 0.0_wp

!  derivatives of the covalent coordination number
!  ∂CN/∂rij = ∂/∂rij k4·exp(-(∂EN+k5)²/(2·k6²))/(1+exp(-k1(-1)))
!  print'(" * Calculating CN")'
   call covncoord(nat,at,xyz,cn,cn_thr)

!  precalc
!  print'(" * Setting up index table")'
   k = 0
   do i = 1, nat
      do ii = 1, dispm%nref(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k,norm,dnorm,twf,tgw,dexpw,expw,gwk,dgwk)  &
!$omp shared (nat,at,dispm,itbl,wf,cn,g_a,g_c,q) &
!$omp shared (gw,dgw,zvec,dzvec,dzdq)
!$omp do
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      norm  = 0.0_wp
      dnorm = 0.0_wp
      do ii=1,dispm%nref(ia)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            tgw = cngw(twf,cn(i),dispm%cn(ii,ia))
            norm  =  norm + tgw
            dnorm = dnorm + 2*twf*(dispm%cn(ii,ia)-cn(i))*tgw
         enddo
      enddo
      norm = 1._wp/norm
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         dexpw=0.0_wp
         expw=0.0_wp
         do iii = 1, dispm%ncount(ii,ia)
            twf = wf*iii
            tgw = cngw(twf,cn(i),dispm%cn(ii,ia))
            expw  =  expw + tgw
            dexpw = dexpw + 2*twf*(dispm%cn(ii,ia)-cn(i))*tgw
         enddo

         ! save
         gwk = expw*norm
         if (gwk.ne.gwk) then
            if (maxval(dispm%cn(:dispm%nref(ia),ia)).eq.dispm%cn(ii,ia)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zvec(k) = zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         gw(k) =  zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,iz) * gwk

         dgwk = dexpw*norm-expw*dnorm*norm**2
         if (dgwk.ne.dgwk) then
            dgwk = 0.0_wp
         endif
         dzvec(k) = zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz) * dgwk
         dzdq(k) = dzeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         dgw(k) = zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,iz) * dgwk
      enddo
   enddo
!$omp end do
!$omp end parallel

!  print'(" * Entering second OMP section")'
!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ij,dic6ij,djc6ij,disp,dizij,djzij,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10,door6,door8,door10)  &
!$omp shared(nat,at,xyz,dispm,itbl,zvec,dzvec,c6abns,par,dzdq) &
!$omp shared(r2ab) reduction(+:dc6dr,dc6dq,dc6dcn,ed)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r2ab(ij) = r2
         if (r2.gt.r_thr) cycle
         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         dizij = 0.0_wp
         djzij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
               dizij = dizij + dzdq(k)*zvec(l)*c6abns(k,l)
               djzij = djzij + zvec(k)*dzdq(l)*c6abns(k,l)
            enddo
         enddo

         r = sqrt(r2)

         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2

         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
!        oor6   = fdmpr_bj( 6,r,r0)
!        oor8   = fdmpr_bj( 8,r,r0)
!        oor10  = fdmpr_bj(10,r,r0)
!        door6  = fdmprdr_bj( 6,r,r0)
!        door8  = fdmprdr_bj( 8,r,r0)
!        door10 = fdmprdr_bj(10,r,r0)

         disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
         &    + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10

         ! save
         dc6dq(i) = dc6dq(i) + dizij*disp
         dc6dq(j) = dc6dq(j) + djzij*disp
         dc6dcn(i) = dc6dcn(i) + dic6ij*disp
         dc6dcn(j) = dc6dcn(j) + djc6ij*disp
         dc6dr(ij) = dc6dr(ij) + c6ij*(par%s6*door6 + par%s8*r4r2ij*door8 &
         &                       + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10 )

         ed = ed - c6ij*disp
      enddo
   enddo
!$omp enddo
!$omp end parallel

!  select case(mbd)
!  case(1) ! full RPA-like MBD
!     print'(1x,''* MBD effects calculated by RPA like scheme'')'
!     call raise('W','MBD gradient not fully implemented yet')
!     call mbdgrad(nat,xyz,aw,daw,oor6ab,g,embd)
!  case(1,2) ! Axilrod-Teller-Muto three-body term
!     if(mbd.eq.1) then
!        call raise('W','MBD gradient not fully implemented yet')
!        print'(''MBD gradient not fully implemented yet'')'
!        print'(1x,''* calculate MBD effects with ATM formula instead'')'
!     else
!        print'(1x,''* MBD effects calculated by ATM formula'')'
!     endif
!     call dabcgrad(nat,ndim,at,xyz,par,dcn,gw,dgw,itbl,g,embd)
!  case(3) ! D3-like approximated ATM term
!     print'(1x,''* MBD effects approximated by ATM formula'')'
!  print'(" * Starting MBD gradient calculation")'
   if (mbd.ne.p_mbd_none) &
   &   call dabcappr(nat,ndim,at,xyz,par,  &
           &        r2ab,gw,dgw,c6abns,itbl,dc6dr,dc6dcn,eabc)
!  end select

!  print'(" * Entering third OMP section")'
!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,rij,r2,r,drdx,den,rcovij,  &
!$omp         expterm,dcndr,dtmp) reduction(+:g) &
!$omp shared(nat,at,xyz,dc6dr,dc6dcn)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r = sqrt(r2)
         drdx = rij/r

         if (r2.lt.cn_thr) then
            den=k4*exp(-(abs((en(ia)-en(ja)))+ k5)**2/k6 )
            rcovij=rcov(ia)+rcov(ja)
            !expterm=exp(-k1*(rcovij/r-1._wp))
            !dcndr=(-k1*rcovij*expterm*den)/(r2*((1._wp + expterm)**2))
            dcndr = -den*kn/sqrtpi/rcovij*exp(-kn**2*(r-rcovij)**2/rcovij**2)
         else
            dcndr=0.0_wp
         endif

         dtmp = dc6dr(ij) + (dc6dcn(i)+dc6dcn(j))*dcndr

         g(:,i) = g(:,i) + dtmp * drdx
         g(:,j) = g(:,j) - dtmp * drdx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   if (present(dq)) then
      call dgemv('n',3*nat,nat,-1.0_wp,dq,3*nat,dc6dq,1,1.0_wp,g,1)
   endif

!  print*,ed,eabc

!  print'(" * Dispersion all done, saving variables")'
   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, dispm%nref(ia)
            aout(:,i) = aout(:,i) + zvec(k) * dispm%alpha(:,ii,ia)
         enddo
      enddo
   endif


end subroutine dispgrad

subroutine apprabc(nat,at,xyz,c6ab,par,E)
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: c6ab(nat*(nat+1)/2)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka,ij,ik,jk
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine apprabc

subroutine dispabc(nat,at,xyz,aw,par,E)
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2)
         do k = 1, j-1
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2)
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = (par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = par%s9*thopi*trapzd( aw(:,i)*aw(:,j)*aw(:,k) )
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**par%alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine dispabc

subroutine abcappr(nat,ndim,at,xyz,g_a,g_c,par,gw,r2ab,c6abns,eabc)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: g_a,g_c
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   integer, allocatable :: itbl(:,:)
   real(wp),allocatable :: c6ab(:),zvec(:),c(:)
   real(wp) :: r2ij,r2jk,r2ik,iz
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: sqrt

   allocate( c6ab(nat*(nat+1)/2),zvec(ndim),c(nat*(nat+1)/2),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   eabc = 0.0_wp

!  precalc
   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = k+1
         itbl(ii,i) = k
         ! NEW: q=0 for ATM
         zvec(k) = zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,iz) * gw(k)
      enddo
   enddo

!$OMP parallel private(i,ia,j,ja,ij,r2ij,c6ij)  &
!$omp&         shared (c6ab,c)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         c(ij) = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel private(i,j,ij,ia,ja,k,ka,ik,jk,atm,fdmp,  &
!$omp&                 r2ij,cij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp&                 c9ijk,oor9ijk) &
!$omp&         reduction(+:eabc)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = c(ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = c(ik)
            cjk   = c(jk)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)

            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)

            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   deallocate( c6ab,c,zvec )

end subroutine abcappr

subroutine dabcappr(nat,ndim,at,xyz,par,  &
                &        r2ab,zvec,dzvec,c6abns,itbl,dc6dr,dc6dcn,eout)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: dc6dr(nat*(nat+1)/2)
   real(wp),intent(inout)        :: dc6dcn(nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(nat*(nat+1)/2),dc6ab(nat,nat),  &
   &         source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij,k,l)  &
!$omp shared (nat,at,r2ab,dispm,itbl,c6abns,zvec,dzvec) &
!$omp shared (c6ab,dc6ab)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ij,ia,ja,k,ka,ik,jk,oorjk,oorik,atm,fdmp,  &
!$omp         r2ij,cij,oorij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp         dijatm,djkatm,dikatm,dtmp,dijfdmp,djkfdmp,dikfdmp,  &
!$omp         c9ijk,oor9ijk,dic9ijk,djc9ijk,dkc9ijk) &
!$omp shared(nat,at,r2ab,par,c6ab,dc6ab) &
!$omp reduction(+:eabc,dc6dr) reduction(-:dc6dcn)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2
            cjk   = par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)
            dtmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

            dc6dr(ij) = dc6dr(ij) + (atm*dijfdmp - dijatm*fdmp)*c9ijk
            dc6dr(ik) = dc6dr(ik) + (atm*dikfdmp - dikatm*fdmp)*c9ijk
            dc6dr(jk) = dc6dr(jk) + (atm*djkfdmp - djkatm*fdmp)*c9ijk
            dic9ijk = dc6ab(i,j)/c6ab(ij) + dc6ab(i,k)/c6ab(ik)
            djc9ijk = dc6ab(j,i)/c6ab(ij) + dc6ab(j,k)/c6ab(jk)
            dkc9ijk = dc6ab(k,j)/c6ab(jk) + dc6ab(k,i)/c6ab(ik)
            dc6dcn(i) = dc6dcn(i) - 0.5_wp*c9ijk*oor9ijk*dic9ijk
            dc6dcn(j) = dc6dcn(j) - 0.5_wp*c9ijk*oor9ijk*djc9ijk
            dc6dcn(k) = dc6dcn(k) - 0.5_wp*c9ijk*oor9ijk*dkc9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   if (present(eout)) eout=eabc

end subroutine dabcappr


!* here is the theory for the ATM-gradient (SAW, 180224)
! EABC = WA·WB·WC·DABC
! ∂EABC/∂X = ∂/∂X(WA·WB·WC·DABC)
!          = ∂WA/∂X·WB·WC·DABC + WA·∂WB/∂X·WC·DABC + WA·WB·∂WC/∂X·WC·DABC
!            + WA·WB·WC·∂DABC/∂X
! ∂/∂X =  ∂rAB/∂X·∂/∂rAB +  ∂rBC/∂X·∂/∂rBC +  ∂rCA/∂X·∂/∂rCA
!      = (δAX-δBX)∂/∂rAB + (δBX-δCX)∂/∂rBC + (δCX-δAX)∂/∂rCA
! ∂EABC/∂A = ∑A,ref ∑B,ref ∑C,ref
!            + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!            + WA·∂WB/∂rAB·WC·DABC
!            - WA·WB·∂WC/∂rCA·DABC
!            + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
! ∂EABC/∂B = ∑A,ref ∑B,ref ∑C,ref
!            - ∂WA/∂rAB·WB·WC·DABC
!            + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!            + WA·WB·∂WC/∂rBC·DABC
!            + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
! ∂EABC/∂C = ∑A,ref ∑B,ref ∑C,ref
!            + ∂WA/∂rCA·WB·WC·DABC
!            - WA·∂WB/∂rBC·WC·DABC
!            + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!            + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
! ∂WA/∂rAB = ∂CNA/∂rAB·∂WA/∂CNA w/ ζ=1 and WA=wA
! ATM = 3·cos(α)cos(β)cos(γ)+1
!     = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)/(r²BC·r²CA·r²AB)+1
! ∂ATM/∂rAB = 3/4(2r⁶AB-r⁶BC-r⁶CA-r⁴AB·r²BC-r⁴AB·r²CA+r⁴BC·r²CA+r²BC·r⁴CA)
!             /(r³AB·r²BC·r²CA)
! DABC = C9ABCns·f·ATM/(rAB·rBC·rCA)³
! f = 1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
! ∂(f/(r³AB·r³BC·r³CA)/∂rAB = 
!   ⅓·((6·(16-9)·(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶-9)·f²/(r⁴AB·r³BC·r³CA)
subroutine dabcgrad(nat,ndim,at,xyz,par,dcn,zvec,dzvec,itbl,g,eout)
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: dcn(nat,nat)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt,sum

   eabc = 0._wp

!$omp parallel private(ia,ja,ka,l,m,n, &
!$omp&         rij,rjk,rik,r2ij,r2jk,r2ik, &
!$omp&         cij,cjk,cik,cijk, &
!$omp&         fdmp,dtmp,oor9tmp,c9tmp, &
!$omp&         atm,r2ijk,c9ijk,oor9ijk,rijk, &
!$omp&         drij,drjk,drik,oorij,oorjk,oorik, &
!$omp&         dijfdmp,dikfdmp,djkfdmp, &
!$omp&         dijatm,dikatm,djkatm, &
!$omp&         dijoor9ijk,djkoor9ijk,dikoor9ijk, &
!$omp&         x1,x2,x3,x4,x5,x6,x7,x8,x9) &
!$omp&         reduction(+:g,eabc)
!$omp do schedule(runtime)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
!    --- all distances, cutoff radii ---
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         if(r2ij.gt.r_thr) cycle
         cij  = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2)
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
!           if(k.eq.j) cycle
!           if(i.eq.k) cycle
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            rjk   = xyz(:,k) - xyz(:,j)
            r2ik  = sum(rik**2)
            r2jk  = sum(rjk**2)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2)
            cjk   = (par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            x2 = 0._wp
            x4 = 0._wp
            x6 = 0._wp
            c9ijk = 0._wp

!       --- sum up all references ---
            do ii = 1, dispm%nref(ia) ! refs of A
               l = itbl(ii,i)
               do jj = 1, dispm%nref(ja) ! refs of B
                  m = itbl(jj,j)
                  do kk = 1, dispm%nref(ka) ! refs of C
                     n = itbl(kk,k)
                     if ((zvec(l)*zvec(m)*zvec(n)).lt.gw_thr) cycle
                     c9tmp = par%s9*thopi*trapzd(dispm%alpha(:,ii,ia)*dispm%alpha(:,jj,ja) &
                     &                      *dispm%alpha(:,kk,ka))

                     c9ijk = c9ijk + c9tmp*zvec(n)*zvec(m)*zvec(l)
!                --- intermediates ---
!                    ∂WA/∂CNA·WB·WC
                     x2 = x2 - dzvec(l)*zvec(m)*zvec(n)*c9tmp
!                    WA·∂WB/∂CNB·WC
                     x4 = x4 - dzvec(m)*zvec(l)*zvec(n)*c9tmp
!                    WA·WB·∂WC/∂CNC
                     x6 = x6 - dzvec(n)*zvec(m)*zvec(l)*c9tmp

                  enddo ! refs of k/C
               enddo ! refs of j/B
            enddo ! refs of i/A

!       --- geometrical term and r⁻³AB·r⁻³BC·r⁻³CA ---
!           ATM = 3·cos(α)cos(β)cos(γ)+1
!               = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)
!                 /(r²BC·r²CA·r²AB)+1
            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

!       --- damping function ---
!           1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**par%alp)
            dtmp = -(oth*six*par%alp*(thf*(cijk/rijk)**oth)**par%alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

!       --- intermediates ---
!           ∂WA/∂rAB·WB·WC·DABC = ∂CNA/∂rAB·(∂WA/∂CNA·WB·WC)·DABC
            x1 = x2*dcn(i,j)*( atm*fdmp )
!           ∂WA/∂rCA·WB·WC·DABC = ∂CNA/∂rCA·(∂WA/∂CNA·WB·WC)·DABC
            x2 = x2*dcn(i,k)*( atm*fdmp )
!           WA·∂WB/∂rBC·WC·DABC = ∂CNB/∂rBC·(WA·∂WB/∂rBC·WC)·DABC
            x3 = x4*dcn(j,k)*( atm*fdmp )
!           WA·∂WB/∂rAB·WC·DABC = ∂CNB/∂rAB·(WA·∂WB/∂rAB·WC)·DABC
            x4 = x4*dcn(i,j)*( atm*fdmp )
!           WA·WB·∂WC/∂rCA·DABC = ∂CNC/∂rCA·(WA·WB·∂WC/∂rCA)·DABC
            x5 = x6*dcn(i,k)*( atm*fdmp )
!           WA·WB·∂WC/∂rBC·DABC = ∂CNC/∂rBC·(WA·WB·∂WC/∂rBC)·DABC
            x6 = x6*dcn(j,k)*( atm*fdmp )
!           WA·WB·WC·∂DABC/∂rAB
            x7 = c9ijk*( atm*dijfdmp-dijatm*fdmp )
!           WA·WB·WC·∂DABC/∂rBC
            x8 = c9ijk*( atm*djkfdmp-djkatm*fdmp )
!           WA·WB·WC·∂DABC/∂rCA
            x9 = c9ijk*( atm*dikfdmp-dikatm*fdmp )

!       --- build everything together ---
            eabc = eabc + c9ijk*atm*fdmp

!           ∂rAB/∂A = -∂rAB/∂B
            drij = rij*oorij
!           ∂rBC/∂B = -∂rBC/∂C
            drjk = rjk*oorjk
!           ∂rCA/∂C = -∂rCA/∂A
            drik = rik*oorik

!           ∂EABC/∂A =
!           + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!           + WA·∂WB/∂rAB·WC·DABC
!           - WA·WB·∂WC/∂rCA·DABC
!           + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
            g(:,i) = g(:,i) + ( &
            &        + (x1+x4+x7)*drij &
            &        - (x2+x5+x9)*drik )
!           ∂EABC/∂B =
!           - ∂WA/∂rAB·WB·WC·DABC
!           + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!           + WA·WB·∂WC/∂rBC·DABC
!           + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
            g(:,j) = g(:,j) + ( &
            &        - (x1+x4+x7)*drij &
            &        + (x3+x6+x8)*drjk )
!           ∂EABC/∂C =
!           + ∂WA/∂rCA·WB·WC·DABC
!           - WA·∂WB/∂rBC·WC·DABC
!           + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!           + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
            g(:,k) = g(:,k) + ( &
            &        + (x2+x5+x9)*drik &
            &        - (x3+x6+x8)*drjk )

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp endparallel

   if (present(eout)) eout=eabc

end subroutine dabcgrad

subroutine dispmb(E,aw,xyz,oor6ab,nat)
   use mctc_la, only : syev,gemm
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(out) :: E

   integer  :: i,j,ii,jj,k
   integer  :: info
   real(wp) :: tau(3,3),spur(23),d_,d2,r(3),r2,alpha
   real(wp) :: two(23),atm(23),d3
   real(wp),allocatable :: T (:,:)
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: F (:,:)
   real(wp),allocatable :: F_(:,:)
   real(wp),allocatable :: d (:)
   real(wp),allocatable :: w (:)

   intrinsic :: sum,sqrt,minval,log

   allocate( T(3*nat,3*nat),  A(3*nat,3*nat), AT(3*nat,3*nat), &
   &         F(3*nat,3*nat), F_(3*nat,3*nat),  d(3*nat), &
   &         w(12*nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)*r(ii)-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo
         tau = tau*sqrt(oor6ab(i,j))
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
      enddo
   enddo

   !call prmat(6,T,3*nat,3*nat,'T')

   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0d0 
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
  &             3*nat,0.0_wp,F_,3*nat)
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,A, &
  &             3*nat,0.0_wp,AT,3*nat)

      F_ = F - AT

      d = 0.0d0
      call syev('N','U',3*nat,F_,3*nat,d,w,12*nat,info)
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(1x,''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(1x,''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
  &             3*nat,0.0_wp,F_,3*nat)
!     call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,AT, &
! &             3*nat,0.0_wp,A,3*nat)
       
      d_ = 1.0_wp; d2 = 0.0_wp!; d3 = 0.0_wp
      do i = 1, 3*nat
         d_ = d_ * d(i)
         d2 = d2 - F_(i,i)
!        d3 = d3 - A(i,i)
      enddo
      spur(k) = log(d_) - d2*0.5
!     two(k) = d2/2.0_wp
!     atm(k) = d3/3.0_wp
   enddo

   E = trapzd(spur)*ooTPI
   !print*,'     full contribution', trapzd(spur)*ooTPI
   !print*,' manybody contribution', trapzd(spur-two)*ooTPI
   !print*,'  twobody contribution', trapzd(two)*ootpi
   !print*,'threebody contribution', trapzd(atm)*ootpi

   deallocate(T,A,AT,F,F_,d)
end subroutine dispmb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!! This implementation of the MBD gradient is incorrect, DO NOT USE! !!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mbdgrad(nat,xyz,aw,daw,oor6ab,g,E)
   use mctc_la, only : syev,gemm
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: daw(23,3,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: E

   integer  :: i,ii,j,jj,k,kk,l,ll,m,n
   integer  :: info
   real(wp) :: spur(23)
   real(wp) :: fdmp,dfdmp
   real(wp) :: dum,dum2,ddum(3)
   real(wp) :: r(3),r2,drij(3)
   real(wp) :: alpha
   real(wp) :: tau(3,3),dtau(3,3,3)
   real(wp),allocatable :: dspur(:,:,:)
   real(wp),allocatable :: dA(:,:)
   real(wp),allocatable :: dT(:,:)
   real(wp),allocatable :: dAT(:,:)
   real(wp),allocatable :: nT(:,:,:) ! nabla T
   real(wp),allocatable :: T(:,:)
   real(wp),allocatable :: A(:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: dF(:,:)
   real(wp),allocatable :: F(:,:)
   real(wp),allocatable :: tmp1(:,:)
   real(wp),allocatable :: tmp2(:,:)
   real(wp),allocatable :: d(:)
   real(wp),allocatable :: w(:)

   intrinsic :: sum,sqrt,minval

   ! the MBD calculation in this RPA-like scheme needs quite a lot
   ! of memory, which might be problematic for big systems.
   ! We could use one or maybe to matrices less which would
   ! obfuscate the code
   allocate( T(3*nat,3*nat),A(3*nat,3*nat),AT(3*nat,3*nat), &
   &         F(3*nat,3*nat),dF(3*nat,3*nat), &
   &         tmp1(3*nat,3*nat),tmp2(3*nat,3*nat), &
   &         dT(3*nat,3*nat),dA(3*nat,3*nat), &
   &         dAT(3*nat,3*nat),nT(3*nat,3*nat,3), &
   &         d(3*nat),w(12*nat),dspur(23,3,nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

!-----------------------------------------------------------------------
!  interaction tensor setup                                    SAW 1708
!-----------------------------------------------------------------------
   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)**2-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo

         fdmp = sqrt(oor6ab(i,j))
         !print*, fdmp
         tau = tau*fdmp
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
         
!-----------------------------------------------------------------------
!        derivative of interaction tensor                      SAW 1802
!-----------------------------------------------------------------------
         dfdmp = -3*oor6ab(i,j)**2/fdmp*r2**2 ! *sqrt(r2)

         do ii = 1, 3
            dtau(ii,ii,ii) = (15*r(ii)*(0.6_wp*r2-r(ii)**2))/r2**2
            do jj = ii+1, 3
               dtau(jj,ii,ii) = (15*r(jj)*(0.2_wp*r2-r(ii)**2))/r2**2
               dtau(ii,jj,ii) = dtau(jj,ii,ii)
               dtau(ii,ii,jj) = dtau(jj,ii,ii)
               dtau(jj,jj,ii) = (15*r(ii)*(0.2_wp*r2-r(jj)**2))/r2**2
               dtau(jj,ii,jj) = dtau(jj,jj,ii)
               dtau(ii,jj,jj) = dtau(jj,jj,ii)
               do kk = jj+1, 3
                  dtau(ii,jj,kk) = -(15*r(ii)*r(jj)*r(kk))/r2**2
                  dtau(jj,kk,ii) = dtau(ii,jj,kk)
                  dtau(kk,ii,jj) = dtau(ii,jj,kk)
                  dtau(kk,jj,ii) = dtau(ii,jj,kk)
                  dtau(ii,kk,jj) = dtau(ii,jj,kk)
                  dtau(jj,ii,kk) = dtau(ii,jj,kk)
               enddo
            enddo
         enddo

         !drij = r/sqrt(r2)

         dtau(:,:,1) = ( dtau(:,:,1)*fdmp + tau*dfdmp*r(1) )
         dtau(:,:,2) = ( dtau(:,:,2)*fdmp + tau*dfdmp*r(2) )
         dtau(:,:,3) = ( dtau(:,:,3)*fdmp + tau*dfdmp*r(3) )

         !print*,dtau
         nT(3*i-2:3*i,3*j-2:3*j,1:3) = dtau
         nT(3*j-2:3*j,3*i-2:3*i,1:3) = dtau

      enddo
   enddo
        !call prmat(6,T,3*nat,3*nat,'T')

!-----------------------------------------------------------------------
! RPA-like calculation of MBD energy                           SAW 1708
!-----------------------------------------------------------------------
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²} ! this two-body energy has to be removed
!-----------------------------------------------------------------------
   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0_wp
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
           &     3*nat,0.0_wp,tmp1,3*nat)
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,A, &
           &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,AT,3*nat,3*nat,'AT')

      tmp1 = F - AT
        !call prmat(6,dF,3*nat,3*nat,'1-AT')

      d = 0.0d0
      call syev('V','U',3*nat,tmp1,3*nat,d,w,12*nat,info)
        !call prmat(6,tmp1,3*nat,3*nat,'F eigv.')
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(1x,''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(1x,''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

!     two-body contribution to energy
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
           &     3*nat,0.0_wp,tmp2,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')

       
!     get the trace
      dum = 1.0_wp; dum2 = 0.0_wp
      do i = 1, 3*nat
         dum = dum * d(i)
         dum2 = dum2 + tmp2(i,i)
      enddo
      spur(k) = log(dum) + 0.5_wp * dum2
      !print*, log(dum), 0.5_wp*dum2,spur(k)

!-----------------------------------------------------------------------
! MBD gradient calculation                                     SAW 1803
!-----------------------------------------------------------------------
! Some theory:
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
! ∇EMBD = 1/(2π)∫dω Tr{∇(log[1-AT])} = 1/(2π)∫dω Tr{(1-AT)⁻¹·(∇AT+A∇T)}
!       = 1/(2π)∫dω Tr{(1-A^½TA^½)⁻¹·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))}
! (1-AT)⁻¹ = U·Λ⁻¹·U†  w/ U⁻¹ = U†  if (1-AT) is symmetrized
! Aijkl = δij·δkl·∑(ref)Wi·α
! ∂/∂X Aijkl = δij·δkl·∑(ref)∂Wi/∂CNi·∂CNi/∂rij·∂rij/∂X·α
!-----------------------------------------------------------------------
! eigenvectors are still saved on tmp1, eigenvalues are still on d
! we (over)use tmp2 to hold intermediate results
!-----------------------------------------------------------------------

!     get the inverse of 1-AT by using its eigenvalues
      dF = 0.0_wp
      do i = 1, 3*nat
         dF(i,i) = 1._wp/d(i)
      enddo
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     (1-AT)⁻¹ = U·Λ⁻¹·U†
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,dF, &
           &     3*nat,0.0_wp,tmp2,3*nat)
!     call gemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     unfortunately this might not be enough, we have to substract the
!     twobody gradient for the dipole-dipole interaction, which is in
!     fact not easily accessable from here.
!     If this is correct,
!        E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
!       ∇E(2) = 1/(2π)∫dω Tr{(AT)·(∇AT+A∇T)}
!     then the MBD gradient w/o two-body contrib. could be represented by
!     dF = dF - AT
!     or more easier by the use of dgemm's beta by replacing the last
!     dgemm by:
!     call gemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,1.0_wp,AT,3*nat)
!     But we would have to consider AT instead of dF in the following code
!     which would be at least confusing

!-----------------------------------------------------------------------
!     Efficient MBD-Gradient calculation, O(N³)                SAW 1803
!-----------------------------------------------------------------------
      do l = 1, 3
         dA = 0.0_wp
         do i = 1, nat
            dA(3*i-2,3*i-2) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i-1,3*i-1) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i  ,3*i  ) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
         enddo

!        (∇A^½)TA^½
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,dA,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,0.0_wp,dAT,3*nat)

!        (∇A^½)TA^½+A^½(∇T)A^½ (please note the use of beta=1.0!)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,nT(:,:,l), &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        last term and we have: (∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,dA, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        by multiplying the prefactor we get to:
!        ((1-A^½TA^½)⁻¹+A^½TA^½)·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,dAT, &
              &     3*nat,0.0_wp,tmp2,3*nat)
!        now the other term of
!        [(1-A^½TA^½)⁻¹+A^½TA^½),((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)]
         call gemm('N','N',3*nat,3*nat,3*nat,-1.0_wp,dAT,3*nat,AT, &
              &     3*nat,1.0_wp,tmp2,3*nat)

         do i = 1, nat
            dspur(k,l,i) = dspur(k,l,i) &
            &              + tmp2(3*i-2,3*i-2) &
            &              + tmp2(3*i-1,3*i-1) &
            &              + tmp2(3*i  ,3*i  )
         enddo ! all atoms
         
      enddo ! cartesian parts

   enddo ! k, imaginary frequencies

   E = trapzd(spur)*ootpi
   !print*, E
   do i = 1, nat
      do j = 1, 3
         g(j,i) = g(j,i) + ootpi*trapzd(dspur(:,j,i))
         !print*,ootpi*trapzd(dspur(:,j,i))
      enddo
   enddo

   deallocate( T,A,AT,F,tmp1,tmp2,dT,dA,dAT,dF,nT,d,w,dspur )

end subroutine mbdgrad

! --- PBC
subroutine pbc_d4(nat,ndim,at,wf,g_a,g_c,covcn,gw,refc6)
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
      norm = 1._wp / norm
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),dispm%cn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical precision, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(dispm%cn(:dispm%nref(ia),ia)).eq.dispm%cn(ii,ia)) then
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

! compute D4 gradient under pbc
subroutine edisp_3d(mol,ndim,q,rep,atm_rep,r_thr,atm_thr,par,g_a,g_c,gw,refc6,mbd, &
      &             ed,etwo,embd)
   use mctc_constants
   use tbdef_molecule
   type(tb_molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   integer, intent(in)  :: rep(3)
   integer, intent(in)  :: atm_rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(in)  :: atm_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: ed
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: embd

   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den,dcndr

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   !  print'(" * Allocating local memory")'
   allocate( zetavec(ndim),zerovec(ndim), source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k)  &
!$omp shared (mol,dispm,itbl,g_a,g_c,q) &
!$omp shared (gw,zetavec,zerovec,r_thr)
!$omp do
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         zetavec(k) = zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gw(k)
         ! NEW: q=0 for ATM
         zerovec(k) =  zeta(g_a,gam(ia)*g_c,dispm%q(ii,ia)+iz,iz) * gw(k)

      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,disp,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10, &
!$omp         t,tx,ty,tz)  &
!$omp shared(mol,dispm,itbl,zetavec,refc6,par,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:ed)
!$omp do schedule(runtime)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do concurrent(tx = -rep(1):rep(1), &
            &       ty = -rep(2):rep(2), &
            &       tz = -rep(3):rep(3))
         ! cycle i with i interaction in same cell
         if (tx.eq.0.and.ty.eq.0.and.tz.eq.0) cycle
         rij = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ed = ed - c6ij*disp/2
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do concurrent(tx = -rep(1):rep(1), &
               &       ty = -rep(2):rep(2), &
               &       tz = -rep(3):rep(3))
            t = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
            rij = mol%xyz(:,i) - mol%xyz(:,j) + t
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ed = ed - c6ij*disp
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call abcappr_3d(mol,ndim,par,zerovec,refc6,itbl,atm_rep,atm_thr,eabc)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(etwo)) etwo = ed
   if (present(embd)) embd = eabc
   ed = ed + eabc

end subroutine edisp_3d

!> @brief calculates threebody dispersion gradient from C6 coefficients
subroutine abcappr_3d(mol,ndim,par,zvec,refc6,itbl,rep,r_thr,eabc)
   use tbdef_molecule
   use pbc_tools
   type(tb_molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer  :: ij
   real(wp),allocatable :: c6ab(:)
   real(wp) :: r,r2ij
   real(wp) :: c6ij
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n*(mol%n+1)/2), source = 0.0_wp )

   eabc = 0.0_wp

   !$OMP parallel default(none) &
   !$omp private(i,ia,j,ja,ij,r2ij,c6ij,k,l,r)  &
   !$omp shared (mol,dispm,itbl,refc6,zvec) &
   !$omp shared (c6ab)
   !$omp do schedule(runtime)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(ij) = c6ij

      do j = 1, i-1
         !        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

         !        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
   !$omp enddo
   !$omp end parallel

   call abcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep, &
      &                             mol%lattice,c6ab,eabc)
   !call abcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,eabc)
   !call abcappr_3d_bvk(mol,par,rep,r_thr,c6ab,eabc)

end subroutine abcappr_3d

subroutine abcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,eabc)
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: eabc

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1

   real(wp) rij(3),rik(3),rjk(3)
   real(wp), allocatable,dimension(:,:,:,:) ::  dc6dr  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax

   allocate(dc6dr(-rep(3):rep(3),-rep(2):rep(2), &
      &           -rep(1):rep(1),nat*(nat+1)/2))
   dc6dr = 0.0_wp


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)

   do iat=3,nat
      do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(ij)
         do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(ik)
            c6jk=c6ab(jk)
            c9=-sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            do concurrent (tx=-rep(1):rep(1), &
                  &        ty=-rep(2):rep(2), &
                  &        tz=-rep(3):rep(3))
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij2=SUM((ijvec+t)*(ijvec+t))
               if(rij2.gt.thr)cycle

               !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))


               do concurrent (sx=repmin(1):repmax(1), &
                     &        sy=repmin(2):repmax(2), &
                     &        sz=repmin(3):repmax(3))
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik2=SUM((ikvec+s)*(ikvec+s))
                  if(rik2.gt.thr)cycle

                  dumvec=jkvec+s-t
                  rjk2=SUM(dumvec*dumvec)
                  if(rjk2.gt.thr)cycle
                  !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
                  !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang
                  !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
                  !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  eabc=eabc-r9ijk*c9

               enddo !sz
            enddo !tx
         enddo !kat
      enddo !jat
   enddo !iat

   ! Now the interaction with jat=iat of the triples iat,iat,kat
   DO iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(ij)
      DO kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(ik)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik
         do concurrent (tx=-rep(1):rep(1), &
               &        ty=-rep(2):rep(2), &
               &        tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=t
            rij2=SUM(dumvec*dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
               !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))


               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sx

         enddo !tx
      enddo !kat
   enddo !iat
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

   do iat=2,nat
      do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(ij)
         c6ik=c6ij

         c6jk=c6ab(jk)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-sqrt(c6ij*c6ik*c6jk)
         do concurrent(tx=-rep(1):rep(1), &
               &       ty=-rep(2):rep(2), &
               &       tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=ijvec+t
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5
               IF (tx.eq.sx .and. ty.eq.sy  &
                  .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle
               !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

               !              if (rij*rjk*rik.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6._wp*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sx

         enddo !tx
      enddo !kat
   enddo !iat


   ! and finally the self interaction iat=jat=kat all

   do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(ij)
      c6ik=c6ij
      c6jk=c6ij
      c9=-sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      do concurrent ( tx=-rep(1):rep(1), &
            &         ty=-rep(2):rep(2), &
            &         tz=-rep(3):rep(3))
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle !IF iat and jat are the same then cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         dumvec=t
         dumvec=dumvec*dumvec
         rij2=SUM(dumvec)
         if(rij2.gt.thr)cycle
         !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

         do concurrent (sx=repmin(1):repmax(1), &
               &        sy=repmin(2):repmax(2), &
               &        sz=repmin(3):repmax(3))
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle !IF iat and kat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty)  &
               .and. (sz.eq.tz)) cycle      !If kat and jat are the same then cycle

            ! every result * 1/6 becaues every triple is counted twice due to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            dumvec=s
            dumvec=dumvec*dumvec
            rik2=SUM(dumvec)
            if(rik2.gt.thr)cycle
            !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

            dumvec=jkvec+s-t
            dumvec=dumvec*dumvec
            rjk2=SUM(dumvec)
            if(rjk2.gt.thr)cycle
            !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

            rijk2=(rij2*rjk2*rik2)
            !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
            !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            eabc=eabc-c9*r9ijk

         enddo !sx
      enddo !tx


   enddo !iat


end subroutine abcappr_3d_dftd3_like_style

! compute D4 gradient under pbc
subroutine dispgrad_3d(mol,ndim,q,cn,dcndr,dcndL,rep,atm_rep,r_thr,atm_thr,par, &
      &                wf,g_a,g_c,refc6,mbd,g,sigma,eout,dqdr,dqdL,aout)
   use mctc_constants
   use tbdef_molecule
   use pbc_tools
   type(tb_molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   real(wp),intent(in)  :: cn(mol%n)
   real(wp),intent(in)  :: dcndr(3,mol%n,mol%n)
   real(wp),intent(in)  :: dcndL(3,3,mol%n)
   integer, intent(in)  :: rep(3),atm_rep(3)
   real(wp),intent(in)  :: r_thr,atm_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,mol%n)
   real(wp),intent(inout)        :: sigma(3,3)
   real(wp),intent(out),optional :: eout
   real(wp),intent(in), optional :: dqdr(3,mol%n,mol%n+1)
   real(wp),intent(in), optional :: dqdL(3,3,mol%n+1)
   real(wp),intent(out),optional :: aout(23,mol%n)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:,:)
   real(wp),allocatable :: dzvec(:,:)
   real(wp),allocatable :: gw(:,:)
   real(wp),allocatable :: dgw(:,:)
   real(wp),allocatable :: dc6dq(:)
   real(wp),allocatable :: dzdq(:,:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   !  print'(" * Allocating local memory")'
   allocate( dc6dcn(mol%n),r2ab(mol%n*(mol%n+1)/2),dc6dq(mol%n),dzdq(7,mol%n), &
      &      zvec(7,mol%n),dzvec(7,mol%n),gw(7,mol%n),dgw(7,mol%n), &
      &      source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   call weight_references(mol%n, mol%at, g_a, g_c, wf, q, cn, &
      &                   zvec, gw, dzvec, dgw, dzdq)

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,dic6ii,dic6ij,djc6ij,disp,ddisp,dizii,dizij,djzij,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10,door6,door8,door10, &
!$omp         t,tx,ty,tz,dtmp,drdx)  &
!$omp shared(mol,dispm,itbl,zvec,dzvec,refc6,par,dzdq,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:dc6dq,dc6dcn,ed,g,sigma)
!$omp do schedule(runtime)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      dizij  = 0.0_wp
      djzij  = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zvec(ii,i) *  zvec(jj,i) * refc6(k,l)
            dic6ij = dic6ij + dzvec(ii,i) *  zvec(jj,i) * refc6(k,l)
            djc6ij = djc6ij +  zvec(ii,i) * dzvec(jj,i) * refc6(k,l)
            dizij  = dizij  +  dzdq(ii,i) *  zvec(jj,i) * refc6(k,l)
            djzij  = djzij  +  zvec(ii,i) *  dzdq(jj,i) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do concurrent(tx = -rep(1):rep(1), &
            &       ty = -rep(2):rep(2), &
            &       tz = -rep(3):rep(3), &
            &       tx.ne.0.or.ty.ne.0.or.tz.ne.0)
         ! cycle i with i interaction in same cell
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
         ed = ed - c6ij*disp/2
         ! save this
         dtmp = c6ij*ddisp
         dc6dq(i)  = dc6dq(i)  + (dizij  + djzij )*disp/2
         dc6dcn(i) = dc6dcn(i) + (dic6ij + djc6ij)*disp/2
         drdx = rij/r
         sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)/2
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         dizij  = 0.0_wp
         djzij  = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zvec(ii,i) *  zvec(jj,j) * refc6(k,l)
               dic6ij = dic6ij + dzvec(ii,i) *  zvec(jj,j) * refc6(k,l)
               djc6ij = djc6ij +  zvec(ii,i) * dzvec(jj,j) * refc6(k,l)
               dizij  = dizij  +  dzdq(ii,i) *  zvec(jj,j) * refc6(k,l)
               djzij  = djzij  +  zvec(ii,i) *  dzdq(jj,j) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do concurrent(tx = -rep(1):rep(1), &
               &       ty = -rep(2):rep(2), &
               &       tz = -rep(3):rep(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,i) - mol%xyz(:,j) + matmul(mol%lattice,t)
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            door6 = -6*r2**2*r*oor6**2
            door8 = -8*r2**3*r*oor8**2
            door10 = -10*r2**4*r*oor10**2
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
            ed = ed - c6ij*disp
            ! save this
            dc6dq(i)  = dc6dq(i)  + dizij  *disp
            dc6dq(j)  = dc6dq(j)  + djzij  *disp
            dc6dcn(i) = dc6dcn(i) + dic6ij *disp
            dc6dcn(j) = dc6dcn(j) + djc6ij *disp
            dtmp = c6ij*ddisp
            drdx = rij/r
            g(:,i) = g(:,i) - dtmp * drdx
            g(:,j) = g(:,j) + dtmp * drdx

            sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call dabcappr_3d(mol,ndim,par,gw,dgw,refc6,itbl,atm_rep,atm_thr, &
         &             g,sigma,dc6dcn,eabc)
      !dc6dcn = 0.0_wp
   endif

   if(present(dqdr)) then
      ! handle dqdr  :: gradient is exact e-11
      call dgemv('n',3*mol%n,mol%n,-1.0_wp,dqdr,3*mol%n,dc6dq,1,1.0_wp,g,1)
   endif
   if (mol%npbc > 0) then
      if(present(dqdL)) then
         ! handle dqdL
         call dgemv('n',3*3,mol%n,-1.0_wp,dqdL,3*3,dc6dq,1,1.0_wp,sigma,1)
      endif
   endif

   ! always handle dcndr :: gradient is exact e-11
   call dgemv('n',3*mol%n,mol%n,-1.0_wp,dcndr,3*mol%n,dc6dcn,1,1.0_wp,g,1)
   if (mol%npbc > 0) then
      call dgemv('n',3*3,mol%n,-1.0_wp,dcndL,3*3,dc6dcn,1,1.0_wp,sigma,1)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, mol%n
         ia = mol%at(i)
         do ii = 1, dispm%nref(ia)
            aout(:,i) = aout(:,i) + zvec(ii,i) * dispm%alpha(:,ii,ia)
         enddo
      enddo
   endif
end subroutine dispgrad_3d

!> @brief calculates threebody dispersion gradient from C6 coefficients
subroutine dabcappr_3d(mol,ndim,par,zvec,dzvec,refc6,itbl,rep,r_thr, &
      &                g,sigma,dc6dcn,eabc)
   use tbdef_molecule
   use pbc_tools
   type(tb_molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(:,:)
   real(wp),intent(in)  :: dzvec(:,:)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(inout) :: g(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: dc6dcn(mol%n)
   real(wp),intent(out)   :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:,:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik,r
   real(wp) :: cii,cij,cjk,cik,ciii,ciij,cijk
   real(wp) :: c9iii,c9iij,c9ijk,oor9ijk,rijk
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9iii,dic9iij,djc9iij,dic9ijk,djc9ijk,dkc9ijk
   real(wp),parameter :: zero(3) = [0.0_wp,0.0_wp,0.0_wp]
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n,mol%n),dc6ab(mol%n,mol%n), source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij,k,l,r)  &
!$omp shared (mol,dispm,itbl,refc6,zvec,dzvec) &
!$omp shared (c6ab,dc6ab)
!$omp do schedule(runtime)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(ii,i)*zvec(jj,i)*refc6(k,l)
            dic6ij = dic6ij + dzvec(ii,i)*zvec(jj,i)*refc6(k,l)
            djc6ij = djc6ij + zvec(ii,i)*dzvec(jj,i)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(i,i) = c6ij
      dc6ab(i,i) = dic6ij! + djc6ij

      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(ii,i)*zvec(jj,j)*refc6(k,l)
               dic6ij = dic6ij + dzvec(ii,i)*zvec(jj,j)*refc6(k,l)
               djc6ij = djc6ij + zvec(ii,i)*dzvec(jj,j)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(i,j) = c6ij
         c6ab(j,i) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

   call dabcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep,&
      &    mol%lattice,c6ab,dc6ab,eabc,dc6dcn,g,sigma)
   !call dabcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,dc6ab,g,dc6dcn,eabc)
   !call dabcappr_3d_bvk(mol,par,rep,r_thr,c6ab,dc6ab,g,dc6dcn,eabc)

end subroutine dabcappr_3d

subroutine dabcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,dc6ab, &
      &    eabc,dc6dcn,g,sigma)
   use mctc_constants
   use pbc_tools
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat,nat)
   real(wp),intent(in) :: dc6ab(nat,nat)    !dC6(iat,jat)/cCN(iat) in dc6ab(i,j) for ABC-grad
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: eabc
   real(wp),intent(inout) :: dc6dcn(nat)

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1
   real(wp) :: rco,den,tmp,dtmp

   real(wp) rij(3),rik(3),rjk(3)
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   real(wp) :: r3,g3(3,3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)

!$omp parallel default(none) &
!$omp shared(nat,xyz,c6ab,rep,par,at,dlat,dc6ab,thr) &
!$omp private(iat,jat,ij,ijvec,c6ij,kat,ik,jk,ikvec,jkvec,c6ik,c6jk,c9, &
!$omp&        cij,cik,cjk,cijk,tx,ty,tz,repmin,repmax,t,rij2,sx,sy,sz, &
!$omp&        s,rik2,vec,rjk2,rijk2,rijk,fdmp,rijk3,ang,r9ijk,dfdmp, &
!$omp&        r,dang,tmp1,dc9,rij,rik,rjk,r3,g3) &
!$omp reduction(+:eabc,dc6dcn,sigma,g)
!$omp do schedule(runtime)
   iAt_ijk: do iat=3,nat
      jAt_ijk: do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(iat,jat)
         kAt_ijk: do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(iat,kat)
            c6jk=c6ab(jat,kat)
            c9=-sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            r3 = 0.0_wp
            g3 = 0.0_wp
            do concurrent (tx=-rep(1):rep(1), &
                  &        ty=-rep(2):rep(2), &
                  &        tz=-rep(3):rep(3))
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij=ijvec+t
               rij2=SUM(rij*rij)
               if(rij2.gt.thr)cycle


               do concurrent (sx=repmin(1):repmax(1), &
                     &        sy=repmin(2):repmax(2), &
                     &        sz=repmin(3):repmax(3))
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik = ikvec+s
                  rik2=SUM(rik*rik)
                  if(rik2.gt.thr)cycle

                  rjk = jkvec+s-t
                  rjk2=SUM(rjk*rjk)
                  if(rjk2.gt.thr)cycle
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  r3 = r3 + r9ijk
                  !
                  !start calculating the gradient components dfdmp, dang and dc9

                  !dfdmp is the same for all three distances
                  dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

                  !start calculating the derivatives of each part w.r.t. r_ij
                  r=sqrt(rij2)


                  dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                     +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                     -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  g3(:,1) = g3(:,1) - tmp1*rij
                  g3(:,2) = g3(:,2) + tmp1*rij
                  sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

                  !start calculating the derivatives of each part w.r.t. r_ik

                  r=sqrt(rik2)


                  dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                     +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                     -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  !                 tmp1=-dc9
                  g3(:,1) = g3(:,1) - tmp1*rik
                  g3(:,3) = g3(:,3) + tmp1*rik
                  sigma = sigma + outer_prod_3x3(rik,rik)*tmp1

                  !
                  !start calculating the derivatives of each part w.r.t. r_jk

                  r=sqrt(rjk2)

                  dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                     +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                     -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  g3(:,2) = g3(:,2) - tmp1*rjk
                  g3(:,3) = g3(:,3) + tmp1*rjk
                  sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

               enddo !sz
            enddo !tx
            eabc = eabc - r3*c9
            g(:,iat) = g(:,iat) + g3(:,1)
            g(:,jat) = g(:,jat) + g3(:,2)
            g(:,kat) = g(:,kat) + g3(:,3)
            !calculating the CN derivative dE_disp(ijk)/dCN(i)
            dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
            dc6dcn(iat) = dc6dcn(iat) + r3*dc9
            dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
            dc6dcn(jat) = dc6dcn(jat) + r3*dc9
            dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
            dc6dcn(kat) = dc6dcn(kat) + r3*dc9
         enddo kAt_ijk
      enddo jAt_ijk
   enddo iAt_ijk
!$omp enddo

!$omp do schedule(runtime)
   ! Now the interaction with jat=iat of the triples iat,iat,kat
   iAt_iik: do iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(iat,jat)
      kAt_iik: do kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(iat,kat)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         g3 = 0.0_wp
         r3 = 0.0_wp
         do concurrent (tx=-rep(1):rep(1), &
               &        ty=-rep(2):rep(2), &
               &        tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            rij=t
            rij2=SUM(rij*rij)
            if(rij2.gt.thr)cycle

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               rik=ikvec+s
               rik2=SUM(rik*rik)
               if(rik2.gt.thr)cycle

               rjk=jkvec+s-t
               rjk2=SUM(rjk*rjk)
               if(rjk2.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               r3 = r3 + r9ijk

               !              iat=jat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                  -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rik
               g3(:,3) = g3(:,3) + tmp1*rik
               sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
               !
               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r

               g3(:,2) = g3(:,2) - tmp1*rjk
               g3(:,3) = g3(:,3) + tmp1*rjk
               sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

            enddo !sx

         enddo !tx
         eabc = eabc - r3*c9
         g(:,iat) = g(:,iat) + g3(:,1)
         g(:,jat) = g(:,jat) + g3(:,2)
         g(:,kat) = g(:,kat) + g3(:,3)
         !calculating the CN derivative dE_disp(ijk)/dCN(i)
         dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
         dc6dcn(iat) = dc6dcn(iat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
         dc6dcn(jat) = dc6dcn(jat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
         dc6dcn(kat) = dc6dcn(kat) + r3*dc9
      enddo kAt_iik
   enddo iAt_iik
!$omp enddo
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

!$omp do schedule(runtime)
   iAt_ijj: do iat=2,nat
      jAt_ijj: do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(iat,jat)
         c6ik=c6ij

         c6jk=c6ab(jat,kat)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-sqrt(c6ij*c6ik*c6jk)
         g3 = 0.0_wp
         r3 = 0.0_wp
         do concurrent(tx=-rep(1):rep(1), &
               &       ty=-rep(2):rep(2), &
               &       tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            rij=ijvec+t
            rij2=SUM(rij*rij)
            if(rij2.gt.thr)cycle

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5
               if (tx.eq.sx .and. ty.eq.sy .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               rik=ikvec+s
               rik2=SUM(rik*rik)
               if(rik2.gt.thr)cycle

               rjk=jkvec+s-t
               rjk2=SUM(rjk*rjk)
               if(rjk2.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               r3 = r3 + r9ijk

               !              jat=kat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0_wp*rjk2*rik2+3.0_wp*rik2**2) &
                  -5.0_wp*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rij
               g3(:,2) = g3(:,2) + tmp1*rij
               sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rik
               g3(:,3) = g3(:,3) + tmp1*rik
               sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
               !
               !start calculating the derivatives of each part w.r.t. r_jk
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0_wp*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

            enddo !sx
         enddo !tx
         eabc = eabc - r3*c9
         g(:,iat) = g(:,iat) + g3(:,1)
         g(:,jat) = g(:,jat) + g3(:,2)
         g(:,kat) = g(:,kat) + g3(:,3)
         !calculating the CN derivative dE_disp(ijk)/dCN(i)
         dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
         dc6dcn(iat) = dc6dcn(iat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
         dc6dcn(jat) = dc6dcn(jat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
         dc6dcn(kat) = dc6dcn(kat) + r3*dc9
      enddo jAt_ijj
   enddo iAt_ijj
!$omp enddo

   ! And finally the self interaction iat=jat=kat all

!$omp do schedule(runtime)
   iAt_iii: do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(iat,jat)
      c6ik=c6ij
      c6jk=c6ij
      c9=-sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      r3 = 0.0_wp
      do concurrent ( tx=-rep(1):rep(1), &
            &         ty=-rep(2):rep(2), &
            &         tz=-rep(3):rep(3))
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         ! IF iat and jat are the same then cycle
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         rij=t
         rij2=SUM(rij*rij)
         if(rij2.gt.thr)cycle

         do concurrent (sx=repmin(1):repmax(1), &
               &        sy=repmin(2):repmax(2), &
               &        sz=repmin(3):repmax(3))
            ! if iat and kat are the same then cycle
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle
            ! If kat and jat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty) .and. (sz.eq.tz)) cycle

            ! every result * 1/6 becaues every triple is counted twice due
            ! to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            rik=s
            rik2=SUM(rik*rik)
            if(rik2.gt.thr)cycle

            rjk=jkvec+s-t
            rjk2=SUM(rjk*rjk)
            if(rjk2.gt.thr)cycle

            rijk2=(rij2*rjk2*rik2)

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            r3 = r3 + r9ijk

            !                          iat=jat=kat
            dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
            !start calculating the derivatives of each part w.r.t. r_ij

            r=sqrt(rij2)
            dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
               +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
               -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
               /(r*rijk3*rijk2)


            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

            !start calculating the derivatives of each part w.r.t. r_ik

            r=sqrt(rik2)

            dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
               +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
               -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
            !
            !start calculating the derivatives of each part w.r.t. r_jk

            r=sqrt(rjk2)
            dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
               +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
               -5.0*(rik2-rij2)**2*(rik2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

         enddo !sx
      enddo !tx
      eabc = eabc - r3*c9
      !calculating the CN derivative dE_disp(ijk)/dCN(i)
      dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
      dc6dcn(iat) = dc6dcn(iat) + r3*dc9
      dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
      dc6dcn(jat) = dc6dcn(jat) + r3*dc9
      dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
      dc6dcn(kat) = dc6dcn(kat) + r3*dc9
   enddo iAt_iii
!$omp enddo
!$omp end parallel

end subroutine dabcappr_3d_dftd3_like_style

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(nat, atoms, g_a, g_c, wf, q, cn, &
      &                      zetavec, zerovec, zetadcn, zerodcn, zetadq)
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

   zetavec = 0.0_wp
   zerovec = 0.0_wp
   zetadcn = 0.0_wp
   zerodcn = 0.0_wp
   zetadq  = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)

      zi = zeff(ati)
      gi = g_c * gam(ati)

      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, dispm%nref(ati)
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            norm = norm + gw
            dnorm = dnorm + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo
      end do
      norm = 1.0_wp / norm
      do iref = 1, dispm%nref(ati)
         expw = 0.0_wp
         expd = 0.0_wp
         do icount = 1, dispm%ncount(iref, ati)
            twf = icount * wf
            gw = cngw(twf, cn(iat), dispm%cn(iref, ati))
            expw = expw + gw
            expd = expd + 2*twf*(dispm%cn(iref, ati) - cn(iat)) * gw
         enddo

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(dispm%cn(:dispm%nref(ati), ati)) &
               & == dispm%cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zetavec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerovec(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         zetadcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * dgwk
         zetadq(iref, iat) = dzeta(g_a,gi,dispm%q(iref,ati)+zi,q(iat)+zi) * gwk
         zerodcn(iref, iat) = zeta(g_a,gi,dispm%q(iref,ati)+zi,zi) * dgwk

      end do
   end do

end subroutine weight_references

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(nat, atoms, zetavec, zetadcn, zetadq, c6, dc6dcn, dc6dq)
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

   c6 = 0.0_wp
   dc6dcn = 0.0_wp
   dc6dq = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         dc6dqi = 0.0_wp
         dc6dqj = 0.0_wp
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
end subroutine get_atomic_c6

end module tbmod_dftd4
