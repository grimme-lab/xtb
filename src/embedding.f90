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

module xtb_embedding
   use, intrinsic :: iso_fortran_env, only : output_unit
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data, only : TCoulombData
   implicit none
!! ========================================================================
!  Embedding Stuff -- for QM/MM and subsystem-DFTB
!! ------------------------------------------------------------------------
!  in our tight binding model we can do much more than just point charge
!  embedding. Since we have a multipole electrostatic we can also use
!  these multipoles (dipole and quadrupoles). Also there is the dispersion
!  potential so we can embed the system in polarizibilities (coordination).
!  Given this fact the pcharges input will be more complicated
!  1> <npc>
!  2> <q> <x> <y> <z> [at [dx dy dz [qxx qyy qzz qxy qxz qxy [CN]]]]
!  ...
!  where q is the partial charge and x,y,z are it respective coordinates
!  at the moment the atom type is optional and by default nitrogen (7).
!! ------------------------------------------------------------------------
!  by this embedding scheme we can get in a subsystem-DFTB scheme
!  gradients for:
!  - monopole electrostatic
!  - anisotropic electrostatic
!  - dispersion
!  still missing is
!  - polynomial wavefunction terms
!  - repulsion terms
!  - shellwise electrostatic
!  - CN-shift of H0 (maybe via divide and conquere with P)
!  - nonadditive dispersion (ATM)
!  GBSA might need some special attention
!! ========================================================================

   integer,private,parameter :: p_str_length = 48
   integer,private,parameter :: p_arg_length = 24

contains

subroutine init_pcem
   use xtb_setparam
   implicit none
   select case(pcem_interface)
   case default
      if (.not.allocated(pcem_file)) pcem_file = "pcharge"
      if (.not.allocated(pcem_grad)) pcem_grad = "pcgrad"
   end select
end subroutine init_pcem

subroutine read_pcem(iunit,env,pcem,jData)
   use xtb_mctc_strings
   use xtb_mctc_convert
   use xtb_type_pcem
   use xtb_type_environment, only : TEnvironment
   use xtb_setparam
   use xtb_readin, only : getline => strip_line
   implicit none
   integer,      intent(in)    :: iunit
   type(TEnvironment), intent(inout) :: env
   type(tb_pcem),intent(inout) :: pcem
   type(TCoulombData), intent(in) :: jData

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   character(len=:),allocatable :: line
   integer  :: err
   integer  :: npc,n
   integer  :: i,iat
   real(wp) :: xyz(3),q,gami
   real(wp) :: conv

   conv = 1.0_wp ! input is usually in *Bohr*
   if (pcem_orca) then
      conv = aatoau ! except when we try to read ORCA's pc-files
   endif

   ! file is not open and error was not catched from caller
   if (iunit.eq.-1) return

   ! reset file position
   rewind(iunit)

   ! load first line to buffer and attempt to parse
   call getline(iunit,line,err)
   read(line,*,iostat=err) npc
   if (err.ne.0) then
      return
   endif

   call pcem%allocate(npc)

   ! we now the number of lines that should follow, but we don't care
   ! -> read until EOF and count lines
   n = 0
   do

      ! load a new line into buffer
      ! exit if EOF occurs
      call getline(iunit,line,err)
      if (err.ne.0) exit

      ! we skip empty lines (might been commented out)
      if (len(line).eq.0) cycle

      ! try to parse the lines into argv, this has a hard limit on length and items
      ! let's assume a well-behaving user for now
      call parse(line,' ',argv,narg)

      ! user provided crap input, bad user
      if (narg.lt.4) then
         call raise('E',"Not enough entries for PC, please check!")
      endif

      ! we count this line
      n = n+1
      ! before we start parsing we check the number of lines
      if (n > npc) then ! bad user
         call raise('E',"Wrong dimension input for PC, too many lines provided")
      endif

      ! first position should be the partial charge
      read(argv(1),*,iostat=err) q
      ! get the position from the next three entries
      do i = 1, 3
         read(argv(i+1),*,iostat=err) xyz(i)
      enddo

      ! now our well-behaving user has three options:
      ! -> provide no fifth argument -> fallback to default from $pcem group
      ! -> provide a atom number     -> use gfn_method gam value (aoparam)
      ! -> provide a gam value       -> double check this one!!!
      if (narg.eq.5) then
         read(argv(5),*,iostat=err) gami
         if (err.ne.0) then
            call elem(argv(5),iat)
            if (iat.eq.0) then ! so much for the well-behaved user
               call raise('E',"Invalid PC input: '"//trim(argv(5))//"'")
            endif
            gami = jData%chemicalHardness(iat)
         endif
         ! GFN0-xTB has negative gam-values since they are internally combined
         ! with the atom radius for the full chemical hardness (overall > 0).
         ! This is certainly a problem here, but since GFN0-xTB does
         ! not provide the possibility to use a PC potential, we don't care...
         ! But if YOU are reading this comment and are in charge of implementing
         ! a PC potential for GFN0-xTB, than this is YOUR problem now.
         if (gami < 0.0_wp) then
            call raise('S',"Found negative gam-value in user input: '"//&
               &            trim(argv(5))//"'")
         endif
      else
         ! we trust the dummy atom from xcontrol, since xcontrol could already
         ! have handled potential errors, which is usually not the case...
         gami = jData%chemicalHardness(pcem_dummyatom)
      endif
      pcem%xyz(:,n) = xyz*conv
      pcem%q(n) = q
      pcem%gam(n) = gami
   enddo

   ! user screwed up, we could rescue this by resetting npc
   ! or we rub it in his/her face :)
   if (n /= npc) then
      call raise('E',"Wrong dimension input for PC, too few lines provided")
   endif

end subroutine read_pcem

!! ========================================================================
!  J potentials for GFN1 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn1(jData,n,pcem,nshell,at,xyz,alphaj,Vpc)
   use xtb_type_pcem
   implicit none
   type(TCoulombData), intent(in) :: jData
   integer, intent(in)  :: n
   type(tb_pcem),intent(inout) :: pcem
   integer, intent(in)  :: nshell(:)
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(inout) :: Vpc(:)
   integer  :: ish,ii,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,rab

   ii = 0
   do iat = 1, n
      ati = at(iat)
      do ish = 1, nshell(ati)
         gi = jData%shellHardness(ish,ati)
         eh1 = 0.0_wp
         do kk = 1, pcem%n
            gj = pcem%gam(kk)
            rab = sqrt(sum((pcem%xyz(:,kk) - xyz(:,iat))**2))
            xj = 2.0_wp/(1./gi+1./gj)
            dum = 1.0_wp/(rab**alphaj+1._wp/xj**alphaj)**(1._wp/alphaj)
            eh1 = eh1+pcem%q(kk)*dum
         end do
         Vpc(ii+ish) = eh1
      end do
      ii = ii + nshell(ati)
   end do

end subroutine jpot_pcem_gfn1

!! ========================================================================
!  J potentials for GFN2 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn2(jData,n,pcem,nshell,at,xyz,Vpc)
   use xtb_type_pcem
   implicit none
   type(TCoulombData), intent(in) :: jData
   integer, intent(in)  :: n
   type(tb_pcem),intent(in) :: pcem
   integer, intent(in)  :: nshell(:)
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(inout) :: Vpc(:)
   integer  :: ish,ii,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,r2

   ii = 0
   do iat = 1, n
      ati = at(iat)
      do ish = 1, nshell(ati)
         gi = jData%shellHardness(ish,ati)
         eh1 = 0.0_wp
         do kk = 1, pcem%n
            gj = pcem%gam(kk)
            r2 = sum((pcem%xyz(:,kk) - xyz(:,iat))**2)
            xj = 0.5_wp*(gi+gj)
            dum = 1.0_wp/sqrt(r2+1._wp/xj**2)
            eh1 = eh1+pcem%q(kk)*dum
         enddo
         Vpc(ii+ish) = eh1
      enddo
      ii = ii + nshell(ati)
   enddo

end subroutine jpot_pcem_gfn2

!! ========================================================================
!  point charge embedding ES subroutine for SCC
!! ========================================================================
pure subroutine electro_pcem(nshell,dqsh,Vpc,es,scc)
   use xtb_mctc_convert, only : evtoau
   implicit none
   integer, intent(in)    :: nshell
   real(wp),intent(in)    :: Vpc(nshell),dqsh(nshell)
   real(wp),intent(out)   :: es
   real(wp),intent(inout) :: scc

   integer :: i,j

!  point charge contributions to ES
   es =0.0d0
   do i=1,nshell
      es =es + dqsh(i)*Vpc(i)
   enddo

!  Etotal
   scc = scc + es

end subroutine electro_pcem

!! ========================================================================
!  point charge embedding gradient for GFN1
!! ========================================================================
subroutine pcem_grad_gfn1(jData,g,gpc,n,pcem,at,nshell,xyz,alphaj,qsh)
   use xtb_type_pcem
   implicit none
   type(TCoulombData), intent(in) :: jData
   integer, intent(in) :: n
   type(tb_pcem),intent(in) :: pcem
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,pcem%n)
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell(:)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: alphaj
   real(wp),intent(in) :: qsh(:)

   integer  :: ish,ii,j,iat,ati
   real(wp) :: dr(3)
   real(wp) :: gi,gj,r2,rr,yy,ff

   ii = 0
   do iat = 1, n
      ati = at(iat)
      do ish = 1, nshell(ati)
         gi = jData%shellHardness(ish,ati)
         do j = 1, pcem%n
            dr = xyz(:,iat) - pcem%xyz(:,j)
            r2 = sum(dr**2)
            gj = pcem%gam(j)
            rr = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)
            rr = 1.0_wp/rr**alphaj
            ff = r2**(alphaj/2.0_wp-1.0_wp)*&
               &(r2**(alphaj*0.5_wp)+rr)**(-1.0_wp/alphaj-1.0_wp)
            yy = ff*qsh(ii+ish)*pcem%q(j)
            g(:,iat)=g(:,iat)-dr*yy
            gpc(:,j)=gpc(:,j)+dr*yy
         end do
      end do
      ii = ii + nshell(ati)
   end do

end subroutine pcem_grad_gfn1

!! ========================================================================
!  point charge embedding gradient for GFN2
!! ========================================================================
subroutine pcem_grad_gfn2(jData,g,gpc,n,pcem,at,nshell,xyz,qsh)
   use xtb_type_pcem
   implicit none
   type(TCoulombData), intent(in) :: jData
   integer, intent(in) :: n
   type(tb_pcem),intent(in) :: pcem
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,pcem%n)
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell(:)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: qsh(:)

   integer  :: ish,ii,j,iat,ati
   real(wp) :: dr(3)
   real(wp) :: gi,gj,r2,rr,yy

   ii = 0
   do iat = 1, n
      ati = at(iat)
      do ish = 1, nshell(ati)
         gi = jData%shellHardness(ish,ati)
         do j = 1, pcem%n
            dr = xyz(:,iat) - pcem%xyz(:,j)
            r2 = sum(dr**2)
            gj = pcem%gam(j)
            rr = 0.5_wp*(gi+gj)
            rr = 1.0_wp/rr**2
            yy = qsh(ii+ish)*pcem%q(j)*(r2+rr)**(-1.5_wp)
            g(:,iat) = g(:,iat) - dr*yy
            gpc(:,j) = gpc(:,j) + dr*yy
         end do
      end do
      ii = ii + nshell(ati)
   end do

end subroutine pcem_grad_gfn2


end module xtb_embedding
