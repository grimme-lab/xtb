! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module embedding
   use iso_fortran_env, wp => real64
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

   interface read_pcem
      module procedure read_pcem_legacy
      module procedure read_pcem_orca
   end interface

   interface jpot_pcem_gfn1
      module procedure jpot_pcem_gfn1_new
      module procedure jpot_pcem_gfn1_old
   end interface

   interface jpot_pcem_gfn2
      module procedure jpot_pcem_gfn2_new
      module procedure jpot_pcem_gfn2_old
   end interface

   interface pcem_grad_gfn1
      module procedure pcem_grad_gfn1_new
      module procedure pcem_grad_gfn1_old
   end interface

   interface pcem_grad_gfn2
      module procedure pcem_grad_gfn2_new
      module procedure pcem_grad_gfn2_old
   end interface

contains

subroutine init_pcem
   use setparam
   implicit none
   select case(pcem_interface)
   case default
      if (.not.allocated(pcem_file)) pcem_file = "pcharge"
      if (.not.allocated(pcem_grad)) pcem_grad = "pcgrad"
   end select
end subroutine init_pcem

!! ========================================================================
!  initialize point charges
!! ========================================================================
subroutine read_pcem_legacy(fname,pcem,npc,pc,qpc,apc,pr)
   implicit none
   character(len=*),intent(in) :: fname
   logical, intent(out) :: pcem
   integer, intent(out) :: npc
   real(wp),allocatable,intent(out) :: pc(:,:)
   real(wp),allocatable,intent(out) :: qpc(:)
   integer, allocatable,intent(out) :: apc(:)
   logical, intent(in)  :: pr
   real(wp) :: pcxx(10),xsum
   integer  :: i,pcio,pcnn
   character(len=80) :: atmp
   npc=0
   call open_file(pcio,fname,'r')
   pcem = pcio.ne.-1
   if (pcem) then
      !write(*,*)'reading point charges ...'
      read(pcio,*) npc
      allocate(pc(3,npc),qpc(npc),apc(npc))
      do i=1,npc
         read(pcio,'(a)')atmp
         call readl(atmp,pcxx,pcnn)
         qpc(i)=pcxx(1)
         pc(1:3,i)=pcxx(2:4)
         if(pcnn.le.4)then
            apc(i)=7
         else
            apc(i)=int(pcxx(5))
         endif
      enddo
      call close_file(pcio)
      xsum=sum(qpc)
      if(pr) then
         write(output_unit,'('' # of PC   : '',i6   )')npc
         write(output_unit,'('' sum of qPC: '',f10.5)')xsum
      endif
   endif

end subroutine read_pcem_legacy

subroutine read_pcem_orca(iunit,pcem)
   use tbdef_pcem
   use mctc_strings
   use mctc_econv
   use aoparam
   use setparam
   use readin, only : getline => strip_line
   implicit none
   integer,      intent(in)    :: iunit
   type(tb_pcem),intent(inout) :: pcem

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
         call raise('E',"Not enough entries for PC, please check!",1)
      endif

      ! we count this line
      n = n+1
      ! before we start parsing we check the number of lines
      if (n > npc) then ! bad user
         call raise('E',"Wrong dimension input for PC, too many lines provided",1)
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
               call raise('E',"Invalid PC input: '"//trim(argv(5))//"'",1)
            endif
            gami = gam(iat)
         endif
         ! GFN0-xTB has negative gam-values since they are internally combined
         ! with the atom radius for the full chemical hardness (overall > 0).
         ! This is certainly a problem here, but since GFN0-xTB does
         ! not provide the possibility to use a PC potential, we don't care...
         ! But if YOU are reading this comment and are in charge of implementing
         ! a PC potential for GFN0-xTB, than this is YOUR problem now.
         if (gami < 0.0_wp) then
            call raise('S',"Found negative gam-value in user input: '"//&
               &            trim(argv(5))//"'",1)
         endif
      else
         ! we trust the dummy atom from xcontrol, since xcontrol could already
         ! have handled potential errors, which is usually not the case...
         gami = gam(pcem_dummyatom)
      endif
      pcem%xyz(:,n) = xyz*conv
      pcem%q(n) = q
      pcem%gam(n) = gami
   enddo

   ! user screwed up, we could rescue this by resetting npc
   ! or we rub it in his/her face :)
   if (n /= npc) then
      call raise('E',"Wrong dimension input for PC, too few lines provided",1)
   endif

end subroutine read_pcem_orca

!! ========================================================================
!  J potentials for GFN1 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn1_new(n,pcem,nshell,at,xyz,ash,lsh,alphaj,Vpc)
   use mctc_econv, only : autoev
   use tbdef_pcem
   use aoparam, only : lpar,gam
   implicit none
   integer, intent(in)  :: n
   type(tb_pcem),intent(inout) :: pcem
   integer, intent(in)  :: nshell
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   integer, intent(in)  :: ash(nshell)
   integer, intent(in)  :: lsh(nshell)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(inout) :: Vpc(nshell)
   integer  :: is,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,rab

   do is = 1, nshell
      iat = ash(is)
      ati = at(iat)
      gi  = gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      eh1 = 0.0_wp
      do kk = 1, pcem%n
         gj = pcem%gam(kk)
         rab = norm2(pcem%xyz(:,kk) - xyz(:,iat))
         xj = 2.0_wp/(1./gi+1./gj)
         dum = autoev/(rab**alphaj+1._wp/xj**alphaj)**(1._wp/alphaj)
         eh1 = eh1+pcem%q(kk)*dum
      enddo
      Vpc(is) = eh1
   enddo

end subroutine jpot_pcem_gfn1_new

!! ========================================================================
!  J potentials for GFN2 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn2_new(n,pcem,nshell,at,xyz,ash,lsh,Vpc)
   use mctc_econv, only : autoev
   use tbdef_pcem
   use aoparam, only : lpar,gam
   implicit none
   integer, intent(in)  :: n
   type(tb_pcem),intent(in) :: pcem
   integer, intent(in)  :: nshell
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   integer, intent(in)  :: ash(nshell)
   integer, intent(in)  :: lsh(nshell)
   real(wp),intent(inout) :: Vpc(nshell)
   integer  :: is,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,r2

   do is = 1, nshell
      iat = ash(is)
      ati = at(iat)
      gi  = gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      eh1 = 0.0_wp
      do kk = 1, pcem%n
         gj = pcem%gam(kk)
         r2 = sum((pcem%xyz(:,kk) - xyz(:,iat))**2)
         xj = 0.5_wp*(gi+gj)
         dum = autoev/sqrt(r2+1._wp/xj**2)
!        dum = autoev/sqrt(r2+1._wp/(gi*gj)) !NEWAV
         eh1 = eh1+pcem%q(kk)*dum
      enddo
      Vpc(is) = eh1
   enddo

end subroutine jpot_pcem_gfn2_new

!! ========================================================================
!  J potentials for GFN1 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn1_old(n,npc,nshell,at,apc,xyz,pc,qpc,ash,lsh,alphaj,Vpc)
   use mctc_econv, only : autoev
   use aoparam, only : lpar,gam
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: npc
   integer, intent(in)  :: nshell
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: apc(npc)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: pc(3,npc)
   real(wp),intent(in)  :: qpc(npc)
   integer, intent(in)  :: ash(nshell)
   integer, intent(in)  :: lsh(nshell)
   real(wp),intent(in)  :: alphaj
   real(wp),allocatable,intent(out) :: Vpc(:)
   integer  :: is,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,r2,rab,xa,ya,za

   if (allocated(Vpc)) deallocate(Vpc)
   allocate( Vpc(nshell), source = 0.0_wp )

   do is=1,nshell
      iat=ash(is)
      ati=at(iat)
      gi=gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      eh1=0.0_wp
      do kk=1,npc
         gj=gam(apc(kk))
         r2=(pc(1,kk)-xa)**2+(pc(2,kk)-ya)**2+(pc(3,kk)-za)**2
         rab=sqrt(r2)
         xj=2.0_wp/(1./gi+1./gj)
         dum=autoev/(rab**alphaj+1._wp/xj**alphaj)**(1._wp/alphaj)
         eh1=eh1+qpc(kk)*dum
      enddo
      Vpc(is)=eh1
   enddo

end subroutine jpot_pcem_gfn1_old

!! ========================================================================
!  J potentials for GFN2 including the point charge stuff
!! ========================================================================
subroutine jpot_pcem_gfn2_old(n,npc,nshell,at,apc,xyz,pc,qpc,ash,lsh,Vpc)
   use mctc_econv, only : autoev
   use aoparam, only : lpar,gam
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: npc
   integer, intent(in)  :: nshell
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: apc(npc)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: pc(3,npc)
   real(wp),intent(in)  :: qpc(npc)
   integer, intent(in)  :: ash(nshell)
   integer, intent(in)  :: lsh(nshell)
   real(wp),allocatable,intent(out) :: Vpc(:)
   integer  :: is,iat,ati,kk
   real(wp) :: eh1,dum,gi,gj,xj,r2,rab,xa,ya,za

   if (allocated(Vpc)) deallocate(Vpc)
   allocate( Vpc(nshell), source = 0.0_wp )

   do is=1,nshell
      iat=ash(is)
      ati=at(iat)
      gi=gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      eh1=0.0_wp
      do kk=1,npc
         gj=gam(apc(kk))
         r2=(pc(1,kk)-xa)**2+(pc(2,kk)-ya)**2+(pc(3,kk)-za)**2
         xj=0.5_wp*(gi+gj)
         dum=autoev/sqrt(r2+1._wp/xj**2)
!        dum=autoev/sqrt(r2+1._wp/(gi*gj)) !NEWAV
         eh1=eh1+qpc(kk)*dum
      enddo
      Vpc(is)=eh1
   enddo

end subroutine jpot_pcem_gfn2_old

!! ========================================================================
!  point charge embedding ES subroutine for SCC
!! ========================================================================
pure subroutine electro_pcem(nshell,dqsh,Vpc,es,scc)
   use mctc_econv, only : evtoau
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

!  ES energy in Eh
   es = es*evtoau

!  Etotal
   scc = scc + es

end subroutine electro_pcem

!! ========================================================================
!  point charge embedding gradient for GFN1
!! ========================================================================
subroutine pcem_grad_gfn1_new(g,gpc,n,pcem,at,nshell,xyz,ash,lsh,alphaj,qsh)
   use aoparam, only : lpar,gam
   use tbdef_pcem
   implicit none
   integer, intent(in) :: n
   type(tb_pcem),intent(in) :: pcem
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,pcem%n)
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: alphaj
   real(wp),intent(in) :: qsh(nshell)

   integer  :: is,j,iat,ati
   real(wp) :: dr(3)
   real(wp) :: gi,gj,r2,rr,yy,ff


   do is = 1, nshell
      iat = ash(is)
      ati = at(iat)
      gi  = gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      do j = 1, pcem%n
         dr = xyz(:,iat) - pcem%xyz(:,j)
         r2 = sum(dr**2)
         gj = pcem%gam(j)
         rr = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)
         rr = 1.0_wp/rr**alphaj
         ff = r2**(alphaj/2.0_wp-1.0_wp)*&
            &(r2**(alphaj*0.5_wp)+rr)**(-1.0_wp/alphaj-1.0_wp)
         yy = ff*qsh(is)*pcem%q(j)
         g(:,iat)=g(:,iat)-dr*yy
         gpc(:,j)=gpc(:,j)+dr*yy
      enddo
   enddo

end subroutine pcem_grad_gfn1_new

!! ========================================================================
!  point charge embedding gradient for GFN2
!! ========================================================================
subroutine pcem_grad_gfn2_new(g,gpc,n,pcem,at,nshell,xyz,ash,lsh,qsh)
   use aoparam, only : lpar,gam
   use tbdef_pcem
   implicit none
   integer, intent(in) :: n
   type(tb_pcem),intent(in) :: pcem
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,pcem%n)
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: qsh(nshell)

   integer  :: is,j,iat,ati
   real(wp) :: dr(3)
   real(wp) :: gi,gj,r2,rr,yy

   do is = 1,nshell
      iat = ash(is)
      ati = at(iat)
      gi  = gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      do j = 1, pcem%n
         dr = xyz(:,iat) - pcem%xyz(:,j)
         r2 = sum(dr**2)
         gj = pcem%gam(j)
         rr = 0.5_wp*(gi+gj)
         rr = 1.0_wp/rr**2
!        rr = 1.0d0/(gi*gj) !NEWAV
         yy = qsh(is)*pcem%q(j)*(r2+rr)**(-1.5_wp)
         g(:,iat) = g(:,iat) - dr*yy
         gpc(:,j) = gpc(:,j) + dr*yy
      enddo
   enddo

end subroutine pcem_grad_gfn2_new

!! ========================================================================
!  point charge embedding gradient for GFN1
!! ========================================================================
subroutine pcem_grad_gfn1_old(g,gpc,n,npc,at,apc,nshell,xyz,pc,ash,lsh, &
   &                          alphaj,qsh,qpc)
   use aoparam, only : lpar,gam
   implicit none
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,npc)
   integer, intent(in) :: n
   integer, intent(in) :: npc
   integer, intent(in) :: at(n)
   integer, intent(in) :: apc(npc)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: pc(3,npc)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: alphaj
   real(wp),intent(in) :: qsh(nshell)
   real(wp),intent(in) :: qpc(npc)

   integer  :: is,j,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy,ff


   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=gam(ati)*(1.0d0+lpar(lsh(is),ati))
      do j=1,npc
         dx=xa-pc(1,j)
         dy=ya-pc(2,j)
         dz=za-pc(3,j)
         !atj=at(jat) ! out of bounds
         r2=((pc(1,j)-xyz(1,iat))**2 &
         &  +(pc(2,j)-xyz(2,iat))**2 &
         &  +(pc(3,j)-xyz(3,iat))**2)
         gj=gam(apc(j))
         rr=2.0d0/(1./gi+1./gj)
         rr=1.0d0/rr**alphaj
         ff=r2**(alphaj/2.0d0-1.0d0)*&
         & (r2**(alphaj*0.5d0)+rr)**(-1.0d0/alphaj-1.0d0)
         yy=ff*qsh(is)*qpc(j)
         g(1,iat)=g(1,iat)-dx*yy
         g(2,iat)=g(2,iat)-dy*yy
         g(3,iat)=g(3,iat)-dz*yy
         gpc(1,j)=gpc(1,j)+dx*yy
         gpc(2,j)=gpc(2,j)+dy*yy
         gpc(3,j)=gpc(3,j)+dz*yy
      enddo
   enddo

end subroutine pcem_grad_gfn1_old

!! ========================================================================
!  point charge embedding gradient for GFN2
!! ========================================================================
subroutine pcem_grad_gfn2_old(g,gpc,n,npc,at,apc,nshell,xyz,pc,ash,lsh,qsh,qpc)
   use aoparam, only : lpar,gam
   implicit none
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: gpc(3,npc)
   integer, intent(in) :: n
   integer, intent(in) :: npc
   integer, intent(in) :: at(n)
   integer, intent(in) :: apc(npc)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: pc(3,npc)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: qsh(nshell)
   real(wp),intent(in) :: qpc(npc)

   integer  :: is,j,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy


   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=gam(ati)*(1.0d0+lpar(lsh(is),ati))
      do j=1,npc
         dx=xa-pc(1,j)
         dy=ya-pc(2,j)
         dz=za-pc(3,j)
         !atj=at(jat) ! out of bounds
         r2=((pc(1,j)-xyz(1,iat))**2 &
         &         +(pc(2,j)-xyz(2,iat))**2 &
         &         +(pc(3,j)-xyz(3,iat))**2)
         gj=gam(apc(j))
         rr=0.5d0*(gi+gj)
         rr=1.0d0/rr**2
!        rr=1.0d0/(gi*gj) !NEWAV
         yy=qsh(is)*qpc(j)*(r2+rr)**(-1.5d0)
         g(1,iat)=g(1,iat)-dx*yy
         g(2,iat)=g(2,iat)-dy*yy
         g(3,iat)=g(3,iat)-dz*yy
         gpc(1,j)=gpc(1,j)+dx*yy
         gpc(2,j)=gpc(2,j)+dy*yy
         gpc(3,j)=gpc(3,j)+dz*yy
      enddo
   enddo

end subroutine pcem_grad_gfn2_old

end module embedding
