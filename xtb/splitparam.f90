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

module splitparam
   use iso_fortran_env, wp => real64
   implicit none
   private :: wp
   public

!   integer  :: splitlist,iatf1,iatf2
!   real(wp) :: rcma,massf1,massf2,atmass
!   dimension :: splitlist(10000),atmass(10000)
   integer  :: maxsplit = 0
   integer  :: iatf1 = 0
   integer  :: iatf2 = 0
   integer, allocatable :: splitlist(:)
   integer, allocatable :: iatf(:)
   real(wp) :: rcma = 0.0_wp
   real(wp) :: massf1 = 0.0_wp
   real(wp) :: massf2 = 0.0_wp
   real(wp),allocatable :: atmass(:)

   integer  :: maxfrag = 0

contains

subroutine init_split(nsplit)
   implicit none
   integer,intent(in) :: nsplit
   maxsplit = nsplit
   call clear_split
   allocate( splitlist(nsplit), source = 0 )
   allocate( iatf(nsplit),      source = 0 )
   allocate( atmass(nsplit),    source = 0.0_wp )
end subroutine init_split

subroutine clear_split
   if (allocated(splitlist)) deallocate(splitlist)
   if (allocated(atmass))    deallocate(atmass)
end subroutine clear_split

subroutine splitm(nat,at,xyz,cn)
   use mctc_econv, only : autoaa
   use aoparam
   implicit none
   integer, intent(in) :: nat
   real(wp),intent(in) :: xyz(3,nat)
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: cn(nat)

   integer  :: i,j,k
   real(wp) :: r,f,rco
   real(wp),allocatable :: bond(:,:)
   allocate(bond(nat,nat))

!  determine covalent neighbours
   bond = 0
   do i=1,nat
      do j=1,nat
         r=norm2(xyz(:,i)-xyz(:,j))
         rco=rad(at(i))+rad(at(j))
         if(r.lt.2.5*rco) bond(j,i)=1
      enddo
      bond(i,i)=0
   enddo

   call mrec(i,xyz,cn,bond,nat,at,splitlist)

   iatf1=0
   iatf2=0
   do i=1,nat
      if(splitlist(i).gt.2) splitlist(i)=2
      if(splitlist(i).eq.1) iatf1=iatf1+1
      if(splitlist(i).eq.2) iatf2=iatf2+1
   enddo

end subroutine splitm

subroutine cmafrag(nat,at,xyz,r1,r2)
   use aoparam
   implicit none
   real(wp) xyz(3,nat),r1(3),r2(3)
   integer nat,at(nat)

   integer i
   real(wp) sum1x,atmas,sum1y,sum1z,sum2x,sum2y,sum2z,dr(3)

   sum1x=0._wp
   sum1y=0._wp
   sum1z=0._wp
   sum2x=0._wp
   sum2y=0._wp
   sum2z=0._wp
   do i=1,nat
      atmas=atmass(i)
      if(splitlist(i).eq.1)then
         sum1x=sum1x+atmas*xyz(1,i)
         sum1y=sum1y+atmas*xyz(2,i)
         sum1z=sum1z+atmas*xyz(3,i)
      else
         sum2x=sum2x+atmas*xyz(1,i)
         sum2y=sum2y+atmas*xyz(2,i)
         sum2z=sum2z+atmas*xyz(3,i)
      endif
   enddo
   !
   r1(1)=sum1x/massf1
   r1(2)=sum1y/massf1
   r1(3)=sum1z/massf1
   r2(1)=sum2x/massf2
   r2(2)=sum2y/massf2
   r2(3)=sum2z/massf2

   dr(1:3)=r1(1:3)-r2(1:3)

   rcma=sqrt(sum(dr*dr))

end subroutine cmafrag

subroutine splitprint(nat,at,xyz)
   use scanparam
   implicit none
   real(wp) xyz(3,nat)
   integer nat,at(nat)

   real(wp) ra(3),rb(3)
   integer i

   if(iatf1.eq.0.or.iatf2.eq.0) return

   massf1=0
   massf2=0
   do i=1,nat
      if(splitlist(i).eq.1)then
         massf1=massf1+atmass(i)
      else
         massf2=massf2+atmass(i)
      endif
   enddo
   call cmafrag(nat,at,xyz,ra,rb)

   write(output_unit,'(a)')
   write(output_unit,'(''molecular fragmentation (1/2 indicates fragments):'')')
   write(output_unit,'(72i1)') splitlist(1:nat)
   write(output_unit,'(''# atoms in fragment 1/2:'',2i6)')iatf1,iatf2
   write(output_unit,'('' fragment masses (1/2) :'',2f12.2)')massf1,massf2
   write(output_unit,'(''CMA distance (Bohr)    :'',f8.3)')rcma
   write(output_unit,'(''constraining FC (au)   :'',f8.4)')fcconstr

end subroutine splitprint

subroutine cmaiface(nat,at,xyz)
   implicit none
   real(wp) xyz(3,nat)
   integer nat,at(nat)

   integer i,j
   real(wp) r,r0,rvdw(94),asum,bsum,xsum,ff
   rvdw(1:94)= (/ &
    1.09155,0.86735,1.7478 ,1.5491 ,1.608  ,1.45515,1.31125,1.24085, &
    1.1498 ,1.0687 ,1.8541 ,1.74195,2.0053 ,1.89585,1.75085,1.65535, &
    1.5523 ,1.4574 ,2.12055,2.05175,1.94515,1.8821 ,1.86055,1.7207, &
    1.7731 ,1.72105,1.71635,1.6731 ,1.6504 ,1.61545,1.97895,1.93095, &
    1.83125,1.7634 ,1.6831 ,1.6048 ,2.3088 ,2.2382 ,2.1098 ,2.02985, &
    1.9298 ,1.87715,1.7845 ,1.73115,1.69875,1.67625,1.6654 ,1.731, &
    2.13115,2.0937 ,2.0075 ,1.94505,1.869  ,1.79445,2.52835,2.5907, &
    2.31305,2.31005,2.2851 ,2.26355,2.2448 ,2.22575,2.2117 ,2.06215, &
    2.12135,2.07705,2.1397 ,2.1225 ,2.1104 ,2.0993 ,2.0065 ,2.1225, &
    2.049  ,1.99275,1.94775,1.8745 ,1.7228 ,1.67625,1.6282 ,1.67995, &
    2.15635,2.1382 ,2.05875,2.0027 ,1.9322 ,1.8608 ,2.5398 ,2.4647, &
    2.35215,2.2126 ,2.2297 ,2.19785,2.17695,2.21705/)

   atmass=0
   ff=-5._wp

   asum=0
   do i=1,nat
      if(splitlist(i).eq.1)then
         xsum=0
         do j=1,nat
            if(splitlist(j).ne.1)then
               r=sqrt((xyz(1,i)-xyz(1,j))**2+ &
               &      (xyz(2,i)-xyz(2,j))**2+ &
               &      (xyz(3,i)-xyz(3,j))**2)
               r0=rvdw(at(i))+rvdw(at(j))
               xsum=xsum+exp(ff*r/r0**2)
            endif
         enddo
         asum=asum+xsum
         atmass(i)=xsum
      endif
   enddo
   do i=1,nat
      if(splitlist(i).eq.1)atmass(i)=atmass(i)/asum
   enddo

   bsum=0
   do i=1,nat
      if(splitlist(i).ne.1)then
         xsum=0
         do j=1,nat
            if(splitlist(j).eq.1)then
               r=sqrt((xyz(1,i)-xyz(1,j))**2+ &
               &      (xyz(2,i)-xyz(2,j))**2+ &
               &      (xyz(3,i)-xyz(3,j))**2)
               r0=rvdw(at(i))+rvdw(at(j))
               xsum=xsum+exp(ff*r/r0**2)
            endif
         enddo
         bsum=bsum+xsum
         atmass(i)=xsum
      endif
   enddo
   do i=1,nat
      if(splitlist(i).ne.1)atmass(i)=atmass(i)/bsum
   enddo
   ff=maxval(atmass(1:nat))
   atmass(1:nat)=atmass(1:nat)/ff

   !     do i=1,nat
   !        write(*,*) i,atmass(i)
   !     enddo

end subroutine cmaiface

end module splitparam
