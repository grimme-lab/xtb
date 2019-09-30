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

module write_geometry

interface wrcoord
   module procedure :: wrcoord_res
end interface wrcoord

contains
! ------------------------------------------------------------------------
!> general driver for printout for geometry to any unit
!  this driver respects the geometry preference set in the input section
!  for specific geometry formats use the respective implementations below
subroutine wrcoord_res(iunit,n,at,xyz,sccres,frqres)
   use iso_fortran_env, only : output_unit, wp => real64
   use tbdef_data
   use setparam
   implicit none              
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   type(scc_results), intent(in) :: sccres
   type(freq_results),intent(in) :: frqres

   select case(geometry_inputfile)
   case(p_geo_xmol)
      call write_xyzlog(iunit,n,at,xyz,sccres%e_total,sccres%gnorm)
   case(p_geo_sdf)
      call write_sdf(iunit,n,at,xyz,sccres,frqres)
   case default
      call write_coord(iunit,n,at,xyz)
   end select

end subroutine wrcoord_res

! ------------------------------------------------------------------------
subroutine write_poscar(iunit,n,at,xyz,lat,e_total,gnorm)
   use iso_fortran_env, wp => real64
   use mctc_econv
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: lat(3,3)
   real(wp),intent(in),optional :: e_total
   real(wp),intent(in),optional :: gnorm
   logical :: v5
   integer :: i,j,iat
   integer,allocatable :: kinds(:)
   character(len=2),external :: asym

   allocate(kinds(n), source = 1)

   ! use vasp 5.x format if energy and gradient should be printed
   v5 = present(e_total) .and. present(gnorm)
   if (v5) then
      write(iunit,'(1x,"SCF done",2f24.12)') e_total, gnorm
   else
      j = 0
      iat = 0
      do i = 1, n
         if (iat.eq.at(i)) then
            kinds(j) = kinds(j)+1
         else
            j = j+1
            iat = at(i)
            write(iunit,'(1x,a)',advance='no') asym(iat)
         endif
      enddo
      write(iunit,'(a)')
   endif

   ! scaling factor for lattice parameters is always one
   write(iunit,'(f20.14)') 1.0_wp
   ! write the lattice parameters
   do i = 1, 3
      write(iunit,'(3f20.14)') lat(:,i)*autoaa
   enddo

   if (v5) then
      ! now we assume to write in vasp 5.x format, flush the symbols
      j = 0
      iat = 0
      do i = 1, n
         if (iat.eq.at(i)) then
            kinds(j) = kinds(j)+1
         else
            j = j+1
            iat = at(i)
            write(iunit,'(1x,a)',advance='no') asym(iat)
         endif
      enddo
      write(iunit,'(a)')
   endif

   ! write the count of the consequtive atom types
   do i = 1, j
      write(iunit,'(1x,i0)',advance='no') kinds(i)
   enddo
   write(iunit,'(a)')
   deallocate(kinds)

   ! we write cartesian coordinates
   write(iunit,'("Cartesian")')

   ! now write the cartesian coordinates
   do i = 1, n
      write(iunit,'(3f20.14)') xyz(:,i)*autoaa
   enddo

end subroutine write_poscar

! ------------------------------------------------------------------------
!> printout geometry as structure data file using the implementation
!  of the MDL Molfile V2000 (this sounds so awkward...) below
subroutine write_sdf(iunit,n,at,xyz,sccres,frqres)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use tbdef_data
   use gbobc
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   integer :: nb
   integer,allocatable :: bond(:,:)
   type(scc_results), intent(in) :: sccres
   type(freq_results),intent(in) :: frqres
   character(len=*),parameter :: sdfmt = '("> <",a,">",/,f18.8,/)'

   allocate(bond(n,n),source=0)
   call mdl_bond(n,at,nb,bond)

   if ((nb.gt.999) .or. (n.gt.999)) then
      call raise('S',"Cannot write structure data file, dimension to large!",1)
      return
   endif

   call write_mdl2000(iunit,n,at,xyz,nb,bond)
   write(iunit,sdfmt) "total energy in Eh", sccres%e_total
   if (lgbsa) &
   write(iunit,sdfmt) "GBSA Gsolv in Eh", sccres%g_solv
   write(iunit,sdfmt) "HL gap in eV", sccres%hl_gap
   write(iunit,sdfmt) "dipole moment in D", norm2(sccres%dipole)*autod
   if (frqres%gtot .ne. 0.0_wp) &
   write(iunit,sdfmt) "total free energy in Eh", frqres%gtot
   write(iunit,'("$$$$")')
end subroutine write_sdf

! ------------------------------------------------------------------------
!> implementation of the MDL molfile V2000 printing the CT (connection table)
subroutine write_mdl2000(iunit,n,at,xyz,nb,bond)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use setparam
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   integer, intent(in) :: nb
   integer, intent(in) :: bond(n,n)
   integer :: i,j
   character(len=2),external :: asym
   integer,parameter :: z4(4) = [(0,i=1,4)], z12(12) = [(0,i=1,12)]

   if ((nb.gt.999) .or. (n.gt.999)) then
      call raise('S',"Cannot write structure data file, dimension to large!",1)
      return
   endif

   if (allocated(molnameline)) then
      write(iunit,'(a)') molnameline
   else
      write(iunit,'(a)')
   endif
   write(iunit,'("xtb program, MCTC, University of Bonn, Germany")')
   if (allocated(commentline)) then
      write(iunit,'(a)') commentline
   else
      write(iunit,'(a)')
   endif
   write(iunit,'(2i3,"  0     0  0  0  0  0  0999 V2000")') n,nb
   do i=1,n
      write(iunit,'(3F10.5,a3,12i3)') &
         xyz(:,i)*autoaa,asym(at(i)),z12
   enddo
   do i=1,n
      do j=1,i-1
         if (bond(j,i).gt.0) &
         write(iunit,'(7i3)') j,i,bond(j,i),z4
      enddo
   enddo
   if (ichrg .ne. 0) &
   write(iunit,'("M  CRG",1x,i0)') ichrg
   write(iunit,'("M  END")') 

end subroutine write_mdl2000

! ------------------------------------------------------------------------
!> implementation of the MDL molfile V3000 printing the CT (connection table)
subroutine write_mdl3000(iunit,n,at,xyz,nb,bond)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use setparam
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   integer, intent(in) :: nb
   integer, intent(in) :: bond(n,n)
   integer :: i,j,k
   character(len=2),external :: asym
   character(len=*),parameter :: v30="M  V30"
   integer,parameter :: z3(3) = [(0,i=1,3)]

   if (allocated(molnameline)) then
      write(iunit,'(a)') molnameline
   else
      write(iunit,'(a)')
   endif
   write(iunit,'("xtb program, MCTC, University of Bonn, Germany")')
   if (allocated(commentline)) then
      write(iunit,'(a)') commentline
   else
      write(iunit,'(a)')
   endif
   write(iunit,'("  0  0  0     0  0            999 V3000")')
   write(iunit,'(a,1x,"BEGIN CTAB")') v30
   write(iunit,'(a,1x,"COUNTS",5(1x,i0))') v30,n,nb,z3
   write(iunit,'(a,1x,"BEGIN ATOM")') v30
   do i=1,n
      write(iunit,'(a,1x,i0,1x,a,3F10.5,1x,i0)') &
         v30,i,asym(at(i)),xyz(:,i)*autoaa,0
   enddo
   if (ichrg.ne.0) & ! no separate CHG anymore -> "star"-atom for charge...
   write(iunit,'(a,1x,i0,1x,a,3F10.5,1x,i0,1x,"CHG=",i0)') &
      v30,n+1,'*',0.0_wp,0.0_wp,0.0_wp,0,ichrg
   write(iunit,'(a,1x,"END ATOM")') v30
   if (nb.gt.0) then
   write(iunit,'(a,1x,"BEGIN BOND")') v30
   k = 0
   do i=1,n
      do j=1,i-1
         if (bond(j,i).eq.0) cycle
         k = k+1
         write(iunit,'(a,4(1x,i0))') v30,k,bond(j,i),j,i
      enddo
   enddo
   write(iunit,'(a,1x,"END BOND")') v30
   endif
   write(iunit,'(a,1x,"END CTAB")') v30
   write(iunit,'("M  END")') 

end subroutine write_mdl3000

! ------------------------------------------------------------------------
!> print geometry as Turbomole coord data group
subroutine write_coord(iunit,n,at,xyz)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   integer :: i
   character(len=2),external :: asym

   write(iunit,'(a)') "$coord"
   do i = 1, n
      write(iunit,'(3f20.14,6x,a2)') xyz(:,i), asym(at(i))
   enddo
   write(iunit,'(a)') "$end"

end subroutine write_coord

! ------------------------------------------------------------------------
!> coordinate output in xmol/xyz format with additional information in
!  the commentline obtained from an SCC calculation
subroutine write_xyzlog(iunit,n,at,xyz,e_total,gnorm)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use tbdef_data
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: e_total
   real(wp),intent(in) :: gnorm
   integer :: i
   character(len=2),external :: asym

   write(iunit,'(i0)') n
   write(iunit,'(1x,"SCF done",2f24.12)') e_total, gnorm
   do i = 1, n
      write(iunit,'(a2,6x,3f20.14)') asym(at(i)), xyz(:,i)*autoaa
   enddo

end subroutine write_xyzlog

! ------------------------------------------------------------------------
!> plain xyz-printout, retains commentline from input file
subroutine write_xyz(iunit,n,at,xyz)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use setparam
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   integer :: i
   character(len=2),external :: asym

   write(iunit,'(i0)') n
   if (allocated(commentline)) then
      write(iunit,'(a)') commentline
   else
      write(iunit,'(a)')
   endif
   do i = 1, n
      write(iunit,'(a2,6x,3f20.14)') asym(at(i)), xyz(:,i)*autoaa
   enddo

end subroutine write_xyz

! ------------------------------------------------------------------------
!> bond analysis using wbo file, this file might or might not already be
!  present from the calculation, so this can produce complete bogus
subroutine mdl_bond(n,at,nb,bond)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(out) :: nb
   integer, intent(out) :: bond(n,n)
   integer  :: i,j,bo,err
   real(wp) :: w
   integer  :: ich

   bond = 0
   nb=0
   call open_file(ich,'wbo','r')
   if (ich.eq.-1) return
   do
      read(ich,*,iostat=err) i,j,w 
      if (err.ne.0) exit
      if (w.lt.0.1_wp) cycle
      nb = nb + 1
      bo = 1
      if (w.gt.1.3_wp) bo=2
      if (w.gt.2.3_wp) bo=3
      if (w.gt.1.2_wp.and.w.lt.1.5_wp &
         .and.at(i).eq.6.and.at(j).eq.6) bo=4
      bond(i,j)=bo
      bond(j,i)=bo
   enddo
   call close_file(ich)
end subroutine mdl_bond

end module write_geometry
