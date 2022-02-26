! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Supported file types in this program
module xtb_mctc_filetypes
   use mctc_io_filetype, only : filetype, getFileType => get_filetype
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_chartools, only : toLowercase
   implicit none
   private

   public :: fileType, getFileType, generateFileMetaInfo, generateFileName
   public :: hessType


   type :: hessEnum
      integer :: tmol = 1
      integer :: orca = 2
      integer :: dftbplus = 3
   end type hessEnum
   type(hessEnum), parameter :: hessType = hessEnum()


   !> Default file type
   integer, parameter :: defaultFileType = fileType%tmol


contains


!> Split file name into path, basename and extension
subroutine generateFileMetaInfo(Name, path, basename, extension)

   !> File name
   character(len=*), intent(in) :: name

   !> Path to the file, empty if no path found in name
   character(len=:), allocatable, intent(out) :: path

   !> Base name of the file
   character(len=:), allocatable, intent(out) :: basename

   !> Extension of the file name, empty if no extension present
   character(len=:), allocatable, intent(out) :: extension

   integer :: iDot, iSlash

   iDot = index(name, '.', back=.true.)
   iSlash = index(name, '/', back=.true.)

   ! When given a name like "Gau-3333.Ein", the basename should be "Gau-3333",
   ! so if there is no '/', we assume there is an implicit './' at the
   ! beginning of name.
   if (iSlash > 0) then
      path = name(:iSlash)
   else
      path = ''
   endif

   if (iDot > iSlash .and. iDot > 0) then
      extension = name(iDot+1:)
   else ! means point is somewhere in the path or absent
      iDot = len(name)+1
      extension = ''
   endif

   if (iDot > iSlash) then
      basename = name(iSlash+1:iDot-1)
   else
      basename = ''
   endif


end subroutine generateFileMetaInfo


subroutine generateFileName(fname, basename, extension, ftype)
   character(len=:), allocatable, intent(out) :: fname
   character(len=*), intent(in) :: basename
   character(len=*), intent(in) :: extension
   integer, intent(in) :: ftype

   if (len(basename) > 0) then
      fname = basename
   else
      fname = 'xtb'
   endif

   if (len(extension) > 0) then
      fname = fname//'.'//extension
   else
      select case(ftype)
      case(fileType%xyz)
         fname = fname//'.xyz'
      case(fileType%tmol)
         fname = fname//'.coord'
      case(fileType%molfile)
         fname = fname//'.mol'
      case(fileType%sdf)
         fname = fname//'.sdf'
      case(fileType%vasp)
         fname = fname//'.poscar'
      case(fileType%pdb)
         fname = fname//'.pdb'
      case(fileType%gen)
         fname = fname//'.gen'
      case(fileType%gaussian)
         fname = fname//'.ein'
      end select
   endif
end subroutine generateFileName


end module xtb_mctc_filetypes
