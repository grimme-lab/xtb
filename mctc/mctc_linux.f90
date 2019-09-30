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

!! ------------------------------------------------------------------------
!  this module implements system calls (calls to functions implemented by
!  the Linux kernel, *not* call system()!) in a Fortran accessable way.
!  All system calls can be looked up by `man 2 <name>'.
!! ------------------------------------------------------------------------
module mctc_linux
   !  sys/types.h
   !  > typedef unsigned short mode_t; /* file type and permissions bits */
   !    Note: rwxr-xr-x -> o'0755' -> b'0000000111101101'
   use iso_c_binding, mode_t => c_short
   interface c_mkdir
      !  int mkdir(const char *pathname, mode_t mode);
      integer(c_int) function mkdir(pathname, mode) bind(c)
         import c_int,mode_t,c_char
         character(kind=c_char),intent(in) :: pathname
         integer(mode_t),value             :: mode
      end function mkdir
   end interface c_mkdir
   interface c_rmdir
      !  int rmdir(const char *pathname);
      integer(c_int) function rmdir(pathname) bind(c)
         import c_int,c_char
         character(kind=c_char),intent(in) :: pathname
      end function rmdir
   end interface c_rmdir
   interface c_link
      !  int link(const char *oldpath, const char *newpath);
      integer(c_int) function link(oldpath,newpath) bind(c)
         import c_int,c_char
         character(kind=c_char),intent(in) :: oldpath,newpath
      end function link
   end interface c_link
   interface c_unlink
      !  int unlink(const char *pathname);
      integer(c_int) function unlink(pathname) bind(c)
         import c_int,c_char
         character(kind=c_char),intent(in) :: pathname
      end function unlink
   end interface c_unlink
   interface c_symlink
      !  int symlink(const char *oldpath, const char *newpath);
      integer(c_int) function symlink(oldpath,newpath) bind(c)
         import c_int,c_char
         character(kind=c_char),intent(in) :: oldpath,newpath
      end function symlink
   end interface c_symlink
   interface c_chdir
      !  int chdir(const char *pathname);
      integer(c_int) function chdir(pathname) bind(c)
         import c_int,c_char
         character(kind=c_char),intent(in) :: pathname
      end function chdir
   end interface c_chdir
   interface c_umask
      !  mode_t umask(mode_t mask);
      integer(mode_t) function umask(pathname, mode) bind(c)
         import mode_t
         integer(mode_t),value :: mode
      end function umask
   end interface c_umask
   interface c_chmod
      !  int chmod(const char *pathname, mode_t mode);
      integer(c_int) function chmod(pathname, mode) bind(c)
         import c_int,mode_t,c_char
         character(kind=c_char),intent(in) :: pathname
         integer(mode_t),value             :: mode
      end function chmod
   end interface c_chmod

contains
!! ------------------------------------------------------------------------
function sys_mkdir(pathname) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)                :: pathname
   character(kind=c_char,len=len(pathname)+1) :: c_name
   integer(mode_t),parameter :: umask = int(O'0755',mode_t)
   c_name = pathname // c_null_char

   err = c_mkdir(c_name,umask)

end function sys_mkdir
!! ------------------------------------------------------------------------
function sys_rmdir(pathname) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)                :: pathname
   character(kind=c_char,len=len(pathname)+1) :: c_name
   c_name = pathname // c_null_char

   err = c_rmdir(c_name)

end function sys_rmdir
!! ------------------------------------------------------------------------
function sys_link(oldpath,newpath) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)               :: oldpath
   character(len=*),intent(in)               :: newpath
   character(kind=c_char,len=len(oldpath)+1) :: c_old
   character(kind=c_char,len=len(newpath)+1) :: c_new
   c_old = oldpath // c_null_char
   c_new = newpath // c_null_char

   err = c_link(c_old,c_new)
end function sys_link
!! ------------------------------------------------------------------------
function sys_unlink(pathname) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)                :: pathname
   character(kind=c_char,len=len(pathname)+1) :: c_name
   c_name = pathname // c_null_char

   err = c_unlink(c_name)

end function sys_unlink
!! ------------------------------------------------------------------------
function sys_symlink(oldpath,newpath) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)               :: oldpath
   character(len=*),intent(in)               :: newpath
   character(kind=c_char,len=len(oldpath)+1) :: c_old
   character(kind=c_char,len=len(newpath)+1) :: c_new
   c_old = oldpath // c_null_char
   c_new = newpath // c_null_char

   err = c_symlink(c_old,c_new)
end function sys_symlink
!! ------------------------------------------------------------------------
function sys_chdir(pathname) result(err)
   use iso_c_binding
   integer :: err
   character(len=*),intent(in)                :: pathname
   character(kind=c_char,len=len(pathname)+1) :: c_name
   c_name = pathname // c_null_char
   err = c_chdir(c_name)
end function sys_chdir

end module mctc_linux
