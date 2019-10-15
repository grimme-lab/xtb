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

subroutine wrgbw(n,iat,coord,z,basis,wfn)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding
   use tbdef_basisset
   use tbdef_wavefunction
   use aoparam
   implicit none
!! ------------------------------------------------------------------------
!  xtb input
   integer, intent(in) :: n
   integer, intent(in) :: iat(n)
   real(wp),intent(in) :: coord(3,n)
   real(wp),intent(in) :: z(n)
   type(tb_basisset),    intent(in) :: basis
   type(tb_wavefunction),intent(in) :: wfn
!! ------------------------------------------------------------------------
!  ORCA input
!  file pointer
   integer(c_long) :: ptr_info
   integer(c_long) :: ptr_geometry
   integer(c_long) :: ptr_basisset
   integer(c_long) :: ptr_wavefunction
   integer(c_long),parameter :: ptr_ecp = 0
!  Geometry
   integer(c_int) :: nat
   integer(c_int),allocatable :: at(:)
   integer(c_int),allocatable :: frag(:)
   integer(c_int),allocatable :: fix(:,:)
   real(c_double),allocatable :: xyz(:,:)
   real(c_double),allocatable :: za(:)
   real(c_double),allocatable :: xi(:)
   real(c_double),allocatable :: am(:)

!  Basisset
   integer(c_int),parameter :: maxgauss = 32
   type :: BFNGauss ! defined exactly as the ORCA equivalent
      integer(c_int) :: l = 0       ! angular momentum quantum number for the shell
      integer(c_int) :: lstart = 0  ! where the shell starts in the basis
      integer(c_int) :: ng = 0      ! Number of components
      real(c_double) :: a(maxgauss) = 0.0_c_double ! orbital exponents
      real(c_double) :: d(maxgauss) = 0.0_c_double ! expansion coeffients
   end type BFNGauss
   logical(c_bool),parameter :: Have_STO     = .false.
   logical(c_bool),parameter :: Have_GTO     = .true.
   logical(c_bool),parameter :: Have_GTOAUX  = .false.
   logical(c_bool),parameter :: Have_GTOCABS = .false.
   logical(c_bool),parameter :: Have_Label   = .false.
   integer(c_int),allocatable :: nshell(:)
   type(BFNGauss),allocatable :: BG(:,:)

!  Wavefunction
   integer(c_int),parameter :: nop = 1
   integer(c_int) :: nao
   real(c_double),allocatable :: C(:,:)
   real(c_double),allocatable :: occ(:)
   real(c_double),allocatable :: emo(:)
   integer(c_int),allocatable :: irrep(:)
   integer(c_int),allocatable :: core(:)

!  local
   integer :: iunit
   logical :: exist
   integer :: i,ish,ati,ishtyp,icao,iao,ip,iprim,icount
   character(kind=c_char),allocatable :: info(:)
   integer(int64) :: iend,istart,ival(3)
   integer(int64) :: total_size = 0
   integer(int64) :: size_info = 0
   integer(int64) :: size_pointer = 0
   integer(int64) :: size_geometry = 0
   integer(int64) :: size_basisset = 0
   integer(int64) :: size_wavefunction = 0
   integer(c_int) :: ishell

!  check if we can provide the necessary C-types
   if (c_long.eq.0) &
      call raise('E','Fortran and C kinds are incompatible!',1)
   if (c_int.eq.0 .or. c_int.ne.kind(n)) &
      call raise('E','Fortran and C kinds are incompatible!',1)
   if (c_double.eq.0 .or. c_double.ne.wp) &
      call raise('E','Fortran and C kinds are incompatible!',1)

!  copy the important dimensions
   nat = int(n,c_int)
   nao = int(wfn%nao,c_int)
!  and allocate some memory for the printout
!  geometry
   allocate( xyz(3,nat) ); xyz = real(coord,c_double)
   allocate( za(nat) );     za = real(z,c_double)
   allocate( xi(nat) );     xi = real(0.0_wp,c_double)
   allocate( am(nat) );     am = real(0.0_wp,c_double)
   allocate( at(nat) );     at = int(iat,c_int)
   allocate( frag(nat) ); frag = int(1,c_int)
   allocate( fix(3,nat) ); fix = int(0,c_int)
!  basis set
   allocate( BG(3,nat) )
   allocate( nshell(nat)); nshell = int(0,c_int)
!  wavefunction
   allocate( C(nao,nao) ); C = real(wfn%C,c_double)
   allocate( emo(nao) ); emo = real(wfn%emo,c_double)
   allocate( occ(nao) ); occ = real(wfn%focc,c_double)
   allocate( irrep(nao)); irrep = int(0,c_int)
   allocate( core(nao));  core  = int(0,c_int)

!  BFNGauss** is a ragged array of BFNGauss[iat][ishell]
!  but we don't use ragged pointer arrays, so we will fake is...
   do i = 1, nat
      ati = iat(i)
      nshell(i) = int(ao_n(ati),c_int)
      do ish = 1, ao_n(ati)
         ishtyp = ao_l(ish,ati)
         icao = basis%caoshell(ish,ati)
         iao  = basis%saoshell(ish,ati)
         BG(ish,i)%lstart = int(iao,c_int)
         BG(ish,i)%l      = int(ishtyp,c_int)
         BG(ish,i)%ng     = int(basis%nprim(icao+1),c_int)
         do ip = 1, basis%nprim(icao+1)
            iprim = ip + basis%primcount(icao+1)
            BG(ish,i)%a(ip) = real(basis%alp(iprim),c_double)
            BG(ish,i)%d(ip) = real(basis%cont(iprim),c_double)
         enddo
      enddo
   enddo

!  the file pointer section are 5 words of long
   size_pointer = 5*8
!  the geometry section are one + 5*nat words of int and 6*nat words of double precision
   size_geometry = (1+5*n)*4 + (6*n)*8
!  the basis set section are one word of int, 5 words of bool, nat words of int
!  and nshell times 4 words of int and 64 words of double precision
   size_basisset = (1+n)*4 + 5 + (64*8 + 4*4)*wfn%nshell
!  wavefunction
   size_wavefunction = (1+1+wfn%nao+wfn%nao)*4 + (wfn%nao**2+wfn%nao+wfn%nao)*8

   inquire(file='orca.gbw',exist=exist)
   if (.not.exist) &
      call raise('E','Unfortunately, we need an template gbw-file write the info block, otherwise this is not working.',1)
   open(newunit=iunit,file='orca.gbw',form='unformatted',access='stream')
   read(iunit) istart
   read(iunit) ival
   iend = minval(ival)
   size_info = iend-istart
   print*,istart,ival,iend,size_info
   allocate( info(size_info) )
   read(iunit,pos=istart) info
   close(iunit)

   ptr_info = int(size_pointer,c_long)
   ptr_geometry = size_info + ptr_info
   ptr_wavefunction = ptr_geometry + int(size_geometry,c_long)
   ptr_basisset = ptr_wavefunction + int(size_wavefunction,c_long)

   total_size = size_pointer+size_geometry+size_basisset+size_wavefunction+size_info
   write(istdout,'(/,"writing xtb.gbw (",i0,1x,"bytes)")') total_size
   write(istdout,'("file pointer I at",1x,z16)') ptr_info
   write(istdout,'("file pointer G at",1x,z16)') ptr_geometry
   write(istdout,'("file pointer B at",1x,z16)') ptr_basisset
   write(istdout,'("file pointer W at",1x,z16)') ptr_wavefunction

!  from here on this looks like ORCA printout, it behaves like ORCA printout,
!  but it is *not* ORCA printout (but ORCA will never get the difference).

   open(newunit=iunit,file='xtb.gbw',form='unformatted',access='stream')

!  write file pointer
   write(iunit) ptr_info          ! 1× int64
   write(iunit) ptr_geometry      ! 1× int64
   write(iunit) ptr_basisset      ! 1× int64
   write(iunit) ptr_wavefunction  ! 1× int64
   write(iunit) ptr_ecp           ! 1× int64
   !write(istdout,'(z16,1x,z16)') ptr_info,ptr_geometry
   !write(istdout,'(z16,1x,z16)') ptr_basisset,ptr_wavefunction
   !write(istdout,'(z16)') ptr_ecp

   write(iunit) info

!  write geometry
   write(iunit) nat          ! 1×  int32
   !write(istdout,'(z8)') nat
   do i = 1, nat
      !write(istdout,'(z16,1x,z16,/,z16,1x,z16)') xyz(:,i),za(i)
      write(iunit) xyz(:,i)  ! 3× real64
      write(iunit) za(i)     ! 1× real64
      write(iunit) xi(i)     ! 1× real64
      write(iunit) am(i)     ! 1× real64
      write(iunit) at(i)     ! 1×  int32
      write(iunit) frag(i)   ! 1×  int32
      write(iunit) fix(:,i)  ! 3×  int32
   enddo

!  write wavefuntion
   write(iunit) nop          !  1× int32
   write(iunit) nao          !  1× int32
   write(iunit) C            !  nao*nao × real64
   write(iunit) occ          !  nao × real64
   write(iunit) emo          !  nao × real64
   write(iunit) irrep        !  nao ×  int32
   write(iunit) core         !  nao ×  int32

!  write basis set
   write(iunit) nat          ! 1× int32
   write(iunit) Have_STO     ! 1× bool
   write(iunit) Have_GTO     ! 1× bool
   write(iunit) Have_GTOAUX  ! 1× bool
   write(iunit) Have_GTOCABS ! 1× bool
   write(iunit) Have_Label   ! 1× bool

   write(iunit) nshell       ! nat × int32
   do i = 1, nat
      ishell = nshell(i)
      write(iunit) BG(1:ishell,i) ! ishell × (3× int32 + 64× real64)
   enddo

   close(iunit)

end subroutine wrgbw
