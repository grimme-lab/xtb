! This file is part of xtb.
!
! Copyright (C) 2021 Sebastian Ehlert
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

module xtb_freq_io
   use xtb_mctc_accuracy, only : wp
   use xtb_lin, only : lin
   implicit none
   private

   public :: writeHessianOut, rdhess, wrhess, write_tm_vibspectrum
   public :: g98fake, g98fake2


contains


!> Write the second derivative matrix
subroutine writeHessianOut(fname, hessian)

   !> File name
   character(len=*), intent(in) :: fname

   !> Dynamical (Hessian) matrix
   real(wp), intent(in) :: hessian(:,:)

   !> Format string for energy second derivative matrix
   character(len=*), parameter :: fmt = '(4f16.10)'

   integer :: ii, id

   call open_file(id, fname, 'w')
   do ii = 1, size(hessian, dim=2)
      write(id, fmt) hessian(:, ii)
   end do
   call close_file(id)

end subroutine writeHessianOut


subroutine wrhess(nat3,h,fname)
   integer, intent(in) :: nat3
   real(wp),intent(in) :: h(nat3*(nat3+1)/2)
   character(len=*),intent(in) :: fname
   integer iunit,i,j,mincol,maxcol,k
   character(len=5)  :: adum
   character(len=80) :: a80

   adum='   '
   call open_file(iunit,fname,'w')
   a80='$hessian'
   write(iunit,'(a)')a80
   do i=1,nat3
      maxcol = 0
      k=0
      200    mincol = maxcol + 1
      k=k+1
      maxcol = min(maxcol+5,nat3)
      write(iunit,'(a5,5f15.10)')adum,(h(lin(i,j)),j=mincol,maxcol)
      if (maxcol.lt.nat3) goto 200
   enddo
   call close_file(iunit)

end subroutine wrhess

subroutine rdhess(nat3,h,fname)
   integer, intent(in)  :: nat3
   real(wp),intent(out) :: h(nat3,nat3)
   character(len=*),intent(in) :: fname
   integer  :: iunit,i,j,mincol,maxcol
   character(len=5)  :: adum
   character(len=80) :: a80

   !     write(*,*) 'Reading Hessian <',trim(fname),'>'
   call open_file(iunit,fname,'r')
   50  read(iunit,'(a)')a80
   if(index(a80,'$hessian').ne.0)then
      do i=1,nat3
         maxcol = 0
         200       mincol = maxcol + 1
         maxcol = min(maxcol+5,nat3)
         read(iunit,*)(h(j,i),j=mincol,maxcol)
         if (maxcol.lt.nat3) goto 200
      enddo
      call close_file(iunit)
      goto 300
   endif
   goto 50

   300 return
end subroutine rdhess


subroutine write_tm_vibspectrum(ich, n3, freq, ir_int, raman_int)
   use xtb_setparam
   integer, intent(in) :: ich ! file handle
   integer, intent(in) :: n3
   real(wp), intent(in) :: freq(n3)
   real(wp), intent(in) :: ir_int(n3)
   real(wp), intent(in), optional :: raman_int(n3)
   integer :: i
   logical, allocatable :: iract(:), ramanact(:)
   real(wp) :: thr = 0.01_wp
   real(wp) :: thr_int = 1.0e-4_wp
   allocate(iract(n3),ramanact(n3))

   if (set%elprop == p_elprop_alpha) then
      write(ich,'("$vibrational spectrum")')
      write(ich,'("#  mode     symmetry     wave number   IR intensity   Raman activity    selection rules")')
      write(ich,'("#                         cm**(-1)      (km*mol⁻¹)      (Ä⁴*amu⁻¹)       IR     RAMAN")')
      do i = 1, n3
         if (raman_int(i) > thr_int ) then
            ramanact(i) = .true.
         else
            ramanact(i) = .false.
         endif
         if (ir_int(i) > thr_int ) then
            iract(i) = .true.
         else
            iract(i) = .false.
         endif
         if (abs(freq(i)).lt.thr) then
            write(ich,'(i6,9x,    f18.2,f16.5,f16.5,9x," - ",5x," - ")') &
               &  i, freq(i), 0.0_wp, 0.0_wp
         else
            if (iract(i) .and. ramanact(i)) then
               write(ich,'(i6,8x,"a",f18.2,f16.5,f16.5,9x,"YES",5x,"YES")') &
                  &  i,freq(i),ir_int(i),raman_int(i)
            elseif ((.not.iract(i)) .and. ramanact(i)) then
               write(ich,'(i6,8x,"a",f18.2,f16.5,f16.5,9x,"NO",5x,"YES")') &
                  &  i,freq(i),ir_int(i),raman_int(i)
            elseif (iract(i) .and. (.not. ramanact(i))) then
               write(ich,'(i6,8x,"a",f18.2,f16.5,f16.5,9x,"YES",5x,"NO")') &
                  &  i,freq(i),ir_int(i),raman_int(i)
            else
               write(ich,'(i6,8x,"a",f18.2,f16.5,f16.5,9x,"NO",5x,"NO")') &
                  &  i,freq(i),ir_int(i),raman_int(i)
            endif
         endif
      enddo
   else
      write(ich,'("$vibrational spectrum")')
      write(ich,'("#  mode     symmetry     wave number   IR intensity    selection rules")')
      write(ich,'("#                         cm**(-1)      (km*mol⁻¹)        IR     ")')
      do i = 1, n3
         if (ir_int(i) > thr_int ) then
             iract(i) = .true.
         else
             iract(i) = .false.
         endif
         if (abs(freq(i)).lt.thr) then
             write(ich,'(i6,9x,    f18.2,f16.5,9x," - ")') &
                 i,freq(i),0.0_wp
         else
             if (iract(i)) then
                 write(ich,'(i6,8x,"a",f18.2,f16.5,9x,"YES")') &
                     i,freq(i),ir_int(i)
             else
                 write(ich,'(i6,8x,"a",f18.2,f16.5,9x,"NO")') &
                     i,freq(i),ir_int(i)
             endif
         endif
      enddo
   endif

   write(ich,'("$end")')
end subroutine

subroutine g98fake2(fname,n,at,xyz,freq,red_mass,ir_int,u2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: freq(3*n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: u2(3*n,3*n)
   character(len=*),intent(in) :: fname
   real(wp),intent(in)  :: red_mass(3*n)
   real(wp),intent(in)  :: ir_int  (3*n)

   integer  :: gu,i,j,ka,kb,kc,la,lb,k
   character(len=2) :: irrep
   real(wp),allocatable :: u(:,:)
   real(wp),allocatable :: red(:)
   real(wp),allocatable :: f2 (:)
   real(wp),allocatable :: ir (:)
   real(wp) :: zero

   allocate( u(3*n,3*n), red(3*n), f2(3*n), ir(3*n), source = 0.0_wp )

   irrep='a'
   zero    =0.0

   k=0
   do i=1,3*n
      if(abs(freq(i)).gt.1.d-1)then
         k=k+1
         u(1:3*n,k)=u2(1:3*n,i)
         f2(k)=freq(i)
         ir(k)=ir_int(i)
         red(k)=red_mass(i)
      endif
   enddo

   gu=55
   call open_file(gu,fname,'w')
   write (gu,'('' Entering Gaussian System'')')
   write (gu,'('' *********************************************'')')
   write (gu,'('' Gaussian 98:'')')
   write (gu,'('' frequency output generated by the xtb code'')')
   write (gu,'('' *********************************************'')')

   write (gu,*) '                        Standard orientation:'
   write (gu,*) '---------------------------------------------', &
      & '-----------------------'
   write (gu,*) ' Center     Atomic     Atomic', &
      & '              Coordinates (Angstroms)'
   write (gu,*) ' Number     Number      Type ', &
      & '             X           Y           Z'
   write (gu,*) '-----------------------', &
      & '---------------------------------------------'
   j=0
   do i=1,n
      write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
   enddo
   write (gu,*) '----------------------', &
      & '----------------------------------------------'
   write (gu,*) '    1 basis functions        1 primitive gaussians'
   write (gu,*) '    1 alpha electrons        1 beta electrons'
   write (gu,*)
   111 format(i5,i11,i14,4x,3f12.6)

   write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities', &
      & ' (km*mol⁻¹),'
   write (gu,*) 'Raman scattering activities (A**4/amu),', &
      & ' Raman depolarization ratios,'
   write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
      & ' and normal coordinates:'

   ka=1
   kc=3
   60  kb=min0(kc,k)
   write (gu,100) (j,j=ka,kb)
   write (gu,105) (irrep,j=ka,kb)
   write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
   write (gu,110) ' Red. masses --',(red(j),j=ka,kb)
   write (gu,110) ' Frc consts  --',(zero,j=ka,kb)
   write (gu,110) ' IR Inten    --',(ir(j),j=ka,kb)
   write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
   write (gu,110) ' Depolar     --',(zero,j=ka,kb)
   write (gu,*)'Atom AN      X      Y      Z        X      Y', &
      & '      Z        X      Y      Z'
   la=1
   70  lb=n
   do  i=la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),  u(i*3-1,j),  u(i*3  ,j),j=ka,kb)
   enddo
   if (lb.eq.n) go to 90
   go to 70
   90  if (kb.eq.k) then
      return
   endif
   ka=kc+1
   kc=kc+3
   go to 60

   100 format (3(20x,i3))
   105 format (3x,3(18x,a5))
   110 format (a15,f11.4,12x,f11.4,12x,f11.4)
   130 format (2i4,3(2x,3f7.2))

   write(gu,'(''end of file'')')
   call close_file(gu)
   return

end subroutine g98fake2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine g98fake(fname,n,at,xyz,freq,u2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: freq(3*n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: u2(3*n,3*n)
   character(len=*),intent(in) :: fname

   integer  :: gu,i,j,ka,kb,kc,la,lb,k
   character(len=2) :: irrep
   real(wp),allocatable :: u(:,:)
   real(wp),allocatable :: red_mass(:)
   real(wp),allocatable :: force   (:)
   real(wp),allocatable :: ir_int  (:)
   real(wp),allocatable :: f2      (:)
   real(wp) :: zero

   allocate( u(3*n,3*n), red_mass(3*n), force(3*n), ir_int(3*n), f2(3*n), &
             source = 0.0_wp )

   irrep='a'
   red_mass=99.0
   force   =99.0
   ir_int  =99.0
   zero    =0.0

   k=0
   do i=1,3*n
      if(abs(freq(i)).gt.1.d-1)then
         k=k+1
         u(1:3*n,k)=u2(1:3*n,i)
         f2(k)=freq(i)
      endif
   enddo

   gu=55
   call open_file(gu,fname,'w')
   write (gu,'('' Entering Gaussian System'')')
   write (gu,'('' *********************************************'')')
   write (gu,'('' Gaussian 98:'')')
   write (gu,'('' frequency output generated by the xtb code'')')
   write (gu,'('' *********************************************'')')

   write (gu,*) '                        Standard orientation:'
   write (gu,*) '---------------------------------------------', &
      & '-----------------------'
   write (gu,*) ' Center     Atomic     Atomic', &
      & '              Coordinates (Angstroms)'
   write (gu,*) ' Number     Number      Type ', &
      & '             X           Y           Z'
   write (gu,*) '-----------------------', &
      & '---------------------------------------------'
   j=0
   do i=1,n
      write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
   enddo
   write (gu,*) '----------------------', &
      & '----------------------------------------------'
   write (gu,*) '    1 basis functions        1 primitive gaussians'
   write (gu,*) '    1 alpha electrons        1 beta electrons'
   write (gu,*)
   111 format(i5,i11,i14,4x,3f12.6)

   write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',' (km*mol⁻¹),'
   write (gu,*) 'Raman scattering activities (A**4/amu),', &
      & ' Raman depolarization ratios,'
   write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
      & ' and normal coordinates:'

   ka=1
   kc=3
   60  kb=min0(kc,k)
   write (gu,100) (j,j=ka,kb)
   write (gu,105) (irrep,j=ka,kb)
   write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
   write (gu,110) ' Red. masses --',(red_mass(j),j=ka,kb)
   write (gu,110) ' Frc consts  --',(force(j),j=ka,kb)
   write (gu,110) ' IR Inten    --',(ir_int(j),j=ka,kb)
   write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
   write (gu,110) ' Depolar     --',(zero,j=ka,kb)
   write (gu,*)'Atom AN      X      Y      Z        X      Y', &
      & '      Z        X      Y      Z'
   la=1
   70  lb=n
   do  i=la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),  u(i*3-1,j),  u(i*3  ,j),j=ka,kb)
   enddo
   if (lb.eq.n) go to 90
   go to 70
   90  if (kb.eq.k) then
      return
   endif
   ka=kc+1
   kc=kc+3
   go to 60

   100 format (3(20x,i3))
   105 format (3x,3(18x,a5))
   110 format (a15,f11.4,12x,f11.4,12x,f11.4)
   130 format (2i4,3(2x,3f7.2))

   write(gu,'(''end of file'')')
   call close_file(gu)
   return

end subroutine g98fake

end module xtb_freq_io
