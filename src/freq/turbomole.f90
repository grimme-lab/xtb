! This file is part of xtb.
!
! Copyright (C) 2021 Christoph Plett
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

!> Implementation of aoforce call for hessian

module xtb_freq_turbomole
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule

   implicit none

   character(len=*), parameter :: source = "aoforce_hessian"

   private

   public :: aoforce_hessian 

contains

subroutine aoforce_hessian(env, mol, h, dipd)
   use xtb_freq_io, only : writeHessianOut
   use xtb_io_reader, only : readHessian
   use xtb_type_reader, only : TReader
   use xtb_mctc_filetypes, only : fileType
   use xtb_setparam

   implicit none
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TReader) :: reader
   real(wp), intent(inout) :: h(:,:)
   real(wp), intent(inout) :: dipd(:,:)
   integer :: idipd, error

   call wrtm(mol%n,mol%at,mol%xyz) !Overwrite coord with RAM-xyz file
   call execute_command_line('exec aoforce > job.last2>> /dev/null')
   call reader%open('hessian')
   call readHessian(env, mol, h, reader, fileType%tmol)
   call reader%close

   call open_file(idipd,'dipgrad','r')
   if(idipd == -1) then
           call env%error("No dipolegradient found", source)
           return
   end if

   call read_dipgrad(idipd, mol%n, dipd, error)

   if(error /=0 ) then
           call env%error('An error occurred while reading the dipolegradient', source)
           return
   end if

   call close_file(idipd)

end subroutine aoforce_hessian


subroutine read_dipgrad(idipd, n, dipd, error)
   use xtb_mctc_systools

   implicit none

   integer, intent(out) :: error
   character(len=:), allocatable :: line
   integer, intent(in) :: n
   real(wp), intent(inout) :: dipd(:,:)
   integer, intent(in) :: idipd
   integer :: i,j

   error = 0

   do while(error == 0)
      call getline(idipd, line, error)
         if (index(line, '$dipgrad          cartesian dipole gradients') == 1) then
            do i=1, 3*n
               read(idipd,*) (dipd(j,i),j=1,3)
            end do
            exit 
         end if
   enddo
end subroutine read_dipgrad
 

end module xtb_freq_turbomole 
