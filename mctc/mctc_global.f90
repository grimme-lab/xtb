module mctc_global
   character(len=:),allocatable :: name  !< name of the currently running program
   character(len=:),allocatable,target :: msgbuffer !< error message buffer
   integer :: msgid !< number of generated errors
   integer :: maxmsg = 100 !< size of error buffer
   logical :: is_release = .false.
   ! set strict to turn all warnings into hard errors
   logical :: strict = .false.
   ! mute the error messages and generate warnings instead
   logical :: mute = .false.
   ! sanity status of the mctc environment, beware if it goes insane
   logical :: good = .true.

   type :: errormsg
      character(len=:),allocatable :: msg
      integer :: len
   endtype errormsg

   type(errormsg),allocatable :: errorbuffer(:)

contains

subroutine init_errorbuffer
   implicit none
   if (allocated(errorbuffer)) deallocate ( errorbuffer )
   allocate ( errorbuffer(maxmsg) )
   msgid = 0
   msgbuffer = ' '
end subroutine init_errorbuffer

end module mctc_global

