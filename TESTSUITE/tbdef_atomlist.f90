subroutine test_atomlist
   use iso_fortran_env
   use assertion
   use tbdef_atomlist
   implicit none
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)
   integer, parameter :: atoms(*) = [3,1,1,5,8,1,1,2,5]
   logical, parameter :: lpar(*) = [.true., .false., .true., .true., .true., &
      &                             .false., .false., .true., .false.]
   integer, parameter :: ipar(*) = [1, 3, 4, 5, 8]
   character(len=*), parameter :: cpar = '1,3-5,8'

   write(error_unit,'(a)') " * Testing defaults"
   call assert(atl%get_truth() .eqv. .true.)
   call atl%switch_truth
   call assert(atl%get_truth() .eqv. .false.)
   call assert_eq(size(atl), 0)
   write(error_unit,'("-> Done:",1x,i0,1x,"fails")') afail

   write(error_unit,'(a)') " * Testing constructors"
   atl = tb_atomlist(list=lpar, truth=.true.)
   call assert_eq(size(atl), 9)
   call assert_eq(len(atl), 5)

   atl = tb_atomlist(list=lpar, truth=.false.)
   call assert_eq(size(atl), 9)
   call assert_eq(len(atl), 4)

   atl = tb_atomlist(list=ipar, truth=.true.)
   call assert_eq(size(atl), 8)
   call assert_eq(len(atl), 5)
   call atl%to_list(list)
   call assert_eq(list, ipar)

   atl = tb_atomlist(list=cpar, truth=.false.)
   call assert_eq(size(atl), 8)
   call assert_eq(len(atl), 5)
   call atl%to_string(string)
   call assert_eq(string, cpar)
   write(error_unit,'("-> Done:",1x,i0,1x,"fails")') afail
   call atl%new

   write(error_unit,'(a)') " * Testing data manipulation"
   call atl%new(lpar)
   call atl%resize(9)
   call atl%switch_truth
   call atl%to_string(string)
   call assert_eq(string, '2,6-7,9')
   call atl%gather(atoms, list)
   call assert_eq(list, [1,1,1,5])
   call assert_eq(size(list), len(atl))
   call atl%new

   atl = tb_atomlist(list=lpar, truth=.false., delimiter=' ', skip=':')
   call atl%to_string(string)
   call assert_eq(string, '2 6:7 9')
   call atl%switch_truth
   write(string, *) atl
   call assert_eq(string, '1 3:5 8')
   call atl%to_list(list)
   call atl%switch_truth
   call atl%add(list)
   call assert_eq(len(atl), size(atl))
   call atl%to_string(string)
   call assert_eq(string, '1:9')
   write(error_unit,'("-> Done:",1x,i0,1x,"fails")') afail

   call terminate(afail)
end subroutine test_atomlist
