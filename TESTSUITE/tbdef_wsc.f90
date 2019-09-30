subroutine test_wigner_seitz

   ipar = 0
   do i = 1, mol%n
      do j = 1, i
         ipar = ipar + mol%wsc%itbl(j,i)
      enddo
   enddo
   call assert_eq(ipar,mol%wsl%np)

end subroutine test_wigner_seitz
