Subroutine get_schoenflies (n, iat, xyz, sfsym, paramar)
   Use iso_c_binding
   Implicit None
   integer, parameter :: wp = selected_real_kind(15, 307)

   Interface c_interface
      !Interface to c routine for symmetry recognition
      !attypes are Atom types as integers (e.g 6 for Carbon etc...)
      !coord must be ``one dimensional'' sequential(!) arrays of doubles
      !symbol is the recognized schoenflies symbol
      Subroutine schoenflies (natoms, attypes, coord, symbol, paramar) &
            &    bind (C, name="schoenflies")
         Use iso_c_binding
         import
         Implicit None
         Integer (c_int), Intent (In), value :: natoms
         integer(c_int), intent(in) :: attypes(*)
         real(c_double), intent(in) :: coord(3,*)
         Character (kind=c_char), Intent (out)  :: symbol(*)
         real(c_double), intent(in) :: paramar(*)
      End Subroutine schoenflies
   End Interface c_interface

   !Dummy Arguments
   Character (Len=*), target :: sfsym
   Integer, Intent (In)  :: n
   Integer, Intent (In)  :: iat (n)
   Real(wp), Intent (In) :: xyz (3, n)
   Real(wp), Intent (In) :: paramar (11)

   !local variables for passing to c routine:
   Integer (c_int) :: natoms
   Integer (c_int), Allocatable, Dimension (:) :: attypes
   Real (c_double), Allocatable, Dimension (:,:) :: coord
   Real (c_double), Allocatable, Dimension (:) :: c_paramar
   character(kind=c_char) :: symbol(6)

   !local stack:
   Integer :: i

   Allocate (attypes(n))
   Allocate (coord(3,n))
   Allocate (c_paramar(11))

   !now, copy contents
   natoms = n
   attypes = iat
   coord = xyz
   c_paramar = paramar
   symbol = C_NULL_CHAR

   Call schoenflies (natoms, attypes, coord, symbol, c_paramar)

   sfsym = ""
   do i = 1, size(symbol)
      if (symbol(i).eq.c_null_char) exit
      sfsym(i:i) = symbol(i)
   enddo

   !deallocate arrays:
   Deallocate (attypes, coord, c_paramar)
End Subroutine get_schoenflies
