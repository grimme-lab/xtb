Subroutine get_schoenflies (n, iat, xyz, sfsym, paramar)
   Use iso_c_binding
   Implicit None
   integer, parameter :: wp = selected_real_kind(15, 307)

   Interface c_interface
      !Interface to c routine for symmetry recognition
      !attypes are Atom types as integers (e.g 6 for Carbon etc...)
      !x,y,z must be one dimensional sequential(!) arrays of real*8
      !symbol is the recognized schoenflies symbol
      Subroutine schoenflies (natoms, attypes, x, y, z, symbol, paramar) &
            &    bind (C, name="schoenflies")
         Use iso_c_binding
         import
         Implicit None
         Integer (c_int), Intent (In) :: natoms
         Type (c_ptr), value :: attypes
         Type (c_ptr), value :: x
         Type (c_ptr), value :: y
         Type (c_ptr), value :: z
         Type (c_ptr), Intent (InOut) :: symbol
         Type (c_ptr), value :: paramar
         !Character (kind=c_char), Intent (out)  :: symbol(*)
      End Subroutine schoenflies
   End Interface c_interface

   !Dummy Arguments
   Character (Len=*)     :: sfsym
   Integer, Intent (In)  :: n
   Integer, Intent (In)  :: iat (n)
   Real(wp), Intent (In) :: xyz (3, n)
   Real(wp), Intent (In) :: paramar (11)


   !local variables for passing to c routine:
   Integer (c_int) :: natoms
   Integer (c_int), Allocatable, Target, Dimension (:) :: attypes
   Real (c_double), Allocatable, Target, Dimension (:) :: x
   Real (c_double), Allocatable, Target, Dimension (:) :: y
   Real (c_double), Allocatable, Target, Dimension (:) :: z
   Real (c_double), Allocatable, Target, Dimension (:) :: c_paramar

   !passing strings to c and get back values requires some nasty pointer
   !stuff, as ISO-C treats strings as pointer-to-char
   Type (c_ptr) :: csptr !C-Compatible Pointer
   Character, Pointer, Dimension (:) :: fsptr !Fortran Pointer
   Integer :: sflen (1) !Needed for correct call of c_f_pointer

   !local stack:
   Integer :: i

   !Allocate memory for copies of iat, xyz...
   !As they will be passed to C as Pointers, they will be actually the working
   !memory for the C code. They must be allocated with "target" attribute.

   Allocate (attypes(n))
   Allocate (x(n))
   Allocate (y(n))
   Allocate (z(n))
   Allocate (c_paramar(11))

   !now, copy contents
   natoms = n
   attypes (:) = iat (:)
   x (:) = xyz (1, :)
   y (:) = xyz (2, :)
   z (:) = xyz (3, :)
   c_paramar = paramar

   sfsym = sfsym // C_NULL_CHAR
   sflen (1) = 3 !3 Characters soule be enough for schoenflies symbol
   csptr = c_loc (sfsym)! Get address of dummy string

   Call schoenflies (natoms, c_loc(attypes), c_loc(x), c_loc(y), &
      & c_loc(z), csptr, c_loc(c_paramar))

   Call c_f_pointer (csptr, fsptr, sflen)
   Write (sfsym,*) fsptr !Write symbol back to dummy argument
   Nullify (fsptr)

   !deallocate arrays:
   Deallocate (attypes, x, y, z,c_paramar)
End Subroutine get_schoenflies
