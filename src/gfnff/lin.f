      integer function lin(i1,i2)
         integer,intent(in) :: i1,i2
         integer :: idum1,idum2
         idum1=max(i1,i2)
         idum2=min(i1,i2)
         lin=idum2+idum1*(idum1-1)/2        
         return
      end function lin
