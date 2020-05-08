!-----------------------------------------------------------------------------------
      subroutine trfp2xyz(nvar,nat3,p,xyz0,h,dspl)
      implicit none 
      integer nat3,nat,nvar
      real*8 xyz0(3,nat3/3),dspl(3,nat3/3),h(nat3,nat3),p(nvar)

      dspl = xyz0 
      call dgemv('n',nat3,nvar,1.0d0,h,nat3,p,1,1.d0,dspl,1)

      return
      end  

