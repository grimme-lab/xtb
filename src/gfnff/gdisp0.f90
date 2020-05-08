! somewhat stripped pure D3(BJ) routine which is faster by a factor of two compared to gdisp
    subroutine egdisp_gfnff(n,iz,xyz,sqrab,srab,npair,pairlist,&
     &                        disp,g,cn,dcnij)
      use gff_d3com
      use gff_param, only:d3r0,zetac6
      use tbpar_dftd3
      implicit none  

      integer,intent(in) :: n,iz(n),npair,pairlist(2,n*(n+1)/2)
      real*8,intent(in) :: xyz(3,n)
      real*8,intent(in) :: sqrab(n*(n+1)/2) 
      real*8,intent(in) ::  srab(n*(n+1)/2) 
      real*8,intent(inout) :: g(3,n)
      real*8,intent(in) :: dcnij(3,n,n)
      real*8,intent(in) :: cn(n)
 
      integer iat,jat,i,j,m,linij,ati,atj,lina
      real*8 R0,c6,disp,R06,R08
      real*8 r2,r,r4,r5,r6,r7,r8,t6,t8
      real*8 dc6_rest,r423
      real*8 dc6iji,dc6ijj
      real*8 rij(3)
      real*8 drij(n*(n+1)/2)  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
      real*8 dc6i(n)          ! dE_disp/dCN(iat) in dc6i(iat)
                              !dCN(iat)/d(r_ij) is equal to
                              !dCN(jat)/d(r_ij)    
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2        

      disp=0
      drij=0.0d0
      dc6i=0.0d0

!$omp parallel default(none) &
!$omp shared(npair,pairlist,srab,sqrab,r2r4,iz,cn,d3r0,drij,dc6i,disp,zetac6) &
!$omp private(m,iat,jat,linij,r,r2,r4,r5,r6,r7,r8,r423,R0,t6,t8,c6,dc6iji,dc6ijj,dc6_rest, &
!$omp&        ati,atj,r06,r08)
!$omp do REDUCTION (+:disp,dc6i,drij)
      do m=1,npair
          iat=pairlist(1,m)
          jat=pairlist(2,m)
          linij=jat+iat*(iat-1)/2        
          ati=iz(iat)
          atj=iz(jat)
          call getdC6gfnff(number_of_references(ati),number_of_references(atj),&
     &                     cn(iat),cn(jat),&
     &                     iz(iat),iz(jat),iat,jat,c6,dc6iji,dc6ijj)
          c6 = c6 * zetac6(linij)
          r2=sqrab(linij)
          r =srab (linij)
          r4=r2*r2
          r5=r4*r
          r6=r5*r 
          r7=r6*r
          r8=r7*r
          R0=d3r0(lina(ati,atj)) ! is R0^2!
!         R06=R0**6
!         R08=R06*R0*R0
          R06=R0*R0*R0
          R08=R06*R0
          t6=(r6+R06)
          t8=(r8+R08)
          r423=r2r4(ati)*r2r4(atj)*3.0d0*2.0d0 ! factor 2 = s8
          drij(linij)=drij(linij) &
     &             -C6*6.0d0*r5/(t6*t6) &
     &        -r423*C6*8.0d0*r7/(t8*t8)  
          dc6_rest=1.0d0/t6+r423/t8
          disp=disp+dc6_rest*c6  ! calculate E_disp
          dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
          dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
      enddo ! pair
!$omp enddo
!$omp end parallel

      disp=-disp

! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i       
! this part is much faster and OMP makes it slower in wall time

!!$omp parallel shared(npair,pairlist,srab,xyz,dcnij,drij,dc6i) private(m,iat,jat,linij,r,rij)
!!$omp do schedule(dynamic) REDUCTION (+:g)
      do m=1,npair
         iat=pairlist(1,m)
         jat=pairlist(2,m)
         linij=jat+iat*(iat-1)/2        
         r=srab(linij)
         rij=(xyz(:,jat)-xyz(:,iat))/r
         g(:,iat)=g(:,iat)+drij(linij)*rij
         g(:,jat)=g(:,jat)-drij(linij)*rij
      enddo ! pair
!!$omp enddo
!!$omp end parallel

      call dgemv('n',3*n,n,-1.0d0,dcnij,3*n,dc6i,1,1.0d0,g,1)

      end 


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      The   N E W   gradC6 routine    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      subroutine getdC6gfnff(mxci,mxcj,cni,cnj,&
     &                       izi,izj,iat,jat,c6,dc6i,dc6j) 
      use gff_d3com
      use tbpar_dftd3
      IMPLICIT NONE
      integer mxci,mxcj 
      integer iat,jat,izi,izj
      real*8 cni,cnj
      real*8 dc6i,dc6j,c6       

      real*8 k3,k32
      parameter (k3     =-4)
      parameter (k32    =-8)
      real*8 zaehler,nenner,dzaehler_i,dnenner_i,dzaehler_j,dnenner_j
      real*8 expterm,cn_refi,cn_refj,c6ref,r,term,dri,drj
      real*8 c6mem,r_save
      integer a,b

      c6mem=-1.d99
      r_save=9999.0
      zaehler=0.0d0
      nenner=0.0d0
      dzaehler_i=0.d0
      dnenner_i=0.d0
      dzaehler_j=0.d0
      dnenner_j=0.d0

      DO a=1,mxci
        DO b=1,mxcj
            c6ref = get_c6(a,b,izi,izj) !c6ab(1,a,b,izi,izj)
            if (c6ref.gt.0) then
              cn_refi=reference_cn(a,izi) !c6ab(2,a,b,izi,izj)
              cn_refj=reference_cn(b,izj) !c6ab(3,a,b,izi,izj)
              dri = cni-cn_refi
              drj = cnj-cn_refj
              r=dri**2+drj**2
              if (r.lt.r_save) then
                r_save=r
                c6mem=c6ref
              end if
              expterm=exp(k3*r)
              zaehler=zaehler+c6ref*expterm
              nenner=nenner+expterm
              expterm=expterm*k32
              term=expterm*dri           
              dzaehler_i=dzaehler_i+c6ref*term
              dnenner_i =dnenner_i +      term
              term=expterm*drj             
              dzaehler_j=dzaehler_j+c6ref*term
              dnenner_j =dnenner_j +      term
            end if
        ENDDO !b
      ENDDO !a

      if (nenner.gt.1.0d-9) then
        c6  =zaehler/nenner
        dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler)) /(nenner*nenner)
        dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler)) /(nenner*nenner)
      else
        c6=c6mem
        dc6i=0.0d0
        dc6j=0.0d0
      end if
!      nenner=nenner+1.d-9   ! can not be smaller than 1.d-10 in highly coordinated systems like metal clusters for cnmax=6
!      c6  =zaehler/nenner
!      dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler)) /(nenner*nenner)
!      dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler)) /(nenner*nenner)

      end 

module gffmod_dftd3
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: d3_gradient
   private

contains

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(nat, atoms, wf, cn, gwvec, gwdcn)
   use tbpar_dftd3
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk

   gwvec = 0.0_wp
   gwdcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, number_of_references(ati)
         gw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         norm = norm + gw
         dnorm = dnorm + 2*wf*(reference_cn(iref, ati) - cn(iat)) * gw
      end do
      norm = 1.0_wp / norm
      do iref = 1, number_of_references(ati)
         expw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         expd = 2*wf*(reference_cn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(reference_cn(:number_of_references(ati), ati)) &
               & == reference_cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references

!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references_d4(nat, atoms, wf, cn, gwvec, gwdcn)
   use dftd4
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk

   gwvec = 0.0_wp
   gwdcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, refn(ati)
         gw = weight_cn(wf, cn(iat), refcn(iref, ati))
         norm = norm + gw
         dnorm = dnorm + 2*wf*(refcn(iref, ati) - cn(iat)) * gw
      end do
      norm = 1.0_wp / norm
      do iref = 1, refn(ati)
         expw = weight_cn(wf, cn(iat), refcn(iref, ati))
         expd = 2*wf*(refcn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(refcn(:refn(ati), ati)) &
               & == refcn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references_d4

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(nat, atoms, gwvec, gwdcn, c6, dc6dcn)
   use tbpar_dftd3
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         do iref = 1, number_of_references(ati)
            do jref = 1, number_of_references(atj)
               refc6 = get_c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
      end do
   end do
end subroutine get_atomic_c6

!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6_d4(nat, atoms, gwvec, gwdcn, c6, dc6dcn)
   use tbpar_dftd3, only : get_c6
   use dftd4
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         do iref = 1, refn(ati)
            do jref = 1, refn(atj)
               refc6 = get_c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
      end do
   end do
end subroutine get_atomic_c6_d4

subroutine d3_gradient(nat, at, xyz, npair, pairlist, zeta_scale, radii, weighting_factor, &
      &                cn, dcndr, energy, gradient)
   use tbpar_dftd3
   use dftd4
   use mctcpar_r4r2, only: r4r2 => sqrt_z_r4_over_r2

   integer, intent(in) :: nat
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: xyz(:, :)
   integer, intent(in) :: npair
   integer, intent(in) :: pairlist(:, :)
   real(wp), intent(in) :: zeta_scale(:)
   real(wp), intent(in) :: radii(:)

   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp) :: sigma(3, 3)

   integer :: max_ref
   integer :: iat, jat, ati, atj, ij, img

   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   max_ref = maxval(refn(at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references_d4(nat, at, weighting_factor, cn, gw, dgwdcn)

   !gw = gw*zeta_scale
   !dgwdcn = dgwdcn*zeta_scale
   call get_atomic_c6_d4(nat, at, gw, dgwdcn, c6, dc6dcn)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(nat, at, xyz, npair, pairlist, zeta_scale, radii, c6, dc6dcn) &
   !$omp private(ij, img, iat, jat, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do img = 1, npair
      iat = pairlist(1, img)
      jat = pairlist(2, img)
      ij = jat + iat*(iat-1)/2
      ati = at(iat)
      atj = at(jat)
      rij = xyz(:, iat) - xyz(:, jat)
      r2 = sum(rij**2)

      r4r2ij = 3*r4r2(ati)*r4r2(atj)
      r0 = radii(lin(ati, atj))

      t6 = 1._wp/(r2**3+r0**3)
      t8 = 1._wp/(r2**4+r0**4)

      d6 = -6*r2**2*t6**2
      d8 = -8*r2**3*t8**2

      disp = (t6 + 2*r4r2ij*t8) * zeta_scale(ij)
      ddisp= (d6 + 2*r4r2ij*d8) * zeta_scale(ij)

      dE = -c6(iat, jat)*disp * 0.5_wp
      dG = -c6(iat, jat)*ddisp*rij
      dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

      energies(iat) = energies(iat) + dE
      dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
      sigma = sigma + dS
      if (iat /= jat) then
         energies(jat) = energies(jat) + dE
         dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         sigma = sigma + dS
      endif

   enddo
   !$omp end parallel do

   call dgemv('n', 3*nat, nat, 1.0_wp, dcndr, 3*nat, dEdcn, 1, 1.0_wp, gradient, 1)
   !call dgemv('n', 9, nat, 1.0_wp, dcndL, 9, dEdcn, 1, 1.0_wp, sigma, 1)

   energy = sum(energies)

end subroutine d3_gradient

pure elemental integer function lin(i1,i2)
   integer,intent(in) :: i1,i2
   integer :: idum1,idum2
   idum1=max(i1,i2)
   idum2=min(i1,i2)
   lin=idum2+idum1*(idum1-1)/2        
end function lin

real(wp) pure elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

real(wp) pure elemental function pair_scale(iat, jat) result(scale)
   integer, intent(in) :: iat, jat
   if (iat == jat) then
      scale = 0.5_wp
   else
      scale = 1.0_wp
   endif
end function pair_scale

end module gffmod_dftd3
