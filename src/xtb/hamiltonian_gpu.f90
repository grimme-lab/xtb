! This file is part of xtb.
!
! Copyright (C) 2020, NVIDIA CORPORATION. All rights reserved.
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

module xtb_xtb_hamiltonian_gpu
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert, only : evtoau
   use xtb_xtb_data, only : THamiltonianData
   use xtb_intgrad
   use xtb_lin
   use xtb_scc_core, only : shellPoly, h0scal
   use xtb_grad_core, only : dshellPoly
   implicit none
   private

   public :: build_SDQH0_gpu

   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer,parameter :: lx(84) = (/ &
      & 0, &
      & 1,0,0, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/)
   integer,parameter :: ly(84) = (/ &
      & 0, &
      & 0,1,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2/)
   integer,parameter :: lz(84) = (/ &
      & 0, &
      & 0,0,1, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2/)

   real(wp), parameter :: trafo(6,6) = reshape((/ & ! copied from scf.f, simplyfied
      ! --- dS
      & sqrt(1.0_wp/5.0_wp), &
      & sqrt(1.0_wp/5.0_wp), &
      & sqrt(1.0_wp/5.0_wp), &
      & 0.0_wp,0.0_wp,0.0_wp, &
      ! --- dx²-y²
      & 0.5_wp*sqrt(3.0_wp), &
      &-0.5_wp*sqrt(3.0_wp), &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
      ! --- dz²
      & 0.5_wp,0.5_wp,-1.0_wp, &
      & 0.0_wp,0.0_wp, 0.0_wp, &
      ! --- rest
      & 0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp, &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp, &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp /), shape(trafo))

contains

pure subroutine build_sdq_ints_gpu(a,b,c,alpi,alpj,la,lb,kab,t,e,lx,ly,lz,v)
   !$acc routine seq
   implicit none
   !     aufpunkte,ref point,intarray
   integer,intent(in)  :: la,lb
   real(wp), intent(in)  :: alpi,alpj
   real(wp), intent(in)  :: a(3),b(3),c(3)
   real(wp), intent(in)  :: kab,t(0:8),e(3)
   integer, intent(in)  :: lx(84),ly(84),lz(84)
   real(wp), intent(out) :: v(10)
   !     local variables
   real(wp)  :: dd(0:8),va(3),val(3,3)
   real(wp)  :: aa(0:3,3),bb(0:3,3)
   integer :: i,j,ij(3),ii(3),jj(3),lmax

   val = 0

   aa = 0
   bb = 0
   dd = 0
   ii = [lx(la),ly(la),lz(la)]
   jj = [lx(lb),ly(lb),lz(lb)]
   ij = ii+jj
   aa(ii(1),1)=1.0_wp
   aa(ii(2),2)=1.0_wp
   aa(ii(3),3)=1.0_wp
   bb(jj(1),1)=1.0_wp
   bb(jj(2),2)=1.0_wp
   bb(jj(3),3)=1.0_wp

   do i = 1, 3
      ! calculate cartesian prefactor for first gaussian
      call build_hshift2(aa(:,i),a(i),e(i),ii(i))    ! <a|
      ! calculate cartesian prefactor for second gaussian
      call build_hshift2(bb(:,i),b(i),e(i),jj(i))    ! |b>
      ! form their product
      call prod3(aa(:,i),bb(:,i),dd(:),ii(i),jj(i))
      lmax = ij(i)
      do j = 0, lmax
         !    <a|b> <a|x|b>             <a|x²|b>
         va = [t(j), t(j+1) + e(i)*t(j), t(j+2) + 2*e(i)*t(j+1) + e(i)*e(i)*t(j)]
         val(i,1:3) = val(i,1:3) + dd(j)*va(1:3)
      enddo
   enddo
   v( 1)=kab*(val(1,1)*val(2,1)*val(3,1))
   v( 2)=kab*(val(1,2)*val(2,1)*val(3,1))
   v( 3)=kab*(val(1,1)*val(2,2)*val(3,1))
   v( 4)=kab*(val(1,1)*val(2,1)*val(3,2))
   v( 5)=kab*(val(1,3)*val(2,1)*val(3,1))
   v( 6)=kab*(val(1,1)*val(2,3)*val(3,1))
   v( 7)=kab*(val(1,1)*val(2,1)*val(3,3))
   v( 8)=kab*(val(1,2)*val(2,2)*val(3,1))
   v( 9)=kab*(val(1,2)*val(2,1)*val(3,2))
   v(10)=kab*(val(1,1)*val(2,2)*val(3,2))
end subroutine build_sdq_ints_gpu

!> Computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine build_SDQH0_gpu(nShell, hData, nat, at, nbf, nao, xyz, trans, selfEnergy, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, &
      & sint, dpint, qpint, H0)
   implicit none
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)    :: trans(:, :)
   real(wp), intent(in) :: selfEnergy(:, :)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> Overlap integral matrix
   real(wp),intent(out) :: sint(nao,nao)
   !> Dipole integral matrix
   real(wp),intent(out) :: dpint(3,nao,nao)
   !> Quadrupole integral matrix
   real(wp),intent(out) :: qpint(6,nao,nao)
   !> Core Hamiltonian
   real(wp),intent(out) :: H0(:)


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,ioff,joff
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(6,6,3),qq(6,6,6),tmp(6,6),dd2(3,6,6)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   real(wp) :: gama,kab,t(0:8),e(3)
   real(wp),parameter ::  sqrtpi  = sqrt(pi)

   real(wp) s2(6,6),dum(6,6),sspher

   !$acc enter data create(H0(:),sint(nao,nao),dpint(3,nao,nao),qpint(6,nao,nao))

   !$acc kernels default(present)
   ! integrals
   H0(:) = 0.0_wp
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   !$acc end kernels
   ! --- Aufpunkt for moment operator
   point = 0.0_wp

   !$acc enter data copyin(xyz, at, nShell, selfEnergy, caoshell, saoshell, &
   !$acc& nprim, primcount, alp, cont, intcut, trans, point, hData, &
   !$acc& hData%principalQuantumNumber(:, :), hData%angShell(:, :), hData%valenceShell(:, :),  &
   !$acc& hData%numberOfPrimitives(:, :), hData%slaterExponent(:, :), hData%selfEnergy(:, :), &
   !$acc& hData%referenceOcc(:, :), hData%kCN(:, :), hData%electronegativity(:), hData%atomicRad(:), &
   !$acc& hData%shellPoly(:, :), hData%pairParam(:, :), hData%kQShell(:, :), hData%kQAtom(:))

   !$acc parallel vector_length(32) private(ss,dd,qq,rb,iat,jat,izp,ci,ra,& 
   !$acc& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$acc& jsh,jshmax,jshtyp,jcao,naoj,jptyp,shpoly, &
   !$acc& est,alpi,alpj,ab,iprim,jprim,ip,jp,il,jl,hii,hjj,km,zi,zj,zetaij,hav, &
   !$acc& mli,mlj,iao,jao,ii,jj,k,ij,itr,ioff,joff,t,e,saw,s2,dum)

   !$acc loop gang collapse(2)
   do iat = 1, nat
      do jat = 1, nat
         if (jat.ge.iat) cycle
         ra(1:3) = xyz(1:3,iat)
         izp = at(iat)
         jzp = at(jat)
         do ish = 1, nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, izp)
               zj = hData%slaterExponent(jsh, jzp)
               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish, izp).ne.0, &
                  & hData%valenceShell(jsh, jzp).ne.0,km)

               hav = 0.5_wp * km * (hii + hjj) * zetaij

               !$acc loop private(ss,dd,qq,s2,dum)
               do itr = 1, size(trans, dim=2)
                  !$acc cache(ss,dd,qq,s2,dum)
                  rb(1:3) = xyz(1:3,jat) + trans(:, itr)
                  rab2 = sum( (rb-ra)**2 )

                  ! distance dependent polynomial
                  shpoly=shellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                     &             hData%atomicRad(izp),hData%atomicRad(jzp),ra,rb)

                  ss = 0.0_wp
                  dd = 0.0_wp
                  qq = 0.0_wp

                  !call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ra,rb,point, &
                  !   &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
                  !$acc loop vector private(saw,t,e) independent collapse(2)
                  do ip = 1,nprim(icao+1)
                    do jp = 1,nprim(jcao+1)

                      iprim = ip+primcount(icao+1)
                      jprim = jp+primcount(jcao+1)
                      alpi = alp(iprim)
                      alpj = alp(jprim)

                      gama = alpi+alpj
                      ab = 1.0_wp/(gama)
                      e = (alpi * ra + alpj * rb)/(alpi + alpj)
                      t(0:8) = olapp([(j,j=0,8)],gama)
                      gama = alpi+alpj
                      ab = 1.0_wp/(gama)
                      est = rab2*alpi*alpj*ab
                      kab = exp(-est)*(sqrtpi*sqrt(ab))**3

                      do mli = 1,naoi
                        do mlj = 1,naoj
                          saw = 0.0_wp

                          ! prim-prim quadrupole and dipole integrals
                          call build_sdq_ints_gpu(ra,rb,point,alpi,alpj, &
                            & iptyp+mli,jptyp+mlj,kab,t,e,lx,ly,lz,saw)

                          iprim = ip+primcount(icao+mli)
                          jprim = jp+primcount(jcao+mlj)
                          ci = cont(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                          cc = cont(jprim)*ci
                          ! from primitive integrals fill CAO-CAO matrix for ish-jsh block
                          !                             ! overlap
                          !$acc atomic
                          ss(mlj,mli) = ss(mlj,mli)+saw(1)*cc
                          !$acc end atomic
                          ! dipole
                          do k=1,3
                            !$acc atomic
                            dd(mlj,mli,k) = dd(mlj,mli,k)+saw(k+1)*cc
                            !$acc end atomic
                          enddo
                          do k = 1,6
                            ! quadrupole
                            !$acc atomic
                            qq(mlj,mli,k) = qq(mlj,mli,k)+saw(k+4)*cc
                            !$acc end atomic
                          enddo
                        enddo ! mlj
                      enddo ! mli
                    enddo ! jp
                  enddo ! ip

                  !transform from CAO to SAO
                  !call dtrf2(ss,ishtyp,jshtyp)
                  ! DTRF2(S)
                  !select case(ishtyp)
                  !case(0) ! s-d
                  if (ishtyp.ge.2.or.jshtyp.ge.2) then
                    if (ishtyp.eq.0) then
                      !ss
                      do jj=1,6
                        sspher=0
                        do m=1,6
                          sspher=sspher+trafo(m,jj)*ss(m,1)
                        enddo
                        s2(jj,1)=sspher
                      enddo
                      ss(1:5,1) = s2(2:6,1)

                      !dd
                      do k = 1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*dd(m,1,k)
                          enddo
                          s2(jj,1)=sspher
                        enddo
                        dd(1:5,1,k) = s2(2:6,1)
                      enddo

                      !qq
                      do k = 1,6
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*qq(m,1,k)
                          enddo
                          s2(jj,1)=sspher
                        enddo
                        qq(1:5,1,k) = s2(2:6,1)
                      enddo

                    !case(1) ! p-d
                    else if (ishtyp.eq.1) then
                      !ss
                      do ii=1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*ss(m,ii)
                          enddo
                          s2(jj,ii)=sspher
                        enddo
                        ss(1:5,ii) = s2(2:6,ii)
                      enddo

                      ! dd
                      do k = 1,3
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*dd(m,ii,k)
                            enddo
                            s2(jj,ii)=sspher
                          enddo
                          dd(1:5,ii,k) = s2(2:6,ii)
                        enddo
                      enddo

                      !qq
                      do k = 1,6
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*qq(m,ii,k)
                            enddo
                            s2(jj,ii)=sspher
                          enddo
                          qq(1:5,ii,k) = s2(2:6,ii)
                        enddo
                      enddo
                      !return
                      !end select
                      !     wasn't there, then try iat ...
                      !select case(jshtyp)
                      !case(0) ! d-s
                    else if (jshtyp.eq.0) then
                      !ss
                      do jj=1,6
                        sspher=0
                        do m=1,6
                          sspher=sspher+trafo(m,jj)*ss(1,m)
                        enddo
                        s2(1,jj)=sspher
                      enddo
                      ss(1,1:5) = s2(1,2:6)

                      !dd
                      do k = 1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*dd(1,m,k)
                          enddo
                          s2(1,jj)=sspher
                        enddo
                        dd(1,1:5,k) = s2(1,2:6)
                      enddo

                      !qq
                      do k = 1,6
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*qq(1,m,k)
                          enddo
                          s2(1,jj)=sspher
                        enddo
                        qq(1,1:5,k) = s2(1,2:6)
                      enddo
                      !return
                      !case(1) ! d-p
                    else if (jshtyp.eq.1) then

                      !ss
                      do ii=1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*ss(ii,m)
                          enddo
                          s2(ii,jj)=sspher
                        enddo
                        ss(ii,1:5) = s2(ii,2:6)
                      enddo

                      !dd
                      do k = 1,3
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*dd(ii,m,k)
                            enddo
                            s2(ii,jj)=sspher
                          enddo
                          dd(ii,1:5,k) = s2(ii,2:6)
                        enddo
                      enddo

                      !qq
                      do k = 1,6
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*qq(ii,m,k)
                            enddo
                            s2(ii,jj)=sspher
                          enddo
                          qq(ii,1:5,k) = s2(ii,2:6)
                        enddo
                      enddo

                    else
                      !  return
                      !end select
                      !     if not returned up to here -> d-d
                      ! CB: transposing s in first dgemm is important for integrals other than S

                      !CALL mctc_gemm(trafo, s, dum, transa='T')

                      !ss -----------
                      do i=1,6
                        do j=1,6
                          sspher = 0.0_wp
                          do k=1,6
                            sspher = sspher + trafo(k,i) * ss(k,j)
                          enddo
                          dum(i,j) = sspher
                        enddo
                      enddo

                      !CALL mctc_gemm(dum, trafo, s2, transb='N')

                      do i=1,6
                        do j=1,6
                          sspher = 0.0_wp
                          do k=1,6
                            sspher = sspher + dum(i,k) * trafo(k,j)
                          enddo
                          s2(i,j) = sspher
                        enddo
                      enddo

                      do i=1,5
                        do j=1,5
                          ss(j,i) = s2(j+1,i+1)
                        enddo
                      enddo
                      ! end ss -----------

                      !dd --------
                      do kk = 1,3
                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + trafo(k,i) * dd(k,j,kk)
                            enddo
                            dum(i,j) = sspher
                          enddo
                        enddo

                        !CALL mctc_gemm(dum, trafo, s2, transb='N')

                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + dum(i,k) * trafo(k,j)
                            enddo
                            s2(i,j) = sspher
                          enddo
                        enddo

                        do i=1,5
                          do j=1,5
                            dd(j,i,kk) = s2(j+1,i+1)
                          enddo
                        enddo
                      enddo
                      !end dd --------------

                      !qq ---------------
                      do kk = 1,6
                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + trafo(k,i) * qq(k,j,kk)
                            enddo
                            dum(i,j) = sspher
                          enddo
                        enddo

                        !CALL mctc_gemm(dum, trafo, s2, transb='N')

                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + dum(i,k) * trafo(k,j)
                            enddo
                            s2(i,j) = sspher
                          enddo
                        enddo

                        do i=1,5
                          do j=1,5
                            qq(j,i,kk) = s2(j+1,i+1)
                          enddo
                        enddo
                      enddo
                      !end qq -----------------

                    endif
                  endif

                  ioff = saoshell(ish,iat)
                  joff = saoshell(jsh,jat)

                  ! would be good to change data order of sint,dpint,qpint for 
                  ! memory coalescence
                  !$acc loop vector collapse(2)
                  do ii = 1,llao2(ishtyp)
                     do jj = 1,llao2(jshtyp)
                       iao = ii + ioff
                       jao = jj + joff
                       ij = lin(iao, jao)
                       !$acc atomic
                       H0(ij) = H0(ij) + hav * shpoly * ss(jj, ii)
                       !$acc end atomic
                       !$acc atomic
                       sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                       !$acc end atomic
                       !$acc atomic
                       sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                       !$acc end atomic
                       do k = 1,3
                         !$acc atomic
                         dpint(k, iao, jao) = dpint(k, iao, jao) + dd(jj, ii, k)
                         !$acc end atomic
                         !$acc atomic
                         dpint(k, jao, iao) = dpint(k, jao, iao) + dd(jj, ii, k)
                         !$acc end atomic
                       enddo
                       do k = 1,6
                         !$acc atomic
                         qpint(k, iao, jao) = qpint(k, iao, jao) + qq(jj, ii, k)
                         !$acc end atomic
                         !$acc atomic
                         qpint(k, jao, iao) = qpint(k, jao, iao) + qq(jj, ii, k)
                         !$acc end atomic
                       enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$acc end parallel

   !$acc exit data copyout(H0(:),sint(nao,nao),dpint(3,nao,nao),qpint(6,nao,nao))

   ! dd array is transposed in ACC region, so dd2 is created for the diagonal
   ! elements, which is still running on CPU

   ! diagonal elements
   do iat = 1, nat
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            ii = i*(1+i)/2
            sint(i,i) = 1.0_wp + sint(i,i)
            H0(ii) = H0(ii) + selfEnergy(ish, iat)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = hData%angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd2 = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd2,qq)
            !transform from CAO to SAO
            !call dtrf2(ss,ishtyp,jshtyp)
            do k = 1,3
               tmp(1:6, 1:6) = dd2(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd2(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do k = 1,6
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do ii = 1, llao2(ishtyp)
               iao = ii + saoshell(ish,iat)
               do jj = 1, llao2(jshtyp)
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd2(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd2(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qq(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qq(1:6, jj, ii)
                  end if
               end do
            end do
         end do
      end do
   end do

   !$acc exit data delete(xyz, at, nShell, selfEnergy, caoshell, saoshell, &
   !$acc& nprim, primcount, alp, cont, intcut, trans, point)

   !$acc exit data delete(hData, &
   !$acc& hData%principalQuantumNumber(:, :), hData%angShell(:, :), hData%valenceShell(:, :),  &
   !$acc& hData%numberOfPrimitives(:, :), hData%slaterExponent(:, :), hData%selfEnergy(:, :), &
   !$acc& hData%referenceOcc(:, :), hData%kCN(:, :), hData%electronegativity(:), hData%atomicRad(:), &
   !$acc& hData%shellPoly(:, :), hData%pairParam(:, :), hData%kQShell(:, :), hData%kQAtom(:))

end subroutine build_SDQH0_gpu

end module xtb_xtb_hamiltonian_gpu
