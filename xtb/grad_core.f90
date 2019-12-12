! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module grad_core
   use iso_fortran_env, only : wp => real64

   use mctc_la

!! ========================================================================
!  we use the complete scc_core, so we inherit all interfaces and all
!  parameters. Therefore, we don't need to declare them again, note that
!  by including the grad_core you also inherit the scc_core!
   use scc_core

   implicit none

contains

!! ========================================================================
!  derivative of the CM5 additional term for GBSA in GFN1
!! ========================================================================
subroutine cm5_grad_gfn1(gradient, n, q, fgb, fhb, dcm5a, lhb)
   use mctc_econv, only: autoev, evtoau
   implicit none
   real(wp),intent(inout) :: gradient(:, :)
   integer, intent(in)    :: n
   real(wp),intent(in)    :: q(:)
   real(wp),intent(inout) :: fgb(:, :)
   real(wp),intent(inout) :: fhb(:)
   real(wp),intent(in)    :: dcm5a(:, :, :)
   logical, intent(in)    :: lhb

   real(wp), allocatable  :: fgba(:)
   real(wp), allocatable  :: dcm5(:, :)
   integer :: iat,jat

   allocate(fgba(n), dcm5(3, n), source=0.0_wp)

   fgb = fgb * evtoau
   fhb = fhb * evtoau

   do iat = 1, n
      do jat = 1, n
         fgba(iat)=fgba(iat)+q(jat)*fgb(jat,iat)
      enddo
   enddo
   do iat = 1, n
      do jat = 1, n
         gradient(:,iat)=gradient(:,iat)+fgba(jat)*dcm5a(:,jat,iat)
         dcm5(:,jat)=dcm5(:,jat)+dcm5a(:,jat,iat)
      enddo
   enddo
   do iat = 1, n
      gradient(:,iat)=gradient(:,iat)-fgba(iat)*dcm5(:,iat)
   enddo
   if(lhb) then
      do iat = 1, n
         fgba(iat)=fhb(iat)*2.0_wp*q(iat)
      enddo
      do iat = 1, n
         do jat = 1, n
            gradient(:,iat)=gradient(:,iat)+fgba(jat)*dcm5a(:,jat,iat)
         enddo
      enddo
      do iat = 1, n
         gradient(:,iat)=gradient(:,iat)-fgba(iat)*dcm5(:,iat)
      enddo
   endif

end subroutine cm5_grad_gfn1

!! ========================================================================
!  repulsion gradient of GFN1
!! ========================================================================
subroutine rep_grad_gfn1(mol, neighs, neighlist, kexp, rexp, &
      &                  energy, gradient, sigma)
   use aoparam, only : rep
   use tbdef_molecule
   use tbdef_neighbourlist
   type(tb_molecule), intent(in) :: mol
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out) :: energy
   real(wp), intent(in) :: kexp
   real(wp), intent(in) :: rexp

   integer  :: iat, jat, ati, atj, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, repab
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, rep, kexp, rexp)&
   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, alpha, repab, &
   !$omp&        t16, t26, t27, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)

         alpha = sqrt(rep(1, ati)*rep(1, atj))
         repab = rep(2, ati)*rep(2, atj)
         t16 = r1**kexp
         t26 = exp(-alpha*t16)
         t27 = r1**rexp
         dE = repab * t26/t27
         dG = -(alpha*t16*kexp + rexp) * dE * rij/r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iat) = energies(iat) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
         endif
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo
   !$omp end parallel do

   energy = sum(energies)

end subroutine rep_grad_gfn1

!! ========================================================================
!  repulsion gradient of GFN2
!! ========================================================================
subroutine rep_grad_gfn2(mol, neighs, neighlist, rexp, energy, gradient, sigma)
   use aoparam, only : rep
   use tbdef_molecule
   use tbdef_neighbourlist
   type(tb_molecule), intent(in) :: mol
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp),intent(inout) :: gradient(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(out) :: energy
   real(wp),intent(in) :: rexp

   integer  :: iat, jat, ati, atj, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, repab
   real(wp) :: r1, r2, rij(3), dE, dG(3), dS(3, 3)
   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, rep, rexp)&
   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, alpha, repab, &
   !$omp&        t16, t26, t27, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)

         alpha = sqrt(rep(1,ati)*rep(1,atj))
         repab = rep(2,ati)*rep(2,atj)
         t16 = r1**kexp(ati, atj)
         t26 = exp(-alpha*t16)
         t27 = r1**rexp
         dE = repab * t26/t27
         dG = -(alpha*t16*kexp(ati, atj) + rexp) * dE * rij / r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iat) = energies(iat) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
         endif
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo
   !$omp end parallel do

   energy = sum(energies)

contains

real(wp) pure elemental function kexp(ati, atj)
   integer, intent(in) :: ati, atj
   if(ati <= 2 .and. atj <= 2) then
      kexp = 1.0_wp
   else
      kexp = 1.5_wp
   endif
end function kexp

end subroutine rep_grad_gfn2

!! ========================================================================
!  shellwise electrostatic gradient for GFN1
!! ========================================================================
subroutine get_gfn_coulomb_derivs(mol, nshell, ash, gam, gtype, cf, lqpc, qsh, &
      &                           gradient, sigma)
   use tbdef_molecule
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Number of shells in the system.
   integer, intent(in) :: nshell
   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)
   !> GFN-method identifying the averaging function.
   integer, intent(in) :: gtype
   !> Chemical hardness of every shell.
   real(wp), intent(in) :: gam(:)
   !> Convergence for the Ewald summation (only used under PBC).
   real(wp), intent(in) :: cf
   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   select case(gtype)
   case(tb_gam_type%gfn1)
      if (mol%npbc > 0) then
         call coulomb_derivs_3d_impl(mol, nshell, ash, gam, gfn1_gam_average, &
            &                        cf, lqpc, qsh, gradient, sigma)
      else
         call coulomb_derivs_0d_impl(mol, nshell, ash, gam, gfn1_gam_average, &
            &                        qsh, gradient)
      endif
   case(tb_gam_type%gfn2)
      if (mol%npbc > 0) then
         call coulomb_derivs_3d_impl(mol, nshell, ash, gam, gfn2_gam_average, &
            &                        cf, lqpc, qsh, gradient, sigma)
      else
         call coulomb_derivs_0d_impl(mol, nshell, ash, gam, gfn2_gam_average, &
            &                        qsh, gradient)
      endif
   end select

end subroutine get_gfn_coulomb_derivs

subroutine coulomb_derivs_0d_impl(mol, nshell, ash, gam, gav, qsh, gradient)
   use tbdef_molecule
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Number of shells in the system.
   integer, intent(in) :: nshell
   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)
   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav
   !> Chemical hardness of every shell.
   real(wp),intent(in) :: gam(:)
   real(wp),intent(in) :: qsh(:)
   real(wp),intent(inout) :: gradient(:,:)

   integer  :: is, js, iat, jat
   real(wp) :: gi, gj, r2, rij(3), xij, dG(3)

   do is = 1, nshell
      iat = ash(is)
      gi = gam(is)
      do js = 1, is-1
         jat = ash(js)
         if(jat >= iat) cycle
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(rij**2)
         gj = gam(js)
         xij = gav(gi, gj)
         dG = -qsh(is)*qsh(js)/sqrt(r2 + xij**2)**3 * rij
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo

end subroutine coulomb_derivs_0d_impl

subroutine coulomb_derivs_3d_impl(mol, nshell, ash, gam, gav, cf, lqpc, qsh, &
      &                           gradient, sigma)
   use tbdef_molecule
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Number of shells in the system.
   integer, intent(in) :: nshell
   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)
   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav
   !> Chemical hardness of every shell.
   real(wp),intent(in) :: gam(:)
   !> Convergence for the Ewald summation (only used under PBC).
   real(wp),intent(in) :: cf
   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc
   real(wp),intent(in) :: qsh(:)
   real(wp),intent(inout) :: gradient(:, :)
   real(wp),intent(inout) :: sigma(:, :)

   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp

   integer :: is, js, iat, jat, img
   real(wp) :: gi, gj, r2, ri(3), rj(3), rw(3), riw(3), qpc, xij, wqq
   real(wp) :: dG(3), dS(3, 3), tG(3), tS(3, 3)

   qpc = 0.0_wp

   do is = 1, nshell
      iat = ash(is)
      ri = mol%xyz(:, iat)
      gi = gam(is)
      do js = 1, is
         jat = ash(js)
         rj = mol%xyz(:, jat)
         gj = gam(js)
         xij = gav(gi, gj)
         if (lqpc) qpc = xij
         dG = 0.0_wp
         dS = 0.0_wp
         do img = 1, mol%wsc%itbl(jat, iat)
            rw = rj + matmul(mol%lattice, mol%wsc%lattr(:, img, jat, iat))
            riw = ri - rw
            wqq = mol%wsc%w(jat,iat)*qsh(is)*qsh(js)
            call gfn_ewald_dx_3d_rec(riw, ewaldCutR, mol%rec_lat, qpc, &
               &                     mol%volume, cf, tG, tS)
            dG = dG + tG * wqq
            dS = dS + tS * wqq
            call gfn_ewald_dx_3d_dir(riw, ewaldCutD, mol%lattice, xij, qpc, cf, &
               &                     tG, tS)
            dG = dG + tG * wqq
            dS = dS + tS * wqq
         enddo
         if (is == js) then
            sigma = sigma + dS*0.5_wp
         else
            if (iat /= jat) then
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
            endif
            sigma = sigma + dS
         endif
      enddo
   enddo

end subroutine coulomb_derivs_3d_impl

pure subroutine gfn_ewald_dx_3d_rec(riw,rep,rlat,qpc,vol,cf,dG,dS)
   use mctc_constants
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(out) :: dG(3) !< element of interaction matrix
   real(wp),intent(out) :: dS(3,3)
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3),dtmp,fpivol
   real(wp) :: expterm,arg
   real(wp), parameter :: unity(3, 3) = reshape(&
      &[1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & shape(unity))
   dG = 0.0_wp
   dS = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3), &
         &       dx/=0 .or. dy/=0 .or. dz/=0)
      rik = matmul(rlat, [dx,dy,dz])
      rik2 = sum(rik**2)
      expterm = exp(-rik2/(4.0_wp*cf**2))/rik2
      arg = dot_product(rik,riw)
      ! d/dx (sin**2 + cos**2) = -2*sin*cos - 2*cos*sin
      dtmp = -sin(arg)*expterm
      dG = dG + rik*dtmp
      dS = dS + cos(arg)*expterm * (spread(rik, 1, 3)*spread(rik, 2, 3) &
         & * (-2.0_wp/rik2 - 0.5_wp/cf**2 + 4*qpc**2) &
         & - unity*(1.0_wp + 2*rik2*qpc**2))
   end do
   fpivol = 4.0_wp*pi/vol
   dG = dG * fpivol
   dS = dS * fpivol
end subroutine gfn_ewald_dx_3d_rec

pure subroutine gfn_ewald_dx_3d_dir(riw,rep,dlat,xij,qpc,cf,dG,dS)
   use mctc_constants
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: xij       !< interaction radius
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dG(3) !< element of interaction matrix
   real(wp),intent(out) :: dS(3,3)
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz
   real(wp) :: r2,r1,rij(3),arg2
   real(wp) :: t(3),dtmp,stmp(3)
   dG = 0.0_wp
   dS = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3))
      ! real contributions
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
      r1 = norm2(rij)
      r2 = r1*r1
      if(r1 < eps) cycle ! self-interaction handled elsewhere
      arg2 = cf**2 * r2
      dtmp = - 2*cf*exp(-arg2)/(sqrtpi*r2) + erf(cf*r1)/(r2*r1) &
         &   + 0.5_wp*qpc**2 * ((4*cf**3*r2 + 6*cf)*exp(-arg2)/(sqrtpi*r2*r2) &
         &                      - 3*erf(cf*r1)/(r2*r2*r1)) &
         &   - 1.0_wp/sqrt(r1**2 + xij**2)**3
      dG = dG + dtmp*rij
      dS = dS + dtmp*spread(rij, 1, 3)*spread(rij, 2, 3)
   enddo

end subroutine gfn_ewald_dx_3d_dir

!! ========================================================================
!  dH0/dxyz part from S(R) enhancement
!! ========================================================================
pure subroutine dhdr(n,at,xyz,iat,i,j,ati,atj,ishell,jshell,rab2,rf,dHdxyz)
   implicit none
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   real(wp), intent(in)  :: xyz(3,n)
   integer,intent(in)  :: iat
   integer,intent(in)  :: i
   integer,intent(in)  :: j
   integer,intent(in)  :: ishell
   integer,intent(in)  :: jshell
   integer,intent(in)  :: ati
   integer,intent(in)  :: atj
   real(wp), intent(in)  :: rab2
   real(wp), intent(out) :: rf
   real(wp), intent(out) :: dHdxyz(3)

   call drfactor(ishell,jshell,iat,at(ati),at(atj),rab2, &
      &          xyz(:,ati),xyz(:,atj),rf,dHdxyz)
   !if(atj.eq.iat) dHdxyz=-dHdxyz

end subroutine dhdr

!! ========================================================================
!  derivative of S(R) enhancement factor
!! ========================================================================
pure subroutine drfactor(ish,jsh,iat,ati,atj,rab2,xyz1,xyz2,rf,dxyz)
   use mctc_econv
   use aoparam, only : rad,polyr
   implicit none
   integer,intent(in)  :: ati,atj,ish,jsh,iat
   real(wp), intent(in)  :: rab2
   real(wp), intent(out) :: dxyz(3),rf
   real(wp), intent(in)  :: xyz1(3),xyz2(3)
   real(wp) :: rab,k1,k2,rr,r,dum,rf1,rf2,dr(3)
   real(wp) :: t14,t15,t17,t22,t20,t23,t10,t11,t35,t13
   ! R^a dependence 0.5 in GFN1
   real(wp), parameter :: a = 0.5_wp

   dr=xyz1-xyz2

   rab=sqrt(rab2)

   ! this sloppy conv. factor has been used in development, keep it
   r=(rad(ati)+rad(atj))*aatoau

   rr=rab/r

   k1=polyr(ish,ati)*0.01_wp
   k2=polyr(jsh,atj)*0.01_wp

   t22 = rr**a
   t15 = k1*t22
   t17 = 1/rab2
   t23 = k2*t22
   rf=(1.0_wp+t15)*(1.0_wp+k2*t22)
   dxyz=(t15*a*t17*(1.+t23)+(1.+t15)*k2*t22*a*t17)*dr

end subroutine drfactor

end module grad_core
