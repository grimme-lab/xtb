! This file is part of xtb.
!
! Copyright (C) 2022 Stefan Grimme, Christoph Plett
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

module xtb_iff_iffenergy
   use xtb_mctc_accuracy, only: wp
   use xtb_type_environment, only: TEnvironment
   use xtb_docking_param
   use xtb_sphereparam, only: number_walls
   use xtb_sphereparam, only : sphere_alpha, wpot, polynomial_cavity_list, sphere, cavity_egrad,&
           & cavitye, maxwalls, rabc, sphere
   implicit none

   private
   public :: iff_e, intermole_probe, alignmol

contains

   subroutine iff_e(env, n, n1, n2, at1, at2, neigh, A1, c02, q01, q02, c6ab,&
                            & z1, z2, nl1, nl2, l1, l2, AL1, c0l2,&
                            & qdr1, qdr2,&
                            & cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                            & den1, den2, gab1, gab2,&
                            & pr, mode, e, icoord)
      type(TEnvironment), intent(inout) :: env
      integer, intent(in) :: n, n1, n2
      integer, intent(in) :: nl1, nl2
      integer, intent(in) :: at1(n1)
      integer, intent(in) :: at2(n2)
      integer, intent(in) :: mode
      integer, intent(in) :: neigh(0:n, n)
      !> mol1 xyz
      real(wp), intent(in) :: A1(3, n1)    
      !> mol2 xyz
      real(wp), intent(in) :: c02(3, n2)   
      !> mol1 LMO
      real(wp), intent(in) :: AL1(4, n1*10)
      !> mol2 LMO
      real(wp), intent(in) :: c0l2(4, n2*10)
      real(wp), intent(in) :: q01(n1)
      real(wp), intent(in) :: qdr1(n1)
      real(wp), intent(in) :: q02(n2)
      real(wp), intent(in) :: qdr2(n2)
      real(wp), intent(in) :: cn1(n1)
      real(wp), intent(in) :: z1(n1)
      real(wp), intent(in) :: z2(n2)
      real(wp), intent(in) :: cn2(n2)
      real(wp), intent(in) :: alp1(n1)
      real(wp), intent(in) :: alp2(n2)
      real(wp), intent(in) :: c6ab(n, n)
      real(wp), intent(in) :: alpab(n2, n1)
      real(wp), intent(in) :: qct1(n1, 2), qct2(n2, 2)
      real(wp), intent(in) :: den1(2,4,n1),den2(2,4,n2),gab1(n1,n1),gab2(n2,n2)
      integer, intent(in) :: l1(n1*10)
      integer, intent(in) :: l2(n2*10)
      logical, intent(in) :: pr
      real(wp), optional, intent(in) :: icoord(6)
      real(wp), intent(out) :: e

      real(wp) :: ed,es,ep,esl,ei,eabc,ect,ep2,esph,eadd,gborn,gsolv,dgrrho
      real(wp) :: e_const
      real(wp) :: r0,rr,av,c6,r2,r6,r8,r42,er2,edum,symnum,qq,dhrrho,dzpve
      real(wp) :: R06, rmin, f, dc6, t6, t8, qt1, qt2, tmpr, rrr, p72, nn12
      real(wp) :: d2(3), t1, t2, t3, c9, ang, damp, tmp1, tmp2, sum1, sum2
      real(wp) :: c2(3, n2), cl2(4, 10*n2), dip(3), h, zz1, zz2, nn1, nn2, ees
      real(wp) :: rab(n2, n1), rotm(3, 3), q1(n1), q2(n2), AL2(4, n2*10), A2(3, n2)
      real(wp) :: r2ab(n2, n1), r0tmp(n2, n1), r0tmp2(n2, n1), r, oner, sab(n2, n1)

      !> QCG stuff
      real(wp),allocatable :: rab_solu_solv(:,:)

      !> Thermo stuff
      real(wp) :: aa, bb, cc, avmom, wt, h298, g298, ts298, zp, symn, freq(3*n)
      real(wp) :: xyz(3, n)
      real(wp) :: cavgrad(3) = 0.0_wp
      real(wp), parameter ::autoang = 0.52917726d0
      real(wp), parameter ::au = 627.509541d0
      integer :: pair(2, n1*n2)
      integer :: npair
      integer :: i, j, k, kk, i1, j2, iat, jat, n3, metal(86), counter
      integer :: at(n)
      logical :: linear, ATM, ijh, ijx
      character(len=4) :: pgroup

      !> Initial stuff
      metal = 1
      metal(1:2) = 0
      metal(6:10) = 0
      metal(14:18) = 0
      metal(32:36) = 0
      metal(50:54) = 0
      metal(82:86) = 0
      e_const=0.0_wp

      ! new cartesians if icoord is present
      if (present(icoord)) then
         call rotmat(icoord, rotm) !write(env%unit icoord(4:6) with new angles
         call move(n2, nl2, c02, c0l2, A2, AL2, rotm, icoord)      !A2 new coords
      else
         A2 = c02
         AL2 = c0l2
      end if

      ATM = .false.
      if (pr .or. fulle .or. mode .eq. 0) ATM = .true.

      ed = 0
      ep = 0
      es = 0
      esl = 0
      ei = 0
      eabc = 0
      eadd = 0
      esph = 0
      gsolv = 0
      dgrrho = 0
      npair = 0
      pair = 0
      counter = 0

      !> Distances of solute and screened solvent positions for QCG attractive potential
      if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
         allocate(rab_solu_solv(directedset%n, n2), source=0.0_wp)
      end if

      !> Distance between all atoms of molA and molB
      do i = 1, n1
         !> Check if qcg mode and i is a solute atom
         if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
            if(any(i == directedset%atoms)) then
               counter = counter + 1
            end if
         end if
         do j = 1, n2
            r2 = (A1(1, i) - A2(1, j))**2&
             &+ (A1(2, i) - A2(2, j))**2&
             &+ (A1(3, i) - A2(3, j))**2
            r2ab(j, i) = r2 + 1.d-12
            rab(j, i) = sqrt(r2) + 1.d-12
            if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
               if(any(i == directedset%atoms)) then
                  !Safe distance in Angström only, if solute atom
                  ! (important if molA is solute with already added solvents)
                  rab_solu_solv (counter,j) = sqrt(r2)*autoang 
               end if
            end if
            if (r2 .gt. 12000.0d0 .or. r2 .lt. 5.0d0) cycle  ! induction cut-off
            if (metal(at1(i)) .eq. 1 .or. metal(at2(j)) .eq. 1) cycle
            npair = npair + 1
            pair(1, npair) = i
            pair(2, npair) = j
         end do
      end do

!     atomic density overlap
      call ovlp(n1, n2, r2ab, den1, den2, sab)

!     CT
      q1 = q01
      q2 = q02
      call chrgtransfer(env, n1, n2, at1, at2, q1, q2, qct1, qct2, gab1, gab2,&
                       &alpab, rab, debug .and. pr, ect)

!     interatomic ES part
      p72 = par_es_damp
      do i = 1, n1
         do j = 1, n2
            r = rab(j, i)
            r0tmp(j, i) = p72*alpab(j, i)
            if (r .gt. 8) then                 ! point charge cut-off, good to 0.02 kcal
               qq = q1(i)*q2(j)
               es = es + qq/r
            else
               zz1 = z1(i)
               zz2 = z2(j)
               nn1 = zz1 - q1(i)
               nn2 = zz2 - q2(j)
               call truees(n1, n2, i, j, r, r2ab(j, i), den1, den2,&
                          &zz1, zz2, nn1, nn2, ees)  ! via density
               es = es + ees
            end if
         end do
      end do

!     D3 + Pauli via density
      ep2 = 0
      do i = 1, n1
         iat = at1(i)
         nn1 = z1(i) - q1(i)
         do j = 1, n2
            r6 = r2ab(j, i)**3
            r8 = r2ab(j, i)*r6
            jat = at2(j)
            r42 = rrab(jat, iat)
            dc6 = 1./(r6 + r0ab6(jat, iat)) + r42/(r8 + r0ab8(jat, iat))
            ed = ed + dc6*c6ab(j + n1, i)                          ! disp
            if(directedset%n > 0 .and. directed_type == p_atom_att) then
               ed = ed - directedset%val(i) * (exp(-0.01*((rab(j,i)-2)**2)))
               !- Because ed is substracted at the end
            end if

            if (rab(j, i) .gt. 25) cycle
            nn12 = nn1*(z2(j) - q2(j))
            ep = ep + nn12*sab(j, i)/rab(j, i)      ! ovlp S(S+1) is slightly worse, S(S-1) very slightly better
            if (iat .eq. 1 .or. jat .eq. 1)&
           &ep2 = ep2 + sqrt(nn12)*exp(-par_xh1*rab(j, i)*alpab(j, i)) ! X-H corr.
            if(directedset%n > 0 .and. directed_type == p_atom_pot) then
               ep = ep + directedset%val(i) * (1/rab(j,i)**pot_expo) * pot_scal !Add potential for directed docking
            end if
         end do
      end do

!     D3 ATM term (simplified i.e. no damping on i..k and j..k)
      if (ATM) then
         do i = 1, n1
            iat = at1(i)
            do j = 1, n2
               r2 = r2ab(j, i)
               if (r2 .gt. 500.0d0) cycle  ! ATM cut-off E^(3) accurate to about 1 %
               jat = at2(j)
               d2(1) = r2
               r = rab(j, i)
               rrr = r/alpab(j, i)
               damp = 1.0d0/(1.d0 + 6.d0*rrr**(-8))
               do kk = 1, neigh(0, i)
                  k = neigh(kk, i)   ! k on i
                  c9 = c6ab(j + n1, i)*c6ab(i, k)*c6ab(j + n1, k)
                  d2(2) = r2ab(j, k)
                  d2(3) = ((A1(1, i) - A1(1, k))**2&
                           &+ (A1(2, i) - A1(2, k))**2&
                           &+ (A1(3, i) - A1(3, k))**2)
                  t1 = (d2(1) + d2(2) - d2(3))
                  t2 = (d2(1) + d2(3) - d2(2))
                  t3 = (d2(3) + d2(2) - d2(1))
                  tmp2 = d2(1)*d2(2)*d2(3)
                  ang = (0.375d0*t1*t2*t3/tmp2 + 1.0d0)/tmp2**1.5d0
                  eabc = eabc + damp*ang*sqrt(c9)
               end do
               do kk = 1, neigh(0, j + n1)
                  k = neigh(kk, j + n1)   ! k on j
                  c9 = c6ab(j + n1, i)*c6ab(j + n1, k + n1)*c6ab(i, k + n1)
                  d2(2) = ((A2(1, j) - A2(1, k))**2&
                           &+ (A2(2, j) - A2(2, k))**2&
                           &+ (A2(3, j) - A2(3, k))**2)
                  d2(3) = r2ab(k, i)
                  t1 = (d2(1) + d2(2) - d2(3))
                  t2 = (d2(1) + d2(3) - d2(2))
                  t3 = (d2(3) + d2(2) - d2(1))
                  tmp2 = d2(1)*d2(2)*d2(3)
                  ang = (0.375d0*t1*t2*t3/tmp2 + 1.0d0)/tmp2**1.5d0
                  eabc = eabc + damp*ang*sqrt(c9)
               end do
            end do
         end do
      end if

! LP - LP  corrections
      do i = 1, nl1
         iat = l1(i)
         do j = 1, nl2
            jat = l2(j)
            r2 = (AL1(1, i) - AL2(1, j))**2&
              &+ (AL1(2, i) - AL2(2, j))**2&
              &+ (AL1(3, i) - AL2(3, j))**2
            r = sqrt(r2)
            esl = esl + AL1(4, i)*AL2(4, j)*r/(r2 + r0tmp(jat, iat))
         end do
!        LP atom corrections
         do j = 1, n2
            r2 = (A2(1, j) - AL1(1, i))**2&
              &+ (A2(2, j) - AL1(2, i))**2&
              &+ (A2(3, j) - AL1(3, i))**2
            r = sqrt(r2)
            esl = esl + AL1(4, i)*q2(j)*r/(r2 + r0tmp(j, iat))
         end do
      end do
!     LP atom corrections
      do i = 1, nl2
         iat = l2(i)
         do j = 1, n1
            r2 = (A1(1, j) - AL2(1, i))**2&
              &+ (A1(2, j) - AL2(2, i))**2&
              &+ (A1(3, j) - AL2(3, i))**2
            r = sqrt(r2)
            esl = esl + AL2(4, i)*q1(j)*r/(r2 + r0tmp(iat, j))
         end do
      end do

!     Drude Oscillator energy
      call drudescf(n1, n2, nl1, nl2, at1, at2, npair, pair, &
                    A1, q1, qdr1, A2, q2, qdr2, alpab, r0tmp2, ei)

!     Cavity Energies (For Wall Potentials)
      esph = 0.0_wp
      if(qcg .and. (number_walls > 0)) then
         do i=1, maxwalls
            if(.not. allocated(wpot(i)%list)) then
               rabc(1:3) = wpot(i)%radius(1:3)
               exit
            end if
         end do
         sphere = 2
         call cavitye(n1,A1,esph)
         call cavitye(n2,A2,esph)
      else
         at(1:n1) = at1(1:n1)
         at(n1+1:n) = at2(1:n2)
         xyz(1:3, 1:n1) = A1(1:3, 1:n1)
         xyz(1:3, n1+1:n) = A2(1:3, 1:n2)
         call cavity_egrad(n, at, xyz, esph, cavgrad)
     end if

      !> Final Energies
      ep = ep*par_rep_scal + ep2*par_xh2*0.01

      e = -ed + es + ep + esl + ei + gsolv + eabc + ect + esph

      !For QCG v2, the closer the solvent to the solute, the more attractive the energy is
      if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
         !Energie in Hartree
         !Potential of type (a/2)*erf(-(b*r-c*b))+(a/2)
         !a=E_int solv-solv and is the maximum value of potential
         !b=slope of error function
         !c=mid point of error function
         e = e - &
             & ((directedset%val(1)/2) * &
             & erf(-(directedset%expo(1) * minval(rab_solu_solv)) + (directedset%expo(2) * directedset%expo(1))) &
             & + (directedset%val(1)/2))
!            &directedset%val(1) * (exp(-directedset%expo(1) * (minval(rab_solu_solv))**2))
      end if

      if(isnan(e)) e = 1.0_wp**42

! all done now output
      if (pr .and. e .lt. 10000.2) then ! only if it makes sense to print
         call rcma(n1, A1, at1, n2, A2, at2, r, rmin)
         write (env%unit, '(''R CMA (Angst) :'',F10.3)') r*autoang
         write (env%unit, '(''R min (Angst) :'',F10.3)') rmin*autoang
         call diptot(n1, n2, nl1, nl2, A1, A2, AL1, AL2, l1, l2, q1, q2, dip)
         write(env%unit,'(''dipole moment :'',F10.3)') sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
         write (env%unit, '(''intermolecular energies in kcal/mol'')')
         write (env%unit, '(''E Pauli       :'',F10.3)') ep*au
         write (env%unit, '(''E disp ATM    :'',F10.3)') eabc*au
         write (env%unit, '(''E disp 2B     :'',F10.3)') - ed*au
         write (env%unit, '(''E disp total  :'',F10.3)') - ed*au + eabc*au
         write (env%unit, '(''E ES atom     :'',F10.3)') es*au
         write (env%unit, '(''E ES LMO      :'',F10.3)') esl*au
         write (env%unit, '(''E ES total    :'',F10.3)') (es + esl)*au
         write (env%unit, '(''E induction   :'',F10.3)') ei*au
         write (env%unit, '(''E CT          :'',F10.3)') ect*au
         if (sphere .ne. 0) write (env%unit, '(''E cavity      :'',F10.3)') esph*au
         write (env%unit, '(''Eint total,gas:'',F10.3)') (e - gsolv)*au
        write (env%unit, '(''              '',F14.8,''  <== Gint total'')') e*au
      end if

   end subroutine iff_e

   subroutine drudescf(n1,n2,nl1,nl2,at1,at2,npair,pair,xyz1,qat1,qdr1,xyz2,&
                      &qat2, qdr2, alpab, r0, edr)
      integer, intent(in) :: n1, n2, npair, pair(2, n1*n2), nl1, nl2
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)
      real(wp), intent(in) :: qdr1(n1)
      real(wp), intent(in) :: qat1(n1)
      real(wp), intent(in) :: qdr2(n2)
      real(wp), intent(in) :: qat2(n2)
      real(wp), intent(in) :: alpab(n2, n1)
      real(wp), intent(in) :: r0(n2, n1)
      integer, intent(in):: at1(n1), at2(n2)
      real(wp), intent(out) :: edr

      logical :: echo
      logical :: lp
      integer :: i, k, ii, j, jj, i0
      real(wp) :: alpha, e0
      real(wp) :: q1(n1)
      real(wp) :: q2(n2)
      real(wp) :: gdr1(3, n1)
      real(wp) :: gdr2(3, n2)
      real(wp) :: xyzdr1(3, n1)
      real(wp) :: xyzdr2(3, n2)

      xyzdr1 = xyz1
      xyzdr2 = xyz2

      alpha = 2.0

      do i = 1, n1
         q1(i) = qat1(i) - qdr1(i)
      end do
      do i = 1, n2
         q2(i) = qat2(i) - qdr2(i)
      end do

      do ii = 1, 2
         call indeg(n1, n2, nl1, nl2, npair, pair, xyz1, q1, xyzdr1, qdr1,&
              &xyz2, q2, xyzdr2, qdr2, at1, at2, alpab, r0, ii, edr, gdr1, gdr2)
         if (ii .eq. 1) then
            e0 = edr
            do i = 1, n1
               xyzdr1(1:3, i) = xyzdr1(1:3, i) - alpha*gdr1(1:3, i)
            end do
            do i = 1, n2
               xyzdr2(1:3, i) = xyzdr2(1:3, i) - alpha*gdr2(1:3, i)
            end do
         end if

         edr = edr - e0

      end do
   end subroutine drudescf

! Drude SCF gradient
   subroutine indeg(n1, n2, nl1, nl2, npair, pair, xyz1, q1, xyzdr1, qdr1,&
                   &xyz2, q2, xyzdr2, qdr2, at1, at2, alpab, r0, kode,&
                  & ec, gdr1, gdr2)
      integer, intent(in) :: n1, n2, nl1, nl2, npair, pair(2, n1*n2), kode
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyzdr1(3, n1)
      real(wp), intent(in) :: q1(n1)
      real(wp), intent(in) :: qdr1(n1)
      real(wp), intent(in) :: xyz2(3, n2)
      real(wp), intent(in) :: xyzdr2(3, n2)
      real(wp), intent(in) :: q2(n2)
      real(wp), intent(in) :: qdr2(n2)
      real(wp), intent(in) :: alpab(n2, n1)
      real(wp), intent(in) :: r0(n2, n1)
      integer, intent(in) :: at1(n1), at2(n2)
      real(wp), intent(out) :: gdr1(3, n1)
      real(wp), intent(out) :: gdr2(3, n2)
      real(wp), intent(out) :: ec

      integer :: i, j, k, i1, i2, ij
      real(wp) :: rab, rabd, r2, onedr, dum, dum2, oner
      real(wp) :: dx, dy, dz, edr, espring, qq, damp
      real(wp) :: edr1, edr2, edr3

      dum = par_drude_fc
      dum2 = 2.0d0*dum
      damp = par_drude_damp

      if (kode .eq. 1) then  ! E and G
      ! spring force and E
         espring = 0
         do i = 1, n1
            dx = (xyzdr1(1, i) - xyz1(1, i))
            dy = (xyzdr1(2, i) - xyz1(2, i))
            dz = (xyzdr1(3, i) - xyz1(3, i))
            rab = dx*dx + dy*dy + dz*dz
            dum = par_drude_fc*gam(at1(i))
            dum2 = 2.*dum
            espring = espring + dum*rab
            gdr1(1, i) = dum2*dx
            gdr1(2, i) = dum2*dy
            gdr1(3, i) = dum2*dz
         end do
         do i = 1, n2
            dx = (xyzdr2(1, i) - xyz2(1, i))
            dy = (xyzdr2(2, i) - xyz2(2, i))
            dz = (xyzdr2(3, i) - xyz2(3, i))
            rab = dx*dx + dy*dy + dz*dz
            dum = par_drude_fc*gam(at2(i))
            dum2 = 2.*dum
            espring = espring + dum*rab
            gdr2(1, i) = dum2*dx
            gdr2(2, i) = dum2*dy
            gdr2(3, i) = dum2*dz
         end do

         edr = 0
         do ij = 1, npair
            i1 = pair(1, ij)  ! on mol 1
            i2 = pair(2, ij)  ! on mol 2
            ! Coulomb force and E by other Drudes
            dx = (xyzdr1(1, i1) - xyzdr2(1, i2))
            dy = (xyzdr1(2, i1) - xyzdr2(2, i2))
            dz = (xyzdr1(3, i1) - xyzdr2(3, i2))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = qdr1(i1)*qdr2(i2)
            edr = edr + qq/rabd
            oner = qq/(rabd*rabd*rab)
            gdr1(1, i1) = gdr1(1, i1) - dx*oner
            gdr1(2, i1) = gdr1(2, i1) - dy*oner
            gdr1(3, i1) = gdr1(3, i1) - dz*oner
            gdr2(1, i2) = gdr2(1, i2) + dx*oner
            gdr2(2, i2) = gdr2(2, i2) + dy*oner
            gdr2(3, i2) = gdr2(3, i2) + dz*oner
            ! Coulomb force and E on Drudes by atoms
            dx = (xyz1(1, i1) - xyzdr2(1, i2))
            dy = (xyz1(2, i1) - xyzdr2(2, i2))
            dz = (xyz1(3, i1) - xyzdr2(3, i2))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = q1(i1)*qdr2(i2)
            edr = edr + qq/rabd
            oner = qq/(rabd*rabd*rab)
            gdr2(1, i2) = gdr2(1, i2) + dx*oner
            gdr2(2, i2) = gdr2(2, i2) + dy*oner
            gdr2(3, i2) = gdr2(3, i2) + dz*oner
            ! and reverse
            dx = (xyz2(1, i2) - xyzdr1(1, i1))
            dy = (xyz2(2, i2) - xyzdr1(2, i1))
            dz = (xyz2(3, i2) - xyzdr1(3, i1))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = q2(i2)*qdr1(i1)
            edr = edr + qq/rabd
            oner = qq/(rabd*rabd*rab)
            gdr1(1, i1) = gdr1(1, i1) + dx*oner
            gdr1(2, i1) = gdr1(2, i1) + dy*oner
            gdr1(3, i1) = gdr1(3, i1) + dz*oner
         end do

      elseif (kode .eq. 2) then ! E only
        ! spring force and E
         espring = 0
         do i = 1, n1
            dx = (xyzdr1(1, i) - xyz1(1, i))
            dy = (xyzdr1(2, i) - xyz1(2, i))
            dz = (xyzdr1(3, i) - xyz1(3, i))
            rab = dx*dx + dy*dy + dz*dz
            dum = par_drude_fc*gam(at1(i))
            espring = espring + dum*rab
         end do
         do i = 1, n2
            dx = (xyzdr2(1, i) - xyz2(1, i))
            dy = (xyzdr2(2, i) - xyz2(2, i))
            dz = (xyzdr2(3, i) - xyz2(3, i))
            rab = dx*dx + dy*dy + dz*dz
            dum = par_drude_fc*gam(at2(i))
            espring = espring + dum*rab
         end do

         edr = 0
         do ij = 1, npair
            i1 = pair(1, ij)  ! on mol 1
            i2 = pair(2, ij)  ! on mol 2
            ! Coulomb E by other Drudes
            dx = (xyzdr1(1, i1) - xyzdr2(1, i2))
            dy = (xyzdr1(2, i1) - xyzdr2(2, i2))
            dz = (xyzdr1(3, i1) - xyzdr2(3, i2))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = qdr1(i1)*qdr2(i2)
            edr = edr + qq/rabd
            ! Coulomb E on Drudes by atoms
            dx = (xyz1(1, i1) - xyzdr2(1, i2))
            dy = (xyz1(2, i1) - xyzdr2(2, i2))
            dz = (xyz1(3, i1) - xyzdr2(3, i2))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = q1(i1)*qdr2(i2)
            edr = edr + qq/rabd
            ! and reverse
            dx = (xyz2(1, i2) - xyzdr1(1, i1))
            dy = (xyz2(2, i2) - xyzdr1(2, i1))
            dz = (xyz2(3, i2) - xyzdr1(3, i1))
            r2 = dx*dx + dy*dy + dz*dz
            rab = sqrt(r2)
            ! rabd=rab+r0(i2,i1)/r2
            rabd = rab + damp
            qq = q2(i2)*qdr1(i1)
            edr = edr + qq/rabd
         end do

      end if

      ec = edr + espring

   end subroutine indeg

!************************************************************************
!* true electrostatic interaction for two atoms with squared distance
!* r2, distance r,
!* density exponents a,b
!* nuclear (valence) charges z1,z2
!* electron numbers n1,n2 (i.e. the atomic charge is q=z-n)
!* the density is assumed to be spherical
!************************************************************************
   subroutine truees(m1, m2, i, j, r, r2, den1, den2, z1, z2, n1, n2, es)
      integer, intent(in) :: i, j, m1, m2 ! atoms, atom numbers
      real(wp), intent(in) :: r, r2
      real(wp), intent(in) :: z1, z2, n1, n2  ! val charge, val number of el
      real(wp), intent(in) :: den1(2, 4, m1), den2(2, 4, m2)
      real(wp), intent(out) :: es

      real(wp) :: apb, cc, cc1, cc2, apb1, apb2, t, f02, vv, f1, f2
      real(wp) :: vv1, vv2, v1, v2, abcd, eabcd
      integer :: ip, jp, kp, lp, np
      real(wp), parameter :: pi = 3.14159265358979_wp
      real(wp), parameter :: tpi = 6.28318530717958_wp
      real(wp), parameter :: tpi25 = 34.98683665524963_wp

      ! set 1s density functions with 4 prims = np
      np = 4

      ! nuc-nuc
      es = z1*z2/r
      ! e-nuc, nuc-e
      v1 = 0
      v2 = 0
      do ip = 1, np
         do jp = 1, ip - 1
            apb1 = den1(1, ip, i) + den1(1, jp, i)
            cc1 = den1(2, ip, i)*den1(2, jp, i)
            apb2 = den2(1, ip, j) + den2(1, jp, j)
            cc2 = den2(2, ip, j)*den2(2, jp, j)
            call get_f02(apb1*r2, f02)
            v1 = v1 + f02*cc1*tpi/apb1
            call get_f02(apb2*r2, f02)
            v2 = v2 + f02*cc2*tpi/apb2
         end do
         apb1 = den1(1, ip, i) + den1(1, ip, i)
         cc1 = den1(2, ip, i)*den1(2, ip, i)
         apb2 = den2(1, ip, j) + den2(1, ip, j)
         cc2 = den2(2, ip, j)*den2(2, ip, j)
         call get_f02(apb1*r2, f02)
         v1 = v1 + 0.5*f02*cc1*tpi/apb1
         call get_f02(apb2*r2, f02)
         v2 = v2 + 0.5*f02*cc2*tpi/apb2
      end do
      v1 = 2.0d0*v1
      v2 = 2.0d0*v2
      es = es - v1*z2*n1 - v2*z1*n2
      ! e-e
      vv = 0
      do ip = 1, np
         do jp = 1, ip
            f1 = 2.0d0
            if (ip .eq. jp) f1 = 1.0d0
            apb1 = den1(1, ip, i) + den1(1, jp, i)
            cc1 = den1(2, ip, i)*den1(2, jp, i)*f1
            do kp = 1, np
               do lp = 1, kp - 1
                  apb2 = den2(1, kp, j) + den2(1, lp, j)
                  cc2 = den2(2, kp, j)*den2(2, lp, j)*2.0d0
                  eabcd = apb1*apb2
                  abcd = apb1 + apb2
                  t = r2*eabcd/abcd
                  call get_f02(t, f02)
                  vv = vv + tpi25/(eabcd*dsqrt(abcd))*f02*cc1*cc2
               end do
               apb2 = den2(1, kp, j) + den2(1, kp, j)
               cc2 = den2(2, kp, j)*den2(2, kp, j)
               eabcd = apb1*apb2
               abcd = apb1 + apb2
               t = r2*eabcd/abcd
               call get_f02(t, f02)
               vv = vv + tpi25/(eabcd*dsqrt(abcd))*f02*cc1*cc2
            end do
         end do
      end do
      es = es + n1*n2*vv

   end subroutine truees

   subroutine ovlp(n1, n2, r2ab, den1, den2, sab)
      integer, intent(in) :: n1, n2
      real(wp), intent(in) :: r2ab(n2, n1)
      real(wp), intent(in) :: den1(2, 4, n1), den2(2, 4, n2)
      real(wp), intent(out) :: sab(n2, n1)

      real(wp) :: apb, fn, s, e1, e2, c1, c2
      integer :: i, j, ip, jp, np
      real(wp), parameter :: pi = 3.14159265358979_wp

      np = 4
      sab = 0
      do i = 1, n1
         do j = 1, n2
            if (r2ab(j, i) .gt. 500.) cycle
            s = 0
            do ip = 1, np
               e1 = den1(1, ip, i)
               c1 = den1(2, ip, i)
               do jp = 1, np
                  e2 = den2(1, jp, j)
                  c2 = den2(2, jp, j)
                  apb = 1.0d0/(e1 + e2)
                  s = s + c1*c2*(pi*apb)**1.50d0*dexp(-r2ab(j, i)*e1*e2*apb)
               end do
            end do
            sab(j, i) = s
         end do
      end do

   end subroutine ovlp

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! CT corrections
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine chrgtransfer(env, n1, n2, at1, at2, q1, q2, qct1, qct2,&
                          &gab1, gab2, alpab, rab, pr, ect)
      type(TEnvironment), intent(inout) :: env
      integer, intent(in) :: n1, n2
      integer, intent(in) :: at1(n1), at2(n2)
      real(wp), intent(inout) :: q1(n1), q2(n2)
      real(wp), intent(in) :: gab1(n1, n1), gab2(n2, n2)
      real(wp), intent(in) :: qct1(n1, 2), qct2(n2, 2)
      real(wp), intent(in) :: rab(n2, n1)
      real(wp), intent(in) :: alpab(n2, n1)
      logical, intent(in) :: pr
      real(wp), intent(out) :: ect

      real(wp) :: r, a(6), u(3, 3), ee(3), jij12, jij21, qq, oner, shift, expo2
      real(wp) :: coup1, coup2, u0, u12, u21, ex, ex1, ex2, tmp, scal, expo
      real(wp) :: q01(n1), q02(n2), dq1, dq2, e, e0, qq0, emo, r0, gfac
      integer :: i, j, ati, atj

      q01 = q1
      q02 = q2

      shift = 0.37  ! orbital gap shift a la sTDA-xTB 0.45
      expo = 2.9   ! R coupling dep. 3.2
      scal = 0.02  ! off diag scaling 0.03
      expo2 = 100.  ! gap coupling dep.
      gfac = 1.0

      coup1 = 0
      coup2 = 0
      do i = 1, n1
         do j = 1, n2
            r = rab(j, i)
            if (r .gt. 20.0d0) cycle  ! cut-off
            ex = exp(-expo*(r/alpab(j, i))**0.5d0)
            coup1 = coup1 + qct1(i, 1)*qct2(j, 2)*ex            ! Fock coupling 12
            coup2 = coup2 + qct2(j, 1)*qct1(i, 2)*ex            ! Fock coupling 21
         end do
      end do

      a(1) = 0.0d0
      ex1 = exp(-expo2*(elumo(2) - ehomo(1))**2)
      if (elumo(2) - ehomo(1) .lt. 0) ex1 = 1.
      ex2 = exp(-expo2*(elumo(1) - ehomo(2))**2)
      if (elumo(1) - ehomo(2) .lt. 0) ex2 = 1.
      a(2) = -coup1*scal*100.*ex1
      a(4) = -coup2*scal*100.*ex2
      a(3) = elumo(2) - ehomo(1) + shift  ! single exc. CSF energy CT12
      a(5) = 0.0d0
      a(6) = elumo(1) - ehomo(2) + shift  ! single exc. CSF energy CT21
      call rsp(a, 3, 3, ee, u)
      i = 1
      u0 = u(1, 1)**2
      if (u(1, 2)**2 .gt. u0) i = 2
      if (u(1, 3)**2 .gt. u0) i = 3
      u0 = u(1, i)**2
      u12 = u(2, i)**2
      u21 = u(3, i)**2

      ! corrections to charges
      dq1 = 0
      e = 0
      e0 = 0
      do i = 1, n1
         q1(i) = u0*q1(i) + u12*(q1(i) + qct1(i, 1)) + u21*(q1(i) - qct1(i, 2))
         dq1 = dq1 + q1(i) - q01(i)
      end do
      call esxtb(n1, at1, gab1, q1, e)
      call esxtb(n1, at1, gab1, q01, e0)
      ect = (e - e0)*gfac
      if (pr) then
         write (env%unit, *) ' dq system 1       ', dq1
         write (env%unit, *) ' De system 1 /kcal ', 627.51*(e - e0)
      end if

      dq2 = 0
      e = 0
      e0 = 0
      do i = 1, n2
         q2(i) = u0*q2(i) + u12*(q2(i) - qct2(i, 2)) + u21*(q2(i) + qct2(i, 1))
         dq2 = dq2 + q2(i) - q02(i)
      end do
      call esxtb(n2, at2, gab2, q2, e)
      call esxtb(n2, at2, gab2, q02, e0)
      ect = ect + (e - e0)*gfac

      if (dq1 .gt. 0) then
         emo = -dq1*ehomo(1) + dq1*elumo(2)
      else
         emo = dq1*ehomo(2) - dq1*elumo(1)
      end if
      ect = ect + emo

      if (pr) then
         write (env%unit, *) ' dq system 2       ', dq2
         write (env%unit, *) ' De system 2 /kcal ', 627.51*(e - e0)
         write (env%unit, *) ' De MO       /kcal ', 627.51*emo
      end if

   end subroutine chrgtransfer

! xtb intramolecular ES
   subroutine esxtb(n, at, gab, q, es)
      integer, intent(in) :: n, at(n)
      real(wp), intent(in) :: q(n), gab(n, n)
      real(wp), intent(out) :: es

      real(wp) :: t, gam3(94), qq
      integer :: i, j
      data gam3/& ! PBE values
      &-0.02448, 0.178614, 0.194034, 0.154068, 0.173892, 0.16716, 0.156306,&
      &0.321466, 0.163314, 0.170862, 0.256128, 0.18906, 0.146310, 0.136686,&
      &0.123558, 0.122070, 0.119424, 0.115368, 0.175938, 0.152214, 0.24030,&
      &0.204654, 0.186552, 0.178824, 0.174474, 0.205038, 0.188772, 0.180462,&
      &0.180354, 0.176508, 0.158250, 0.139212, 0.131868, 0.119712, 0.115476,&
      &0.108444, 0.091032, 0.07698, 0.102642, 0.095346, 0.088266, 0.086364,&
      &0.085254, 0.088242, 0.087774, 0.088470, 0.091314, 0.090372, 0.110862,&
      &0.093588, 0.079908, 0.074082, 0.069342, 0.064638, 0.077826, 0.0590214,&
      &0.073614, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,&
      &0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,&
      &0.000000, 0.085236, 0.07782, 0.074112, 0.07206, 0.070188, 0.069246,&
      &0.072042, 0.07038, 0.072570, 0.108096, 0.084150, 0.070350, 0.06954,&
      &0.064866, 0.0596874, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,&
      &0.000000, 0.000000, 0.000000/

      es = 0
      do i = 1, n
         do j = i + 1, n
            es = es + q(i)*q(j)*gab(j, i)
         end do
      end do
      es = es*2.0d0

      t = 0
      do i = 1, n
         qq = q(i)*q(i)
         t = t + gam3(at(i))*q(i)*qq
         es = es + qq*gab(i, i)
      end do
      es = es*0.5d0 + t/3.0d0

   end subroutine esxtb

   subroutine get_f02(ARG, f02)
      real(wp), intent(in) :: ARG
      real(wp), intent(out) :: f02

      real(wp), parameter :: PI = 3.14159265358979_wp

      IF (ARG .GE. 1.0D-10) then
         f02 = 0.50D0*DSQRT(PI/ARG)*DERF(DSQRT(ARG))
         RETURN
      else
         f02 = 1.0D0 - 0.33333333333333D0*ARG
         RETURN
      end if

   end subroutine get_f02

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! RG atom as probe for possible sites in R space
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine intermole_probe(n1, nl1, at1, q1, c1, cl1, c2, l1, cn1, c6xx, e, eqp)
      integer, intent(in) :: n1, nl1
      integer, intent(in) :: at1(n1)
      integer, intent(in) :: l1(n1*10)
      real(wp), intent(in) :: c2(3)
      real(wp), intent(in) :: q1(n1)
      real(wp), intent(in) :: c1(3, n1)
      real(wp), intent(in) :: cl1(4, n1*10)
      real(wp), intent(in) :: c6xx(n1)
      real(wp), intent(in) :: cn1(n1)
      real(wp), intent(out) :: e, eqp

      !> QCG stuff
      real(wp),allocatable :: rab_solu_solv(:)

      integer :: i, j, iat, counter
      real(wp) :: ed, ep
      real(wp) :: r, r0, r2, rr, av, c6, r6, r8, r42, tmpr
      real(wp) :: R06, R08, dc6, t6, t8
      real(wp), parameter ::au = 627.509541d0
      real(wp), parameter ::autoang = 0.52917726d0

      ed = 0
      ep = 0
      eqp = 0
      counter = 0

      !> Distances of solute and screened solvent positions for QCG attractive potential
      if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
         allocate(rab_solu_solv(directedset%n), source=0.0_wp)
      end if

      do i = 1, n1 !cycle over every atom of molA
         iat = at1(i)
         r2 = (c1(1, i) - c2(1))**2 + (c1(2, i) - c2(2))**2 + (c1(3, i) - c2(3))**2
         if (r2 .gt. 5000.0d0) cycle
         r = sqrt(r2)
         eqp = eqp + q1(i)/(r + 0.1)
         r6 = r2*r2*r2
         r42 = rrab(iat, probe_atom_type)
         R06 = r0ab6(iat, probe_atom_type)
         R08 = r0ab8(iat, probe_atom_type)
         if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
            if(any(i == directedset%atoms)) then
               counter=counter + 1
               !Safe distance in Angström only, if solute atom
               ! (important if molA is solute with already added solvents)
               rab_solu_solv (counter) = r *autoang
            end if
         end if

         dc6 = 1./(r6 + R06) + r42/(r6*r2 + R08)
         ed = ed + dc6*c6xx(i)
         if(directedset%n > 0 .and. directed_type == p_atom_att) then
            ed = ed - directedset%val(i) * (exp(-0.1*((r2-2)**2)))
         end if

         ! ed= disperion energy (sum over atompairs)
         ! Pauli : fit to true IFF potential
         rr = -3.85*sqrt(r2/r42)
         ep = ep + val_e(iat)*(erf(rr) + 1.)/(r + 0.6)
         if(directedset%n > 0 .and. directed_type == p_atom_pot) then
            ep = ep + directedset%val(i) * (1/r2**pot_expo) * pot_scal !Add potential for directed docking
         end if
      end do

      do i = 1, nl1
         r2 = (c2(1) - cl1(1, i))**2 + (c2(2) - cl1(2, i))**2 + (c2(3) - cl1(3, i))**2
         if (r2 .gt. 1000.0d0) cycle
         r = sqrt(r2)
         eqp = eqp + cl1(4, i)/(r + 0.1)
      end do

      e = 8.0d0*ep - ed*0.70d0
      if(directedset%n > 0 .and. directed_type == p_atom_qcg) then
         !Energie in Hartree
         !Potential of type (a/2)*erf(-(b*r-c*b))+(a/2)
         !a=E_int solv-solv and is the maximum value of potential
         !b=slope of error function
         !c=mid point of error function
         e = e - &
             & (directedset%val(1)/2) * &
             & erf(-((directedset%expo(1) * minval(rab_solu_solv)) - (directedset%expo(2) * directedset%expo(1)))) &
             & + (directedset%val(1)/2)
      end if

      eqp = eqp*0.1  ! 0.1=probe charge, is added to e in search_nci.f90

   end subroutine intermole_probe

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! align xyz2(mol2) on xyz1(gridpoint cluster) by minimizing a dispersion type model potential
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   subroutine alignmol(n1, n2, at1, at2, xyz1, xyz2, coord)
      integer, intent(in) :: n1, n2
      integer, intent(in) :: at1(n1)
      integer, intent(in) :: at2(n2)
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)
      real(wp), intent(out) :: coord(6)

      real(wp), allocatable :: found(:, :), found2(:, :)
      real(wp) :: rr, av, sig, displ(6), e, r1(3), r2(3), dr(3), c2(3, n2), c6i
      real :: f(6)
      integer :: i, j, ii, jj, kk, icycle, maxparent, ndim, maxgen

      maxgen = 20
      maxparent = 50
      ndim = maxparent**2
      allocate (found(7, ndim), found2(7, ndim))
      found = 0
      found2 = 0
      coord = 0

      call cmadock(n1, n1, at1, xyz1, r1)

      displ = 0
      do i = 1, maxparent
         displ(1:3) = r1(1:3)
         call rand6(0.8d0, 0.10d0, 6.28d0, displ)   ! initial random sample
         found(1:6, i) = displ(1:6)
      end do

      do icycle = 1, maxgen

         kk = 0
         do i = 1, maxparent
            do j = 1, maxparent
               kk = kk + 1
               call random_number(f)
              displ(1:6) = found(1:6, i)*f(1:6) + found(1:6, j)*(1.0d0 - f(1:6))
               if (i .ne. j) call rand6(0.8d0, 1.00d0, 0.10d0, displ)   ! mutation only on childs
               call move2(n2, xyz2, c2, displ)
               e = 0
               c6i = float(at1(1))
               do ii = 1, n1
                  do jj = 1, n2
                     rr = (xyz1(1, ii) - c2(1, jj))**2&
                      &+ (xyz1(2, ii) - c2(2, jj))**2&
                      &+ (xyz1(3, ii) - c2(3, jj))**2
                     e = e - c6i*float(at2(jj))/(rr*rr*rr + 1.d+5)
                  end do
               end do
               found2(1:6, kk) = displ(1:6)
               found2(7, kk) = e
            end do
         end do

         call sort6(kk, found2)                              ! and sort for energy

         found = found2                                     ! new set

         av = sum(found(7, 1:maxparent))/float(maxparent)
         sig = 0
         do i = 1, maxparent
            sig = sig + (found(7, i) - av)**2
         end do
         sig = sqrt(sig/float(maxparent - 1))
         if (sig .lt. 0.001d0) exit
      end do

      coord(1:6) = found(1:6, 1)

      deallocate (found, found2)

   end subroutine alignmol

end module xtb_iff_iffenergy
