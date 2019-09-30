! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! pseudodiagonalization
! STEWART. J.J.P., CSASZAR, P., PULAY, P., J. COMP. CHEM.,3, 227, (1982)
!
! n   : total # MOs
! nocc: # occ MOs
! fmo : fock matrix in MO basis
! eig : replaced eigenvalues
! fmo : on output the new (rotated) eigenvectors in MO basis
! vector: scratch
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine pseudodiag(n,nocc,fmo,eig)
    implicit none
    !Dummy Arguments
    integer, intent(in)  :: n
    integer, intent(in)  :: nocc
    real*8, intent(inout):: fmo   (n,n)
    real*8, intent(out)  :: eig(n)
    !Stack Variables
    integer:: i
    integer:: j
    integer:: m
    integer:: aux_occ,aux_virt, loop_count
    integer:: blocksize, batch, dest,omp_get_max_threads
    integer:: i_virt
    integer:: nvirt
    integer:: i_occ
    integer :: bound
    real*8:: a
    real*8:: b
    real*8:: c
    real*8:: d
    real*8:: e
    !real*8:: alpha
    !real*8:: beta
    !Heap Variables and pointer targets
    integer, allocatable, save :: indizes(:,:)
    real*4, allocatable :: vector(:,:)   ! R4 saves 30 % due to less cache misses for 1000 atoms
    !real*8, allocatable, target :: occ(:,:)
    !real*8, allocatable, target :: virt(:,:)
    real*4,allocatable :: alphaarr(:,:)
    real*4,allocatable :: betaarr(:,:)


    allocate(vector(n,n))
    vector = 0
    nvirt = n - nocc

    allocate (alphaarr(nocc,nvirt))
    allocate (betaarr(nocc,nvirt))
    alphaarr = 0.0d0
    betaarr  = 0.0d0

    do i=1,n
        vector(i,i)=1
    enddo
    do i=1,n
        eig(i) = fmo(i,i)
    enddo

    !$omp parallel SHARED(alphaarr, betaarr, fmo) PRIVATE(a,c,d,e, i_occ, i_virt)
    !$omp do
    do i_occ = 1, nocc
        a=fmo(i_occ, i_occ)
        !DIR$ IVDEP
        do i_virt = 1, nvirt
            c=fmo(i_occ,nocc + i_virt)
            d=a-fmo(nocc + i_virt,nocc + i_virt)
            !if(abs(c/d).lt.1.d-6) cycle
            e=sign(sqrt(4.d0*c*c+d*d),d)
            alphaarr(i_occ,i_virt) = sqrt(0.5d0*(1.d0+d/e))
            betaarr(i_occ,i_virt) = -sign(sqrt(1.d0-alphaarr(i_occ,i_virt)*alphaarr(i_occ,i_virt)),c)
        end do ! End Loop over  i_virt from 1 to nvirt
    end do ! End Loop over  i_occ from 1 to nocc
   !$OMP END DO
   !$OMP END PARALLEL


    !blocksize = 28
    !blocksize = 30 * omp_get_max_threads() !for cluster
    blocksize = 6 * omp_get_max_threads() !for Desktop

    !Batched Loop
    do batch = 0, nvirt/blocksize-1
        do aux_occ = 1, nocc
            !$omp parallel SHARED(alphaarr, betaarr, vector, aux_occ, n, nocc), &
            !$OMP& PRIVATE(i_occ, i_virt, aux_virt)
            !$omp do

            !DIR$ IVDEP
            do aux_virt = 1, blocksize
                i_virt = batch*blocksize + aux_virt
                i_occ  = i_virt+aux_occ-1 - (((i_virt-1+aux_occ)/(nocc+1))*nocc)
                 !      rotation of pseudo-eigenvectors
                call srot(n, vector(:,i_occ), 1, vector(:,nocc+i_virt), 1, alphaarr(i_occ,i_virt), betaarr(i_occ,i_virt) )
            end do ! End Loop over  j from 1 to nocc
            !$OMP END DO
        !$OMP END PARALLEL
        end do ! End Loop over  i from 1 to nvirt
    end do ! End Loop over  batch from 1 to nvirt / blocksize + 1

    !Residual Loop
    do aux_occ = 1, nocc
        !$omp parallel SHARED(alphaarr, betaarr, vector, aux_occ, n, nocc), &
        !$OMP& PRIVATE(i_occ, i_virt)
        !$omp do

        !DIR$ IVDEP
        do aux_virt = 1, mod(nvirt,blocksize)
            i_virt = batch*blocksize + aux_virt
            i_occ  = i_virt+aux_occ-1 - (((i_virt-1+aux_occ)/(nocc+1))*nocc)
            !if ((alphaarr(i_occ,i_virt)==0.0) .and. (betaarr(i_occ,i_virt)==0.0)) cycle
            call srot(n, vector(:,i_occ), 1, vector(:,nocc+i_virt), 1, alphaarr(i_occ,i_virt), betaarr(i_occ,i_virt) )
        end do ! End Loop over  j from 1 to nocc
        !$OMP END DO
    !$OMP END PARALLEL
    end do ! End Loop over  i from 1 to nvirt

    fmo = vector
    deallocate(vector, alphaarr, betaarr)
end
