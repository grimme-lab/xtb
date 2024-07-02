module xtb_gfnff_neighbor
  use xtb_type_latticepoint, only : TLatticePoint, init_ => init_l 
  use xtb_type_molecule, only : TMolecule
  use xtb_gfnff_data, only : TGFFData
  use xtb_mctc_accuracy, only : wp
  use xtb_type_environment, only : TEnvironment
  use xtb_mctc_sort

  implicit none
  private
  
  public :: TNeigh

  !> Neighbor lists
  type :: TNeigh
    
    ! neighbor list nb(j,i,iTr); j is neighbor of i when j is shifted by iTr
    integer, allocatable :: nb(:,:,:)  
    ! full nb
    integer, allocatable :: nbf(:,:,:)
    ! no metal nb or unsually coordinated stuff
    integer, allocatable :: nbm(:,:,:)
    ! bpair(j,i,iTr) number bonds between i and j when j is translated by iTr
    integer, allocatable :: bpair(:,:,:)
    integer, allocatable :: molbpair(:)
    ! distances of all atoms to all atoms in central 27 cells
    !real(wp), allocatable :: distances(:,:,:)
    !> Maximum number of neighbors per atom
    ! where nb(numnb,i,transVec) is the num of neigh atom i has in cell transVec
    integer :: numnb = 42
    ! number of central cells so either 1,3,9,27 depending on the PBCs dimension
    integer :: numctr
    integer, allocatable :: iTrSum(:)
    ! index of translation vector pointing in opposite direction of the vector with the input index
    integer, allocatable :: iTrNeg(:)
    ! list with bool for each iTr saying if a bond to this cell should be considered
    ! this is needed to avoid counting bonds to other cells twice
    logical, allocatable :: iTrUsed(:)
    integer :: nbond
    integer :: nbond_blist
    !> bonded atoms
    integer, allocatable :: blist(:,:)
    ! number of hydrogen bonds
    integer, allocatable :: nr_hb(:)
    !> holds 3 vectors for ebond calculation: 1-shift, 2-steepness, 3-tot prefactor
    real(wp), allocatable :: vbond(:,:)
    integer :: nTrans, iTrDim
    integer, allocatable :: idTr(:) ! gives unique id of local iTr (translation vec index)
    integer, allocatable :: Trid(:) ! gives local index for iteration from unique id
    integer, allocatable :: trVecInt(:,:)  ! holds coeff for linear combination that give transVec
    real(wp), allocatable :: transVec(:,:)  ! translation vectors for generating images of unit cell
    real(wp) :: oldCutOff ! save old cutoff in getTransVec calls to quit calc if same cutoff

  contains
   
    ! initialize neigh
    procedure :: init_n
    ! get and fill the neighbor list nb
    procedure :: get_nb
    ! get the central 27 cells
   ! procedure :: getCCells
    ! get all translation Vectors
    procedure :: getTransVec
    ! get locations of nb
    procedure :: nbLoc
    ! get the j'th neighbor of atom i with corresp iTr (transVec index)
    procedure :: jth_nb
    !
    procedure :: id2v
    ! give the shell (surrounding central cell) num that the vector is in (WHichSHell: whsh)
    !!! procedure :: ws3d ! is not in type
    !return resulting vectors index if vectors indx1 and indx2 are added
    procedure :: fTrSum

    


  end type TNeigh

contains

    subroutine init_n(self, mol, env)
      class(TNeigh), intent(inout) :: self
      type(TMolecule), intent(in) :: mol
      type(TEnvironment), intent(inout) :: env
      integer :: iTrDim
      character(len=*), parameter :: source = 'neighbor'
      !
      if(mol%npbc.eq.3) then
        self%numctr = 27
      elseif(mol%npbc.eq.0) then
        self%numctr = 1
      else
        call env%warning("Periodic boundary conditions are only tested for 3 dimensional systems.", source)
      endif

      iTrDim=min(size(self%trVecInt,dim=2),729)
      self%iTrDim=iTrDim
      call filliTrSum(iTrDim, self%trVecInt, self%iTrSum, self%iTrNeg)
      call filliTrUsed(self)
      !init cutoff for getTransVec
      self%oldCutOff=0.0_wp

    end subroutine init_n


    subroutine get_nb(self, mol, rab, rtmp, mchar, icase, f_in, f2_in, param)
      class(TNeigh), intent(inout) :: self
      !> Molecular structure data         
      type(TMolecule), intent(in) :: mol
      real(wp), intent(in) :: rab(:)
      real(wp), intent(in) :: rtmp(mol%n*(mol%n+1)/2)  ! estimated bond lengths
      real(wp), intent(in) :: mchar(mol%n), f_in, f2_in
      integer, intent(in) :: icase
      !integer, intent(inout) :: nbf(20,mol%n)
      type(TGFFData), intent(in) :: param
      !real(wp) :: latThresh, maxNBr
      real(wp) :: dist(mol%n,mol%n,self%numctr)
      integer :: i,j,iTr
      
      dist=0.0_wp
      !$omp parallel do collapse(3) default(none) shared(dist,mol,self) &
      !$omp private(iTr,i,j)
      do iTr=1, self%numctr
       do i=1, mol%n
         do j=1, mol%n
           dist(j,i,iTr) = NORM2(mol%xyz(:,i)-(mol%xyz(:,j)+self%transVec(:,iTr)))
         enddo
       enddo
      enddo
      !$omp end parallel do
      call fillnb(self,mol%n,mol%at,rtmp,dist,mchar,icase,f_in,f2_in,param)

    end subroutine get_nb


    ! get number of transVec within cutoff -> neigh%nTrans
    ! get all linear comb coeff within 60 bohr (maxCutOff)-> neigh%trVecInt
    ! get vector with indices of the nTrans translation vectors -> neigh%idTr
    ! in periodic case at least central 27 transVec to avoid artifacts for cells > 60 bohr
    subroutine getTransVec(self,env,mol,cutoff)
      !> Instance of the lattice point generator                                    
      class(TNeigh), intent(inout) :: self     
      type(TEnvironment), intent(inout) :: env
      !> Molecular structure data         
      type(TMolecule), intent(in) :: mol
      !> cutoff
      real(wp), intent(in) :: cutoff
      
      integer :: ix,iy,iz,i, nLat, nShell, d,indx,id
      integer, allocatable :: idTrTmp(:)
      real(wp), allocatable :: sVec(:) ! shortest lattice vector
      real(wp), allocatable :: vec(:), trVecTmp(:,:)
      real(wp) :: sh
      character(len=*), parameter :: source = 'neighbor'
      ! check if last calculation was with same cutoff
      if(self%oldCutOff.eq.cutoff) then
        return
      else
        self%oldCutoff=cutoff
      endif  

      ! seperate treatment depending on periodicity
      if (mol%npbc.eq.0) then !#ifmolnpbc
        if(.not.allocated(self%idTr)) allocate(self%idTr(1))
        if(.not.allocated(self%Trid)) allocate(self%Trid(1))
        self%nTrans=1
        self%idTr=1
        self%Trid=1
        if(.not.allocated(self%trVecInt)) allocate(self%trVecInt(3,1),source=0)
        if(.not.allocated(self%transVec)) allocate(self%transVec(3,1),source=0.0_wp)

      ! periodic treatment  
      else  !#ifmolnpbc
        nLat=mol%npbc
        ! find shortest lattice vector
        allocate(sVec(3))
        sVec = 1.0d15
        do i=1, 3 ! go through all lattice vectors
          if(.not.mol%pbc(i)) cycle
          if(NORM2(mol%lattice(:,i)).lt.NORM2(sVec)) then
            sVec=mol%lattice(:,i)
          endif
        enddo
        if(NORM2(sVec).eq.0.0_wp) then
          call env%error('A lattice vector for a periodic axis is zero.', source)
          return
        endif
        ! see how often sVec fits in cutoff -> gives number of needed shells
        ! e.g. in 3D  nShell=1 -> 3x3x3    nShell=4 -> 9x9x9
        nShell = CEILING(cutoff/NORM2(sVec))
        d=(2*nShell+1)**nLat

        ! lattice might change during optimization. Therefore d might change
        !  between different runs and trVecInt needs to be adjusted.
        if(allocated(self%trVecInt)) then
          if(size(self%trVecInt,dim=2).lt.d) then
            deallocate(self%trVecInt)
          endif
        endif

        ! get the linear combination coefficient vector
        ! different loops depending on dimension (nLat)
        if(.not.allocated(self%trVecInt)) then !#ifnotallocatedtrvecint
          if(nLat.eq.3) then
            allocate(self%trVecInt(3,d),source=0)
            indx=1
            do i=1,nShell
              do iz=-i, i
                do iy=-i, i
                  do ix=-i, i
                    !
                    if(abs(ix).ne.i.and.abs(iy).ne.i.and.abs(iz).ne.i) cycle
                    indx=indx+1
                    self%trVecInt(:,indx)=[ix,iy,iz]
                  enddo
                enddo
              enddo  
            enddo

          elseif(nLat.eq.2) then
            allocate(self%trVecInt(3,d),source=0)
            indx=1
            do i=1,nShell
              do iz=-i, i
                do iy=-i, i
                  do ix=-i, i
                    ! one of the lattice vecs is not periodic, only cycle wrt the other vecs
                    if(.not.mol%pbc(1)) then
                      if(abs(iy).ne.i.and.abs(iz).ne.i) cycle
                    elseif(.not.mol%pbc(2))then
                      if(abs(ix).ne.i.and.abs(iz).ne.i) cycle
                    elseif(.not.mol%pbc(3))then
                      if(abs(ix).ne.i.and.abs(iy).ne.i) cycle
                    endif
                    indx=indx+1
                    if(mol%pbc(1)) self%trVecInt(1,indx)= ix
                    if(mol%pbc(2)) self%trVecInt(2,indx)= iy
                    if(mol%pbc(3)) self%trVecInt(3,indx)= iz
                  enddo
                enddo
              enddo  
            enddo

          elseif(nLat.eq.1) then
            allocate(self%trVecInt(3,d),source=0)
            indx=1
            do i=1,nShell
              do iz=-i, i
                do iy=-i, i
                  do ix=-i, i
                    ! only one of the lattice vecs is periodic, only cycle wrt this vector
                    if(mol%pbc(1)) then
                      if(abs(ix).ne.i) cycle
                    elseif(.not.mol%pbc(2))then
                      if(abs(iy).ne.i) cycle
                    elseif(.not.mol%pbc(3))then
                      if(abs(iz).ne.i) cycle
                    endif
                    indx=indx+1
                    if(mol%pbc(1)) self%trVecInt(1,indx)= ix
                    if(mol%pbc(2)) self%trVecInt(2,indx)= iy
                    if(mol%pbc(3)) self%trVecInt(3,indx)= iz
                  enddo
                enddo
              enddo  
            enddo

          else
            call env%error('Dimension of PBCs could not be processed.', source)
            return
          endif

        endif !#ifnotallocatedtrvecint
        ! Now the linear combinations of lattice vectors are given with neigh%trVecInt
        ! Next up: get number of translation vectors within cutoff
        if(allocated(self%Trid)) deallocate(self%Trid)
        allocate(self%Trid(d),source=0)
        allocate(idTrTmp(d),source=0)
        allocate(vec(nLat),source=0.0_wp)
        indx=0

        if(nLat.eq.3) then
          allocate(trVecTmp(3,d))
          do id=1, d
            !
            vec = self%trVecInt(1,id)*mol%lattice(:,1) &
                + self%trVecInt(2,id)*mol%lattice(:,2) &
                + self%trVecInt(3,id)*mol%lattice(:,3)
            if(id.gt.27)then ! always want the central 27 cells
              ! check which shell the transVec is from (they are ordered by shell)  
              if(id.eq.28)   sh = 2.0  ! first shell has 27 cells
              if(id.eq.126)  sh = 3.0  ! second has   125
              if(id.eq.344)  sh = 4.0  ! third  has   343
              if(id.eq.730)  sh = 5.0  ! fourth has   729
              if(id.eq.1332) sh = 6.0  ! fifth  has  1331
              if(id.eq.2198) sh = 7.0  ! sixth  has  2197
              if(id.eq.3376) sh = 8.0  !   7th  has  3375
              if(id.eq.4914) sh = 9.0  !   8th  has  4913
              if(id.eq.6860) sh = 10.0 !   9th  has  6859
              if(id.eq.9262) sh = 11.0  ! 10th  has  9261
              if(id.eq.12168) sh = 12.0 ! 11th  has 12167
              if(id.eq.15626) sh = 13.0 ! 12th  has 15625
              if(id.eq.19684) sh = 14.0 ! 13th  has 19683
              ! the factor ensures that the cutoff is completely covered
              if(NORM2(vec)*(1-1/sh).gt.cutoff) then
              cycle
              endif
            endif
            indx=indx+1
            idTrTmp(indx)=id ! id is ALWAYS the same for one specific linear combination (LC)
            self%Trid(id)=indx !-> if you put indx into trVecInt you will always get the same LC
            trVecTmp(:,indx)=vec
          enddo
          if(id.gt.24389) then
            call env%error('Implementation is not adjusted for such small cells.', source)
            return
          endif
          if(allocated(self%transVec)) deallocate(self%transVec)
          if(allocated(self%idTr)) deallocate(self%idTr)
          allocate(self%transVec(3,indx))
          allocate(self%idTr(indx),source=0)
          self%transVec=trVecTmp(:,1:indx)
          self%idTr=idTrTmp(1:indx)

        elseif(nLat.eq.2)then
          !
          allocate(trVecTmp(2,d))
          do id=1, d
            !
            vec = self%trVecInt(1,id)*mol%lattice(:,1) &
                + self%trVecInt(2,id)*mol%lattice(:,2)
           if(id.gt.9)then ! always want the central 9 cells (2-Dim !)
              if(id.eq.10) sh = 2.0
              if(NORM2(vec).gt.cutoff) cycle
            endif
            indx=indx+1
            self%idTr(indx)=id ! assure that indx is ALWAYS the same
                               ! for one specific linear combination (LC)
                               ! that means: if you put indx into trVecInt
                               ! you will always get the same LC
            trVecTmp(:,indx)=vec
          enddo
          if(allocated(self%transVec)) deallocate(self%transVec)
          allocate(self%transVec(2,indx))
          self%transVec=trVecTmp(:,1:indx)

        elseif(nLat.eq.1)then
          !      
          allocate(trVecTmp(1,d))
          do id=1, d
            !
            vec = self%trVecInt(1,id)*mol%lattice(:,1)
            if(indx.gt.27)then ! always want the central 27 cells
              if(NORM2(vec).gt.cutoff) cycle
            endif
            indx=indx+1
            self%idTr(indx)=id ! assure that indx is ALWAYS the same
                               ! for one specific linear combination (LC)
                               ! that means: if you put indx into trVecInt 
                               ! you will always get the same LC
            trVecTmp(:,indx)=vec
          enddo
          if(allocated(self%transVec)) deallocate(self%transVec)
          allocate(self%transVec(1,indx))
          self%transVec=trVecTmp(:,1:indx)
        else
          call env%error('Given lattice could not be processed.', source)
          return
        endif

      endif !#ifmolnpbc
      self%nTrans=size(self%transVec,dim=2)

    end subroutine getTransVec

    
    ! returns translation vector to unique ID, needs neigh%idTr and neigh%transVec
    function id2v(self,mol,id) result(vector)
      !> Instance of the lattice point generator                                    
      class(TNeigh), intent(inout) :: self
      !> Molecular structure data         
      type(TMolecule), intent(in) :: mol
      ! unique ID of the translation vector
      integer, intent(in) :: id
      ! index in transVec
      integer :: indx
      ! result: corresponding translation vector
      real(wp) :: vector(3)

      if(id.le.size(self%trVecInt,dim=2)) then
        !
        vector = self%trVecInt(1,id)*mol%lattice(:,1) &
               + self%trVecInt(2,id)*mol%lattice(:,2) &
               + self%trVecInt(3,id)*mol%lattice(:,3)  

      endif


    end function id2v


    ! is a modified "getnb" from gfnff_ini2.f90
    subroutine fillnb(self,n,at,rad,dist,mchar,icase,f,f2,param)
      implicit none
      type(TGFFData), intent(in) :: param
      class(TNeigh), intent(inout) :: self
      integer, intent(in) :: n,at(n)
      real(wp), intent(in) :: rad(n*(n+1)/2)
      real(wp), intent(in) :: dist(n,n,self%numctr)
      real(wp), intent(in) :: mchar(n),f,f2
      integer i,j,k,iTr,nn,icase,hc_crit,nnfi,nnfj,lin
      integer tag(n,n,self%numctr)
      real(wp) rco,fm
      if(icase.eq.1) then
        if(.not. allocated(self%nbf)) then
          allocate(self%nbf(self%numnb,n,self%numctr),source=0)
        else
          self%nbf=0
        endif
      endif
      if(icase.eq.2) then
        if(.not. allocated(self%nb)) then
          allocate(self%nb(self%numnb,n,self%numctr),source=0)
        else
          self%nb=0
        endif
      endif
      if(icase.eq.3) then
        if(.not. allocated(self%nbm)) then
          allocate(self%nbm(self%numnb,n,self%numctr),source=0)
        else
          self%nbm=0
        endif
      endif

      tag = 0
      do iTr=1, self%numctr
        do i=1,n
          do j=1,n
        
            k=lin(j,i)
            if (dist(j,i,iTr).le.0.0_wp) cycle
            fm=1.0d0
!           full case                                                           
            if(icase.eq.1)then
              if(param%metal(at(i)).eq.2) fm=fm*f2 !change radius for metal atoms
              if(param%metal(at(j)).eq.2) fm=fm*f2                             
              if(param%metal(at(i)).eq.1) fm=fm*(f2+0.025)                     
              if(param%metal(at(j)).eq.1) fm=fm*(f2+0.025)                     
            endif                                                               
!           no HC atoms
            if(icase.eq.2)then
               nnfi=sum(self%nbf(self%numnb,i,:))               ! full CN of i, only valid for icase > 1
               nnfj=sum(self%nbf(self%numnb,j,:))               ! full CN of i
               hc_crit = 6
               if(param%group(at(i)).le.2) hc_crit = 4
               if(nnfi.gt.hc_crit) then
                       cycle
                endif
               hc_crit = 6
               if(param%group(at(j)).le.2) hc_crit = 4
               if(nnfj.gt.hc_crit) then
             cycle  
             endif
            endif
!           no metals and unusually coordinated stuff                                               
            if(icase.eq.3)then
               nnfi=sum(self%nbf(self%numnb,i,:))               ! full CN of i, only valid for icase > 1
               nnfj=sum(self%nbf(self%numnb,j,:))               ! full CN of i
               if(mchar(i).gt.0.25 .or. param%metal(at(i)).gt.0) cycle   ! metal case TMonly
               if(mchar(j).gt.0.25 .or. param%metal(at(j)).gt.0) cycle   ! metal case
               if(nnfi.gt.param%normcn(at(i)).and.at(i).gt.10)   cycle   ! HC case
               if(nnfj.gt.param%normcn(at(j)).and.at(j).gt.10)   cycle   ! HC case
            endif

            rco=rad(k)

            if(dist(j,i,iTr).lt. fm * f * rco) tag(j,i,iTr)=1
          enddo
        enddo
      enddo

      do iTr=1, self%numctr
      do i=1,n
         nn = 0
         do j=1,n
            if(tag(j,i,iTr).eq.1) then ! .AND.i.ne.j .OR. tag(j,i,iTr).eq.1.AND.iTr.ne.1)
               nn=nn+1
               if (icase.eq.1) self%nbf(nn,i,iTr)=j
               if (icase.eq.2) self%nb(nn,i,iTr) =j
               if (icase.eq.3) self%nbm(nn,i,iTr)=j
              ! if (iTr.eq.1) nb(nn,j,iTr)=i
            endif
         enddo
         if (icase.eq.1) self%nbf(self%numnb,i,iTr)=nn
         if (icase.eq.2) self%nb(self%numnb,i,iTr) =nn
         if (icase.eq.3) self%nbm(self%numnb,i,iTr)=nn
      enddo
      enddo

    end subroutine fillnb


    ! size(locarr, dim=2) is the number of cells that contain neighbors of indexi
    ! the other dim contains the indices of those neighbors 
    ! the LastEntry is the index (iTr) of the cell that the neighbors are in
    ! the SecondLastEntry is the number of neighbors indexi has in this cell
    subroutine nbLoc(self, n, nblist, indexi, locarr)
      class(TNeigh), intent(in) :: self
      integer, intent(in) :: n
      integer, intent(in) :: nblist(self%numnb, n, self%numctr)
      integer, intent(in) :: indexi
      integer, allocatable, intent(out) :: locarr(:,:)
      integer :: j,k,iTr,counter,counter2
      counter=0
      do iTr=1, self%numctr
       if (nblist(self%numnb,indexi,iTr).ne.0) counter = counter + 1
      enddo
      allocate(locarr(self%numnb, counter),source=0)
      ! find out which cells the neighbors are in  
      counter2=0
      do iTr=1, self%numctr
        if(nblist(self%numnb,indexi,iTr).eq.0) cycle
        counter=0           
        counter2=counter2+1
        do j=1, self%numnb-2                                  
          if(nblist(j, indexi, iTr).ne.0) then
            counter=counter+1
            locarr(counter,counter2) = nblist(j, indexi, iTr)
            locarr(self%numnb, counter2) = iTr                             ! LastEntry
            locarr(self%numnb-1, counter2) = nblist(self%numnb,indexi,iTr) ! SecondLastEntry
          endif
        enddo                                      
      enddo                                        

    endsubroutine nbLoc


    ! gives index j of the jth neighbor of i and the corresp transVec index iTr
    subroutine jth_nb(self, n, xyz, j, jth, i, iTr)
      class(TNeigh), intent(inout) :: self
      integer, intent(in) :: n
      real(wp), intent(in) :: xyz(3,n)
      integer, intent(inout) :: j, iTr
      integer, intent(in) :: jth, i
      integer ::  k, l, iTrdum, sizenb, sizeindx, jdum
      integer, allocatable :: indx(:)
      real(wp), allocatable :: distnb(:)

      sizenb = sum(self%nb(self%numnb,i,:))
      allocate(distnb(sizenb), source=0.0_wp)
      allocate(indx(sizenb), source=0)
      ! get distances from i to neighbors k in different cells iTr
      j=0
      iTr=0
      l=0
      do iTrdum=1, self%numctr
        do k=1, self%nb(self%numnb,i,iTrdum)
          jdum = self%nb(k,i,iTrdum)
          l=l+1
          distnb(l) = NORM2(xyz(:,i)-(xyz(:,jdum)+self%transVec(:,iTrDum)))
        enddo
      enddo

      ! get sorting index with heapSort from mctc_sort
      call indexHeapSort(indx, distnb)
      l=0
      do iTrdum=1, self%numctr
        do k=1, self%nb(self%numnb,i,iTrdum)
          jdum = self%nb(k,i,iTrdum)
          l=l+1
          if (l.eq.indx(jth)) then
            iTr = iTrdum
            j = jdum 
          endif
        enddo
      enddo

    end subroutine jth_nb


    function vecNorm(vector) result(vectorNorm)
      real(wp), intent(in) :: vector(3)
      real(wp) :: vectorNorm
      vectorNorm = sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
    end function vecNorm


    subroutine filliTrSum(iTrDim,trVecInt,iTrSum,iTrNeg)
       !class(TNeigh), intent(inout) :: neigh                                     
       integer, intent(in) :: iTrDim
       integer, intent(in) :: trVecInt(3,iTrDim) 
       integer, allocatable, intent(out) :: iTrSum(:), iTrNeg(:)
       integer :: trVecInt2(3,iTrDim) 
       integer :: i,j,s,k,kk, lin
       integer :: vec1(3), vec2(3), vsum(3)                                      

       allocate(iTrSum(iTrDim*(iTrDim+1)/2),source=-1) 
       allocate(iTrNeg(iTrDim),source=-1)
       if (iTrDim.eq.1) then                                               
       else                                                                                
         !                               ! 
         do i=1, iTrDim                                                          
           kk=i*(i-1)/2                                                          
           do j=1, i                                                             
             k=kk+j ! packed index                                                   
             vec1 = trVecInt(:,i)   
             vec2 = trVecInt(:,j)                                                          
             vsum = vec1 + vec2
             do s=1, iTrDim 
               !                                                                           
               if (j.eq.1) then                                                          
                 if (vec1(1).eq.-trVecInt(1,s).AND.vec1(2).eq.-trVecInt(2,s)&              
                                      &       .AND.vec1(3).eq.-trVecInt(3,s)) then   
                   iTrNeg(i) = s                                                 
                   iTrNeg(s) = i   !this might be difficult for shared?
                 endif                                                           
               endif  

               if (vsum(1).eq.trVecInt(1,s).AND.vsum(2).eq.trVecInt(2,s)&                  
                                         & .AND.vsum(3).eq.trVecInt(3,s)) then
                 iTrSum(k) = s                                                   
               endif                                                                     
             enddo                                                               
           enddo                                                                 
         enddo                                                                   
       endif                                                                     
       iTrNeg(1)=1                                                         
       iTrSum(1)=1                                                         
                                                                             
    end subroutine filliTrSum            


    function fTrSum(self,indx1,indx2) result(indxSum)
      class(TNeigh) :: self
      integer :: indx1, indx2
      integer :: indxSum
      integer :: lin, id1, id2, k, idOfSum

      ! if fTrSum=-1 then the atom is shifted out of the cutoff used
      !  in the last neigh%getTransVec call
      indxSum=-1  ! default

      ! check if input is viable
      if(abs(indx1).gt.self%nTrans.or.abs(indx2).gt.self%nTrans) then
        return
      endif 
      ! change to unique index
      id1 = self%idTr(indx1)
      id2 = self%idTr(indx2)

      ! get entry from iTrSum and revert to input indx format
      k = lin(id1,id2)
      if(k.le.size(self%iTrSum)) then
        idOfSum = self%iTrSum(k)  ! unique id of sum
        ! only assign index if it can be handled by neigh%transVec
        ! else default of -1 is used, which is checked before use of %transVec
        if(idOfSum.gt.size(self%Trid)) return
        if(idOfSum.le.0) return
        if(self%Trid(idOfSum).gt.0.and.self%Trid(idOfSum).le.self%nTrans) then
          indxSum = self%Trid(idOfSum) 
        endif
      endif

    end function fTrSum

    subroutine filliTrUsed(neigh)
      class(TNeigh), intent(inout) :: neigh
      integer :: iTr,iTrF
      !
      allocate(neigh%iTrUsed(neigh%numctr))
      neigh%iTrUsed=.true.
      do iTr=neigh%numctr, 2, -1 ! iTr=1 is central cell (0,0,0) which should be used in any case
        if(neigh%iTrUsed(iTr)) then
          iTrF=neigh%iTrNeg(iTr)
          neigh%iTrUsed(iTrF)=.false.
        endif
      enddo

    end subroutine filliTrUsed


end module xtb_gfnff_neighbor

