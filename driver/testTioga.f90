!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!>
!> Test code for Topology Independent Overset Grid Assembler (TIOGA)
!> 
!> Jay Sitaraman
!>
!> 03/05/2014
!>
program testTioga
  !
  use gridtype
  !
  implicit none
  !
  include 'mpif.h'
  !
  type(grid), target :: gr(2)
  type(grid), pointer :: g
  integer :: myid,numprocs,ierr,nsave
  integer :: itime,ntimesteps,iter,nsubiter
  real*8 :: t0
  integer :: blockid
  logical :: iclip
  integer :: ntypes
  integer :: nv1,nv2
  real*8 :: t1,t2
  integer :: i,n,m,ib,j
  real*8 :: xt(3),rnorm
  integer :: dcount,fcount
  integer, allocatable :: receptorInfo(:),inode(:)
  real*8, allocatable :: frac(:)
  !
  ! initialize mpi
  !
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  !
  ! read grid based on that generated by the
  ! strand/Cart generator
  !
  if (myid==0) write(6,*) '# tioga test on ',numprocs,' processes'
  call readGrid_cell(gr(1),myid)
  call readGrid_cell(gr(2),myid+numprocs)
  if (myid==0) write(6,*) '# tioga test : finished reading grids'
  !
  ! initialize tioga
  !
  call tioga_init_f90(mpi_comm_world)
  call mpi_barrier(mpi_comm_world,ierr)
  !
  
  ntypes=1
  nv1=6
  nv2=8
  do ib=1,2
   g=>gr(ib)
   if (g%n6 > 0)  then
    call tioga_registergrid_data(ib,g%bodytag(1),g%nv,g%x,g%iblank,g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
       ntypes,nv1,g%n6,g%ndc6)
   else if (g%n8 > 0) then
    call tioga_registergrid_data(ib,g%bodytag(1),g%nv,g%x,g%iblank,g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
       ntypes,nv2,g%n8,g%ndc8)
   endif
  enddo
  !
  ! example call to tioga_registergrid_data should be like this
  !
  !call tioga_registergrid_data(meshtag,       !< mesh tag for this partition (scalar)
  !                             nnodes,        !< number of nodes in this partition (scalar)
  !                             x,             !< coordinates (1-D real array with size=(3*nnodes))
  !                             iblank,        !< iblank (integer array with size=(nnodes)) 
  !                             nwbc,          !< number of wall boundary nodes   
  !                             nobc,          !< number of outer boundary nodes
  !                             wbcnode,       !< index of wall boundary nodes size=(nwbc)
  !                             obcnode,       !< index of overset boundary nodes size=(nobc)
  !                             ntypes,        !< number of type of cells (scalar)
  !                             nv,            !< number of vertices per cell for the first cell type 
  !                             ncells,        !< number cells of the first type 
  !                             connectivity,  !< connectivity of the first type of cells
  !                             ..,            !< number of vertices per cell for the second cell type
  !                             ..,            !< number of cells of second type
  !                             ..)            !< connectivity of the second type of cells 
  !                                            !< .. third, fourth etc
  call tioga_preprocess_grids                  !< preprocess the grids (call again if dynamic) 
  call cpu_time(t1)         
  call tioga_performconnectivity               !< determine iblanking and interpolation patterns
  call cpu_time(t2)

  call mpi_barrier(mpi_comm_world,ierr)
  if (myid==0) write(6,*) 'connectivity time=',t2-t1
 
  call cpu_time(t1)
  do ib=1,2
  g=>gr(ib)
  m=1
  do i=1,g%nv
     xt(:)=g%x(3*i-2:3*i)
     do n=1,g%nvar
       g%q(m)=(xt(1)+xt(2)+xt(3))
       m=m+1
     enddo
  enddo
  enddo

  do ib=1,2
   call tioga_getdonorcount(ib,dcount,fcount)
   allocate(receptorInfo(4*dcount))
   allocate(inode(fcount),frac(fcount))
   write(6,*) dcount,fcount
   !
   ! use this API if you need to interpolate the field variables yourself
   !
   call tioga_getdonorinfo(gr(ib)%bodytag(1),receptorInfo,inode,frac,dcount)
   !
   !> dcount = number of donors
   !> receptorInfo = {receptorProcess Id, receptor Index, receptor Block id, number of fractions}*dcount
   !> inode = indices for each receptor one group after the other
   !> frac  = weights for each receptor one group after the other
   !>
   m=0
   do i=1,dcount
    write(1000+myid,"(20(1x,I6))") ib,(inode(m+j),j=1,receptorInfo(4*i)+1),receptorInfo(4*i-3:4*i)
    m=m+(receptorInfo(4*i)+1)
   enddo 
   call flush(1000+myid)
   deallocate(receptorInfo,inode,frac)
  end do

  do ib=1,2
    g=>gr(ib)
    call tioga_registersolution(g%bodytag(1),g%q)
  enddo
  call tioga_dataupdate(gr(1)%nvar,'row')    !< update the q-variables (can be called anywhere)
                                             !< nvar = number of field variables per node
                                             !< if fields are different arrays, you can also 
                                             !< call this multiple times for each field
  call mpi_barrier(mpi_comm_world,ierr)
    call cpu_time(t2)
  if (myid==0) write(6,*) 'data update time=',t2-t1
  
  !
  ! compute error in interpolation
  ! should be machine-zero for the linear
  ! problem here
  ! 
  do ib=1,2
   g=>gr(ib)
   m=1
   do i=1,g%nv
     xt(:)=g%x(3*i-2:3*i)
     do n=1,g%nvar
       rnorm=rnorm+(g%q(m)-(xt(1)+xt(2)+xt(3)))**2
       m=m+1
     enddo
   enddo
  enddo
  if (myid==0) write(6,"(A36)") '-- Interpolation error statistics --'
  if (myid==0) write(6,"(A15)") 'ProcId    Error'
  call flush()
  call mpi_barrier(mpi_comm_world,ierr)
  write(6,"(I4,3x,E15.7)") myid,sqrt(rnorm/g%nv/g%nvar)
  call mpi_barrier(mpi_comm_world,ierr)

  call tioga_writeoutputfiles(g%nvar,'row') !< write output files, if need be

200 continue    
  call mpi_barrier(mpi_comm_world,ierr)
  call tioga_delete
  call mpi_finalize(ierr)
  !
end program testTioga
