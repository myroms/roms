!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name: MCT_1_0_12 $ 
!BOP -------------------------------------------------------------------
!
! !ROUTINE: cpl  -- coupler for unit tester
!
! !DESCRIPTION:
! A coupler subroutine to test functionality of MCT.
!
! !INTERFACE:
!
      subroutine cpl (CPL_World)
!
! !USES:
!
      use MPH_all
!---Field Storage DataType and associated methods
      use m_AttrVect,only    : MCT_AtrVt_init => init
      use m_AttrVect,only    : MCT_AtrVt_clean => clean
      use m_AttrVect,only    : MCT_AtrVt_nreals => nRAttr
      use m_AttrVect,only    : MCT_AtrVt_lsize => lsize
      use m_AttrVect,only    : AttrVect
#ifndef SYSOSF1
      use m_AttrVect,only    : AttrVect_exportIListToChar =>exportIListToChar
      use m_AttrVect,only    : AttrVect_exportRListToChar =>exportRListToChar
#endif
!---Coordinate Grid DataType and associated methods
      use m_GeneralGrid,only: GeneralGrid
      use m_GeneralGrid,only: MCT_GGrid_clean => clean
      use m_GeneralGridComms,only: MCT_GGrid_recv => recv
      use m_GeneralGridComms,only: MCT_GGrid_scatter => scatter
!---MCT Spatial Integral services...
      use m_SpatialIntegral,only : MCT_PairedSpatialIntegrals => &
	                                      PairedSpatialIntegrals
      use m_SpatialIntegral,only : MCT_PairedSpatialAverages => &
	                                      PairedSpatialAverages
      use m_SpatialIntegral,only : MCT_PairedMaskedSpatialIntegral => &
	                                      PairedMaskedSpatialIntegrals
      use m_SpatialIntegral,only : MCT_PairedMaskedSpatialAverages => &
	                                      PairedMaskedSpatialAverages
!---Domain Decomposition Descriptor DataType and associated methods
      use m_GlobalSegMap,only: MCT_GSMap_init => init
      use m_GlobalSegMap,only: MCT_GSMap_clean => clean
      use m_GlobalSegMap,only: MCT_GSMap_gsize => gsize
      use m_GlobalSegMap,only: MCT_GSMap_lsize => lsize
      use m_GlobalSegMap,only: MCT_GSMap_ngseg => ngseg
      use m_GlobalSegMap,only: MCT_GSMap_nlseg => nlseg
      use m_GlobalSegMap,only: GlobalSegMap
      use m_GlobalSegMap,only: GlobalSegMap_bcast => bcast
!---Global-to-Local indexing services
      use m_GlobalToLocal,only: MCT_GStoL => GlobalToLocalIndices
!---Component Model Registry
      use m_MCTWorld,only: ThisMCTWorld
      use m_MCTWorld,only: MCTComponentRootRank => ComponentRootRank
      use m_MCTWorld,only: MCTWorld_init => init
      use m_MCTWorld,only: MCTWorld_clean => clean
!---Intercomponent communications scheduler
      use m_Router,only: Router
      use m_Router,only: MCT_Router_init => init
      use m_Router,only: MCT_Router_clean => clean

!---Intercomponent transfer
      use m_Transfer,only : MCT_Send => send
      use m_Transfer,only : MCT_Recv => recv

!---Sparse Matrix DataType and associated methods
      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_clean => clean
      use m_SparseMatrix, only : SparseMatrix_lsize => lsize
      use m_SparseMatrix, only : SMatrix_exportGlobalRowIndices => &
                                                        exportGlobalRowIndices
      use m_SparseMatrix, only : SMatrix_exportGlobalColumnInd => &
                                                        exportGlobalColumnIndices
      use m_SparseMatrix, only : SMatrix_exportMatrixElements => &
                                                        exportMatrixElements

      use m_SparseMatrixComms, only: SparseMatrix_ScatterByRow => ScatterByRow

!---SparseMatrixPlus DataType and associated methods
      use m_SparseMatrixPlus, only : SparseMatrixPlus
      use m_SparseMatrixPlus, only : SparseMatrixPlus_init => init
      use m_SparseMatrixPlus, only : SparseMatrixPlus_clean => clean
      use m_SparseMatrixPlus, only : Xonly ! Decompose matrix by column
      use m_SparseMatrixPlus, only : Yonly ! Decompose matrix by row
      use m_SparseMatrixPlus, only : XandY ! Arbitrary row/column decomp
!---Accumulation data type and methods
      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : MCT_Accumulator_init => init
      use m_Accumulator, only : MCT_Accumulator_clean => clean
      use m_Accumulator, only : MCT_SUM
      use m_Accumulator, only : MCT_AVG
      use m_Accumulator, only : MCT_Accumulate => accumulate
!---Matrix-Vector multiply methods
      use m_MatAttrVectMul, only: MCT_MatVecMul => sMatAvMult
!---mpeu file reading routines
      use m_inpak90
!---mpeu routines for MPI communications
      use m_mpif90               
!---mpeu timers
      use m_zeit
!---mpeu stdout/stderr
      use m_stdio
!---mpeu error handling
      use m_die

      implicit none

! !INPUT PARAMETERS:

      integer,intent(in) :: CPL_World  ! communicator for coupler

! !REVISION HISTORY:
!        Oct00 - Yun (Helen) He and Chris Ding, NERSC/LBNL - initial MPH-only version
!      19Nov00 - R. Jacob <jacob@mcs.anl.gov> -- interface with mct
!      06Feb01 - J. Larson <larson@mcs.anl.gov> - slight mod to
!                accomodate new interface to MCT_GSMap_lsize().
!      08Feb01 - R. Jacob <jacob@mcs.anl.gov> -- use MCT_Recv, new interface
!                to MCT_GSMap_lsize().
!      23Feb01 - R. Jacob <jacob@mcs.anl.gov> -- add check for transfer
!                expand size of AttrVect
!      25Feb01 - R. Jacob <jacob@mcs.anl.gov> - add mpe and mpeu
!      22Mar01 - R. Jacob <jacob@mcs.anl.gov> - use new router init
!      27Apr01 - R. Jacob <jacob@mcs.anl.gov> - use SparseMatrix
!      02May01 - R. Jacob <jacob@mcs.anl.gov> - Router is now built
!		 between atmosphere model and sparsematrix-defined
!                atmosphere globalsegmap.  Recv data in aV and check.
!                Add new argument to MCT_Smat2xGSMap.
!      16May01 - Larson/Jacob <jacob@mcs.anl.gov> - only root
!                needs to call ReadSparseMatrix with new Comms
!      17May01 - R. Jacob <jacob@mcs.anl.gov> - perfrom the sparse
!                matrix multiply on the received dummy data and check
!      19May01 - R. Jacob <jacob@mcs.anl.gov> - verify that matrix
!                multiply works on constant data
!      11Jun01 - Larson/Jacob - receive atmosphere's general grid from
!                the atmosphere.
!      15Feb02 - R. Jacob <jacob@mcs.anl.gov> New MCTWorld argument
!      28Mar02 - R. Jacob <jacob@mcs.anl.gov> Use Rearranger
!      12Jun02 - J. Larson <larson@mcs.anl.gov> - Use SparseMatrix
!                export routines.
!
!EOP ___________________________________________________________________

      character(len=*), parameter :: cplname='cpl.F90'

!----------------------- MPH vars
      integer :: myProc, myProc_global
      integer :: Global_World
      integer :: atmo_id, ocn_id
      integer :: ncomps,mycompid,mySize

!----------------------- MCT and dummy model vars

      integer :: root,stat,status
      integer, dimension(2) :: sMat_src_dims, sMat_dst_dims

!  SparseMatrix dimensions and Processor Layout
      integer :: Nax, Nay                     ! Atmosphere lons, lats
      integer :: Nox, Noy                     ! Ocean lons, lats
      integer :: NPROCS_LATA, NPROCS_LONA     ! Processor layout

!  Arrays used to initialize the MCT GlobalSegMap     
      integer :: asize,asize2,i,j,k
      integer :: osize,osize2
      integer,dimension(1) :: start,length
!     integer,dimension(:),pointer :: lstart,llength

!  Number of accumulation steps and accumulator dummy variables
      integer :: steps    
      integer, parameter :: nsteps = 10       
      character*64 :: ACCA2O_rList
      integer, dimension(:), allocatable :: ACCA2O_rAction

! Dummy arrays used for testing SparseMatrix export routines:
      integer :: Num
      integer, dimension(:), pointer :: DummyI
      real,    dimension(:), pointer :: DummyR

!  Atmosphere and Ocean GSMap
      type(GlobalSegMap) :: AGSMap,OGSMap, DAGSMap

! Router from Atm to Cpl
      type(Router)	 :: Atm2Cpl

! Router from Cpl to Ocn
      type(Router)	 :: Cpl2Ocn

! Accumulator for data from atmosphere to ocean
      type(Accumulator) :: ACCA2O

! AttrVect for data from the atm
      type(AttrVect) :: fromatm

! AttrVect for data from the atm on the ocean grid
      type(AttrVect) :: fromatm_ocn

! AttrVect for data from the ocn
      type(AttrVect) :: fromocn

! AttrVect for data from the ocn on the atmosphere's grid
      type(AttrVect) :: fromocn_atm

! AttrVects for PairedSpatialIntegral services
      type(AttrVect) :: IntegratedAVect, IntegratedOVect

! a2o SparseMatrix elements on root
      type(SparseMatrix) :: DummySMat

! a2o and o2a distributed SparseMatrixPlus variables
      type(SparseMatrixPlus) :: A2OMatPlus, O2AMatPlus

! The atmosphere's grid recieved from the atmosphere
      type(GeneralGrid) :: AtmGrid

! The atmosphere's distributed grid 
      type(GeneralGrid) :: dAtmGrid

! The ocean's grid recieved from the ocean
      type(GeneralGrid) :: OcnGrid

! The ocean's distributed grid 
      type(GeneralGrid) :: dOcnGrid

#ifdef MPE
#include "mpe.h"
#endif

!------------------------------------Begin code

  call MPI_COMM_DUP (MPI_COMM_WORLD, Global_World, ierr)

  call MPI_COMM_RANK (MPI_COMM_WORLD, myProc_global, ierr)
  call MPI_COMM_RANK (CPL_World, myProc, ierr)
! write(*,*) myProc, ' in cpl === ', myProc_global, ' in global'
! write(*,*) 'MPH_local_proc_id()=', MPH_local_proc_id_ME_SE() 
! write(*,*) 'MPH_global_proc_id()=', MPH_global_proc_id() 

  call MPI_COMM_SIZE(CPL_World,mySize,ierr)
  if (myProc==0) call MPH_redirect_output ('cpl')
  ncomps=MPH_total_components()
  mycompid=MPH_component_id_ME_SE()

! Get the atmosphere's component id
  atmo_id = MPH_get_component_id("atmosphere")

! Get the ocean's component id
  ocn_id = MPH_get_component_id("ocean")

!-------------------------------------------------------
!  Begin attempts to use MCT

#ifdef MPE
  call mpe_logging_init(myProc_global,init_s,init_e,gsmi_s,gsmi_e, &
   atri_s,atri_e,routi_s,routi_e,send_s,send_e,recv_s,recv_e, &
   clean_s,clean_e)
#endif


  if(myProc==0)write(stdout,*) cplname, ":: Initializing MCTWorld"
  call zeit_ci('Cworldinit')
   call MCTWorld_init(ncomps,MPI_COMM_WORLD,CPL_World,mycompid)
  call zeit_co('Cworldinit')

! Read in Sparse Matrix dimensions and processor layout

  if(myProc==0) then

     ! Read in SparseMatrix dimensions for atmosphere and ocean
     call I90_LoadF("ut_SparseMatrix.rc", ierr)

     call I90_Label("atmosphere_dimensions:", ierr)
     Nax = I90_GInt(ierr)
     Nay = I90_GInt(ierr)

     call I90_Label("ocean_dimensions:", ierr)
     Nox = I90_GInt(ierr)
     Noy = I90_GInt(ierr)

     call I90_Release(ierr)

     ! Read in processor layout information for atmosphere and ocean
     call I90_LoadF("./processors_map.in", ierr)

     call I90_Label("NPROCS_ATM", ierr)
     NPROCS_LATA = I90_GInt(ierr)
     NPROCS_LONA = I90_GInt(ierr)

     call I90_Release(ierr)
         
  endif

  root = MCTComponentRootRank(mycompid,ThisMCTWorld)
  call MPI_BCAST(Nax,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nay,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nox,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Noy,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NPROCS_LATA,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NPROCS_LONA,1,MP_INTEGER,root,MPI_COMM_WORLD,ierr)

!::::Receive the Atmosphere's General Grid on the root process

 if(myProc==0) then
    write(stdout,*) cplname, ":: Receiving Grid from atmosphere"

    call MCT_GGrid_recv(AtmGrid, atmo_id, 1400, status)

! check that we can make inquiries about the atmosphere's grid.
    write(stdout,*) cplname, ':: AtmGrid%coordinate_list%bf = ', &
	  AtmGrid%coordinate_list%bf
    write(stdout,*) cplname, ':: AtmGrid%index_list%bf = ', &
	  AtmGrid%index_list%bf
#ifndef SYSOSF1
    write(stdout,*) cplname, ':: AtmGrid%data%iList%bf = ', &
	 AttrVect_exportIListToChar(AtmGrid%data)
    write(stdout,*) cplname, ':: AtmGrid%data%rList%bf = ', &
	 AttrVect_exportRListToChar(AtmGrid%data)
#endif
    write(stdout,*) cplname, ':: size(AtmGrid%data%iAttr) = ', &
	  size(AtmGrid%data%iAttr)
    write(stdout,*) cplname, ':: size(AtmGrid%data%rAttr) = ', &
	  size(AtmGrid%data%rAttr)

!!!!!!!!!!!!! Receive the Ocean's General Grid
!
    write(stdout,*) cplname, ":: Receiving Grid from ocean"

    call MCT_GGrid_recv(OcnGrid, ocn_id, 2800, status)

! check that we can make inquiries about the atmosphere's grid.
    write(stdout,*) cplname, ':: OcnGrid%coordinate_list%bf = ', &
	  OcnGrid%coordinate_list%bf
    write(stdout,*) cplname, ':: OcnGrid%index_list%bf = ', &
	  OcnGrid%index_list%bf
#ifndef SYSOSF1
    write(stdout,*) cplname, ':: OcnGrid%data%iList%bf = ', &
	 AttrVect_exportIListToChar(OcnGrid%data)
    write(stdout,*) cplname, ':: OcnGrid%data%rList%bf = ', &
	 AttrVect_exportRListToChar(OcnGrid%data)
#endif
    write(stdout,*) cplname, ':: size(OcnGrid%data%iAttr) = ', &
	  size(OcnGrid%data%iAttr)
    write(stdout,*) cplname, ':: size(OcnGrid%data%rAttr) = ', &
	  size(OcnGrid%data%rAttr)
 endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set a decomposition of the atmosphere in the coupler "by hand"
! For this example, the coupler will split atmosphere points
! evenly between processors.
!
! number of local atmosphere points

  asize = (Nay * Nax)/mySize
  asize2 = asize

! (Nay *Nax)/mySize isnt an integer, give extra points to last proc.
  if(myProc == mySize - 1) then
    asize = asize + mod(Nay*Nax,mySize)
  endif

! find starting point in the numbering scheme
! numbering scheme is same as that used in atmosphere model.
  start(1) = (myProc * asize2) +1
  length(1) = asize

! write(stdout,*)myProc,asize2,asize,start(1)

! describe this information in a Global Map for the atmosphere.
  if(myProc==0)write(stdout,*) cplname, ":: Initializing AGSMap"
  call zeit_ci('Cagsmapinit')
   call MCT_GSMap_init(AGSMap,start,length,0,CPL_World,mycompid)
  call zeit_co('Cagsmapinit')

! Test GlobalSegMap_bcast:

  if(myProc==0) then

     DAGSMap%comp_id = AGSMap%comp_id
     DAGSMap%ngseg = AGSMap%ngseg
     DAGSMap%gsize = AGSMap%gsize

     allocate(DAGSMap%start(DAGSMap%ngseg),DAGSMap%length(DAGSMap%ngseg), &
	  DAGSMap%pe_loc(DAGSMap%ngseg), stat=ierr)
     if(ierr/=0) call die(cplname, "allocate(DAGSMap%start...)", ierr)

     do i=1,DAGSMap%ngseg
	DAGSMap%start(i) = AGSMap%start(i)
	DAGSMap%length(i) = AGSMap%length(i)
	DAGSMap%pe_loc(i) = AGSMap%pe_loc(i)
     end do

  endif

  call GlobalSegMap_bcast(DAGSMap, 0, CPL_World)

! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%comp_id = ",DAGSMap%comp_id
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%ngseg = ",DAGSMap%ngseg
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%gsize = ",DAGSMap%gsize
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%start(1) = ",DAGSMap%start(1)
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%start(last) = ",DAGSMap%start(DAGSMap%ngseg)
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%length(1) = ",DAGSMap%length(1)
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%length(last) = ",DAGSMap%length(DAGSMap%ngseg)
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%pe_loc(1) = ",DAGSMap%pe_loc(1)
! write(stdout,*) "CPL:: myID=",myProc," DAGSMap%pe_loc(last) = ",DAGSMap%pe_loc(DAGSMap%ngseg)

!  test some GlobalSegMap module functions
! write(*,*)myProc,'number of global segs is',MCT_GSMap_ngseg(AGSMap)
! write(*,*)myProc,'local size is',MCT_GSMap_lsize(AGSMap,CPL_World)
! write(*,*)myProc,'global size is',MCT_GSMap_gsize(AGSMap)

! if(myProc ==mySize-1) then
!   do i=0,mySize-1
!     write(*,*)myProc,'number of local segs on',i,'is',MCT_GSMap_nlseg(AGSMap,i)
!   enddo
! endif
! call MCT_GStoL(AGSMap,CPL_World,lstart,llength)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Describe OGSMap, the ocean grid decomposed in the coupler

! number of local oceanpoints
  osize = (Noy * Nox)/mySize
  osize2 = osize

! (Noy *Nox)/mySize isnt an integer, give extra points to last proc.
  if(myProc == mySize - 1) then
    osize = osize + mod(Noy*Nox,mySize)
  endif
! find starting point in the numbering scheme
! numbering scheme is same as that used in ocean model.
  start(1) = (myProc * osize2) +1
  length(1) = osize

! describe this information in a Global Map for the ocean.
  if(myProc==0)write(stdout,*) cplname, ":: Initializing OGSMap"
  call zeit_ci('Cogsmapinit')
   call MCT_GSMap_init(OGSMap,start,length,0,CPL_World,mycompid)
  call zeit_co('Cogsmapinit')

!!! test some GlobalSegMap functions
! write(*,*)myProc,'number of global segs is',MCT_GSMap_ngseg(OGSMap)
! write(*,*)myProc,'local size is',MCT_GSMap_lsize(OGSMap,CPL_World)
! write(*,*)myProc,'global size is',MCT_GSMap_gsize(OGSMap)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SparseMatrix Read
! read in the SparseMatrix elements onto root
!
! This example reads in a2o
!
  if(myProc==0)write(stdout,*) cplname, ":: Reading SparseMatrix elements"

  call zeit_ci('CsmatReadnTest')
if(myProc==0) then
! NOTE: this is a custom routine, will not be part of MCT
   call ReadSparseMatrixAsc(DummySMat,"atmosphere_to_ocean_remap_file:", &
                            sMat_src_dims, sMat_dst_dims)
! Check that the values in the SparseMatrix match the values of the
! POP grid and the Gaussian grid
   if(sMat_src_dims(1) /= Nax) call die(cplname, &
        "sMat_src_dims(1) does not match Nax")
   if(sMat_src_dims(2) /= Nay) call die(cplname, &
        "sMat_src_dims(2) does not match Nay")
   if(sMat_dst_dims(1) /= Nox) call die(cplname, &
        "sMat_dst_dims(1) does not match Nox")
   if(sMat_dst_dims(2) /= Noy) call die(cplname, &
        "sMat_dst_dims(2) does not match Noy")

   nullify(DummyI) ! let first export routine create this
   Num = SparseMatrix_lsize(DummySMat)+1
   allocate(DummyR(Num), stat=ierr) ! try this one pre-created
   if(ierr /= 0) then
      write(stderr,'(2a,i8)') cplname,':: allocate(DummyR(...) failed, ierr=',ierr
      call die(cplname)
   endif

   write(stdout,'(2a)') cplname,':: SparseMatrix export tests.  Compare with'
   call SMatrix_exportGlobalRowIndices(DummySMat, DummyI, Num)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalRowIndices(): Num=',Num
   write(stdout,'(2a,i8)') cplname,':: SparseMatrix_lsize(DummySMat)=',&
	                   SparseMatrix_lsize(DummySMat)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalRowIndices() 1st Row=',DummyI(1)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalRowIndices() last Row=',DummyI(Num)

   call SMatrix_exportGlobalColumnInd(DummySMat, DummyI, Num)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalColumnIndices(): Num=',Num
   write(stdout,'(2a,i8)') cplname,':: SparseMatrix_lsize(DummySMat)=',&
	                   SparseMatrix_lsize(DummySMat)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalColumnIndices() 1st Col=',DummyI(1)
   write(stdout,'(2a,i8)') cplname,':: exportGlobalColumnIndices() last Col=',DummyI(Num)

   call SMatrix_exportMatrixElements(DummySMat, DummyR, Num)
   write(stdout,'(2a,i8)') cplname,':: exportMatrixElements(): Num=',Num
   write(stdout,'(2a,i8)') cplname,':: SparseMatrix_lsize(DummySMat)=',&
	                   SparseMatrix_lsize(DummySMat)
   write(stdout,'(2a,f10.8)') cplname,':: exportMatrixElements() 1st wgt=',&
	DummyR(1)
   write(stdout,'(2a,f10.8)') cplname,':: exportMatrixElements() last wgt=', &
	DummyR(Num)

   deallocate(DummyI, DummyR, stat=ierr)
   if(ierr /= 0) then
      write(stderr,'(2a,i8)') cplname,':: deallocate(DummyR(...) failed, ierr=',&
	                      ierr
      call die(cplname)
   endif

endif

  call zeit_co('CsmatReadnTest')
  if(myProc==0)write(stdout,*) cplname, ":: Done Reading elements"

! put in a barrier call to wait for root
  call MPI_Barrier(CPL_World,ierr)
  root=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Build A2OMatPlus from root-centric sMat.  Specify matrix decomposition 
! to be by row, following the ocean's GlobalSegMap (OGSMap)

  call SparseMatrixPlus_init(A2OMatPlus, DummySMat, AGSMap, OGSMap, Xonly, &
                             root, CPL_World, mycompid)

! Destroy DummySMat on root:
  if(myProc==0) then
     call SparseMatrix_clean(DummySMat)
  endif

  if(myProc==0) write(stdout,*) cplname, ':: Reading in O2A on root.'

! On the root, read in O2A ascii file into DummySMat:
  if(myProc==0) then
     call ReadSparseMatrixAsc(DummySMat,"ocean_to_atmosphere_remap_file:", &
                              sMat_src_dims, sMat_dst_dims)
     if(sMat_src_dims(1) /= Nox) call die(cplname, &
          "sMat_src_dims(1) does not match Nox")
     if(sMat_src_dims(2) /= Noy) call die(cplname, &
          "sMat_src_dims(2) does not match Noy")
     if(sMat_dst_dims(1) /= Nax) call die(cplname, &
          "sMat_dst_dims(1) does not match Nax")
     if(sMat_dst_dims(2) /= Nay) call die(cplname, &
          "sMat_dst_dims(2) does not match Nay")
  endif

  if(myProc==0) write(stdout,*) cplname, ':: Finished reading in O2A on root.'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Build O2AMatPlus from root-centric sMat.  Specify matrix decomposition 
! to be by column, following the ocean's GlobalSegMap (OGSMap)

  call SparseMatrixPlus_init(O2AMatPlus, DummySMat, OGSMap, AGSMap, Yonly, &
                             root, CPL_World, mycompid)

! Destroy DummySMat on root:
  if(myProc==0) then
     call SparseMatrix_clean(DummySMat)
  endif

!!!!!!!!!!!!!!!!!----------Attribute Vector for incoming Atmosphere data
! Build an Attribute Vector to hold data coming in from Atmosphere's
! decomposition to AGSMap 
!
  if(myProc==0)write(stdout,*) cplname, ":: Initializing Attrvect"
  call zeit_ci('Catvecinit')
   call MCT_AtrVt_init(fromatm, &
       iList='gsindex',         &! local GSMap values
       rList=&
! height of first atm level
       "alevh:&
!  u wind 
       &uwind:&
!  v wind
       &vwind:&                 
!  potential temp
       &pottem:&
!  specific humidity
       &s_hum:&
!  density
       &rho:&
!  barometric pressure 
       &barpres:&
! surface pressure
       &surfp:&
!  net solar radiation
       &solrad:&
! downward direct visible radiation
       &dirvis:&
! downward diffuse visible radiation
       &difvis:&
! downward direct near-infrared radiation
       &dirnif:&
! downward diffuse near-infrared radiation
       &difnif:&
! downward longwave radiation
       &lngwv:&
! convective precip
       &precc:&
! large-scale precip
       &precl",&
       lsize=MCT_GSMap_lsize(AGSMap, Cpl_World))
  call zeit_co('Catvecinit')

!!! declare an AttrVect to hold atmosphere data on the ocean grid
! use AtrVect already declared so that it has the same Attributes
!
if(myProc==0)write(stdout,*) cplname, ":: Init output AttrVect"
  call MCT_AtrVt_init(fromatm_ocn, fromatm,MCT_GSMap_lsize(OGSMap, Cpl_World))
if(myProc==0)write(stdout,*) cplname, ":: Done with init of output vector"


!!!!!!!!!!!!!!!!!----------Attribute Vector for incoming Ocean data
! Build an Attribute Vector to hold data coming in from Ocean's Decomp
! decomposition to OGSMap 
!
  if(myProc==0)write(stdout,*)cplname,":: Initializing Incoming Ocean Attrvect"

  call zeit_ci('fromocnAVinit')

  call MCT_AtrVt_init(fromocn,    &
       rList=&
!  East-West Gradient of Ocean Surface Height
       "dhdx:&
!  North-South Gradient of Ocean Surface Height
       &dhdy:&
!  Heat of Fusion of Ocean Water
       &Qfusion:&
!  Sea Surface Temperature
       &SST:&
!  Salinity
       &salinity:&
! East Component of the Surface Current
       &Uocean:&
! East Component of the Surface Current
       &Vocean",&
       lsize=MCT_GSMap_lsize(OGSMap, CPL_World))

  call zeit_co('fromocnAVinit')

!!!!!!!!!!!!!!!!!----------Attribute Vector for Ocean data on ATM grid

  call MCT_AtrVt_init(fromocn_atm,    &
       iList="gsindex", &
       rList=&
!  East-West Gradient of Ocean Surface Height
       "dhdx:&
!  North-South Gradient of Ocean Surface Height
       &dhdy:&
!  Heat of Fusion of Ocean Water
       &Qfusion:&
!  Sea Surface Temperature
       &SST:&
!  Salinity
       &salinity:&
! East Component of the Surface Current
       &Uocean:&
! East Component of the Surface Current
       &Vocean",&
       lsize=MCT_GSMap_lsize(AGSMap, CPL_World))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--Build Router
!
! Intialize 2 routers:
! 1.) Between atmosphere and coupler using AGSMap.
! 2.) Between coupler and ocean using OGSMap

! These calls must be paired with similar calls in atm and ocn
!
  if(myProc==0)write(stdout,*) cplname, ":: Initializing Routers"

  call zeit_ci('CAtmRouterInit')
   call MCT_Router_init(atmo_id,AGSMap,CPL_World,Atm2Cpl)
  call zeit_co('CAtmRouterInit')

  call zeit_ci('COcnRouterInit')
   call MCT_Router_init(ocn_id,OGSMap,CPL_World,Cpl2Ocn)
  call zeit_co('COcnRouterInit')

  if(myProc==0)write(stdout,*) cplname, ":: Done Initializing Routers"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--Build Accumulator
  ACCA2O_rList="solrad:dirvis:difvis:dirnif:difnif:precc:precl"

  allocate(ACCA2O_rAction(7),stat=ierr)
  if(ierr/=0) call die(cplname,"allocate(ACCA20_rAction)",ierr)

  ACCA2O_rAction = (/MCT_SUM,MCT_AVG,MCT_AVG,MCT_AVG, &
                     MCT_AVG,MCT_AVG,MCT_AVG/)

  call MCT_Accumulator_init(aC=ACCA2O,          &
       rList=trim(ACCA2O_rList),                &
       rAction=ACCA2O_rAction,                  &
       lsize=MCT_GSMap_lsize(OGSMap,Cpl_World), &
       num_steps=nsteps)

  deallocate(ACCA2O_rAction,stat=ierr)
  if(ierr/=0) call die(cplname,"deallocate(ACCA20_rAction)",ierr) 

  ! Block all processes for timing purposes
  call MPI_Barrier(Global_World,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		Done with Initialization Phase
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!:::::::BEGIN REMAPPING DATA FROM ATMOSPHERE::::::::!

do steps = 1,nsteps

!!!!!!!!!!!!!!!!!----------MCT_Recv
! Receive data into AGSMap associated aV fromatm
!
if((myProc==0).and.(steps==1)) then
   write(stdout,*) cplname, ":: Doing Distributed Recv"
endif
  call zeit_ci('Cmctrecv')
   call MCT_Recv(fromatm,Atm2Cpl)
  call zeit_co('Cmctrecv')
if((myProc==0).and.(steps==1)) then
   write(stdout,*) cplname, ":: Done with Recv"
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Do the parallel A2O SparseMatrix-AttrVect multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if((myProc==0).and.(steps==1)) then
   write(stdout,*) cplname, ":: Begin A2O sparsematrix mul"
endif
  call zeit_ci('CMatMul')
   call MCT_MatVecMul(fromatm, A2OMatPlus, fromatm_ocn)
  call zeit_co('CMatMul')
if((myProc==0).and.(steps==1)) then
   write(stdout,*) cplname, ":: Completed A2O sparsematrix mul"
endif

! Perform Accumulation
call MCT_Accumulate(fromatm_ocn,ACCA2O)

enddo

 ! Send the accumulator registers to the ocean
 call zeit_ci('Cmctsend')
  call MCT_Send(ACCA2O%data,Cpl2Ocn)
 call zeit_co('Cmctsend') 

 ! Check received globalmap values against expected ones
 j=1
 do i=1,MCT_GSMap_ngseg(AGSMap)
  if(myProc==AGSMap%pe_loc(i)) then
   do k=1,AGSMap%length(i)
      if(fromatm%iAttr(1,j) /= AGSMap%start(i)+k-1) then
         write(*,*) cplname, ':: MCT GSMap mismatch. Expected', &
          AGSMap%start(i)+k-1,'got ',fromatm%iAttr(1,j)
      endif
      j=j+1
   enddo
  endif
 enddo

 ! Lets prepare to do some neat integrals using MCT.
 ! First, we scatter both of the General Grids.
 call MCT_GGrid_scatter(AtmGrid, dAtmGrid, AGSMap, 0, CPL_World)
 call MCT_GGrid_scatter(OcnGrid, dOcnGrid, OGSMap, 0, CPL_World)
! call die(cplname,"STOP HERE!")

 ! unmasked paired integral:
 call MCT_PairedSpatialIntegrals(inAv1=fromatm, outAv1=integratedAVect,    &
                                 GGrid1=dAtmGrid,WeightTag1="grid_area",   &
				 inAv2=fromatm_ocn, outAv2=integratedOVect,&
				 GGrid2=dOcnGrid, WeightTag2="grid_area",  &
				 SumWeights=.true., comm=CPL_World)
 if(myProc==0)then

    j=MCT_AtrVt_nreals(integratedAVect)
    do i=1,j,j-1
       write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired MCT ', &
	     'integral:  integratedAVect%rAttr(',i,',1)=', &
	     integratedAVect%rAttr(i,1)
    enddo

    k=MCT_AtrVt_nreals(integratedOVect)
    do i=1,k,k-1
	write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired MCT ', &
	     'integral:  integratedOVect%rAttr(',i,',1)=', &
	     integratedOVect%rAttr(i,1)
     end do
  endif

  call MCT_AtrVt_clean(integratedAVect)
  call MCT_AtrVt_clean(integratedOVect)

  ! unmasked paired average:
  call MCT_PairedSpatialAverages(inAv1=fromatm, outAv1=integratedAVect,    &
                                 GGrid1=dAtmGrid,WeightTag1="grid_area",   &
				 inAv2=fromatm_ocn, outAv2=integratedOVect,&
				 GGrid2=dOcnGrid, WeightTag2="grid_area",  &
				 comm=CPL_World)

if(myProc==0)then

   i=1
   write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired MCT ',&
	'average:  averagedAVect%rAttr(',i,',1)=', &
	integratedAVect%rAttr(i,1)

   write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired MCT ',&
	'average:  averagedOVect%rAttr(',i,',1)=', &
	integratedOVect%rAttr(i,1)

endif

 call MCT_AtrVt_clean(integratedAVect)
 call MCT_AtrVt_clean(integratedOVect)

 ! masked paired integral:
 call MCT_PairedMaskedSpatialIntegral(inAv1=fromatm, &
                           outAv1=integratedAVect,    &
                           GGrid1=dAtmGrid, &
			   SpatialWeightTag1="grid_area",   &
			   iMaskTags1="grid_imask", &
			   inAv2=fromatm_ocn, &
			   outAv2=integratedOVect, &
			   GGrid2=dOcnGrid, &
			   SpatialWeightTag2="grid_area", &
			   iMaskTags2="grid_imask", &
			   UseFastMethod=.true., &
			   SumWeights=.true., &
			   comm=CPL_World)

if(myProc==0)then

  j=MCT_AtrVt_nreals(integratedAVect)
  do i=1,j,j-1
     write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired masked MCT ',  &
     'integral: integratedAVect%rAttr(',i,',1)=', &
 	  integratedAVect%rAttr(i,1)
  end do

  k=MCT_AtrVt_nreals(integratedOVect)
  do i=1,k,k-1
     write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired masked MCT ',  &
     'integral: integratedOVect%rAttr(',i,',1)=', &
 	  integratedOVect%rAttr(i,1)
  end do

endif

 call MCT_AtrVt_clean(integratedAVect)
 call MCT_AtrVt_clean(integratedOVect)

 ! Masked paired average: 
 call MCT_PairedMaskedSpatialAverages(inAv1=fromatm, &
                           outAv1=integratedAVect,    &
                           GGrid1=dAtmGrid, &
			   SpatialWeightTag1="grid_area",   &
			   iMaskTags1="grid_imask", &
			   inAv2=fromatm_ocn, &
			   outAv2=integratedOVect, &
			   GGrid2=dOcnGrid, &
			   SpatialWeightTag2="grid_area", &
			   iMaskTags2="grid_imask", &
			   UseFastMethod=.true., &
			   comm=CPL_World)

if(myProc==0)then
   
   i=1
   write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired masked MCT ',  &
	'average : averagedAVect%rAttr(',i,',1)=', &
	integratedAVect%rAttr(i,1)

   write(stdout,'(3a,i2,a,f12.6)') cplname,':: Paired masked MCT ',  &
	'average : averagedOVect%rAttr(',i,',1)=', &
	integratedOVect%rAttr(i,1)

endif

 call MCT_AtrVt_clean(integratedAVect)
 call MCT_AtrVt_clean(integratedOVect)

 ! Now, receive Input AV from ocean (fromocn)
  if(myProc==0) write(stdout,*) cplname, ':: Before MCT_RECV from ocean'
  call zeit_ci('RecvFromOcn')
  call MCT_Recv(fromocn,Cpl2Ocn)
  call zeit_co('RecvFromOcn') 
  if(myProc==0) write(stdout,*) cplname, ':: After MCT_RECV from ocean'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Do the parallel O2A SparseMatrix-AttrVect multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(myProc==0) write(stdout,*) cplname, ":: Commencing O2A sparsematrix mul"
   call zeit_ci('O2AMatMul')
   call MCT_MatVecMul(fromocn, O2AMatPlus, fromocn_atm)
   call zeit_co('O2AMatMul')
  if(myProc==0) write(stdout,*) cplname, ":: Completed O2A sparsematrix mul"

  ! Check the interpolated values
  do i=2,MCT_AtrVt_nreals(fromocn_atm)
     do j=1,MCT_AtrVt_lsize(fromocn_atm)
        if(abs(fromocn_atm%rAttr(1,j)-fromocn_atm%rAttr(i,j)) > 1e-4) then
           write(stderr,*) cplname, ":: Interpolation Error", &
                fromocn_atm%rAttr(1,j), fromocn_atm%rAttr(i,j), i, j
           call die(cplname,"Interpolation Error")
        endif
     enddo
  enddo

if(myProc==0)write(stdout,*) cplname, ":: All Done, cleanup"
  call zeit_ci('Ccleanup')

  ! Clean MCT datatypes
  if(myProc==0) then
     call MCT_GGrid_clean(AtmGrid)
     call MCT_GGrid_clean(OcnGrid)
  endif

  call MCT_GGrid_clean(dAtmGrid)
  call MCT_GGrid_clean(dOcnGrid)
  call MCT_GSMap_clean(AGSMap)
  call MCT_GSMap_clean(OGSMap)
  call MCT_GSMap_clean(DAGSMap)
  call MCT_Router_clean(Atm2Cpl)
  call MCT_Router_clean(Cpl2Ocn)
  call SparseMatrixPlus_clean(A2OMatPlus)
  call SparseMatrixPlus_clean(O2AMatPlus)
  call MCT_Accumulator_clean(ACCA2O)
  call MCT_AtrVt_clean(fromatm)
  call MCT_AtrVt_clean(fromatm_ocn)
  call MCT_AtrVt_clean(fromocn)
  call MCT_AtrVt_clean(fromocn_atm)
  call MCTWorld_clean()

  call zeit_co('Ccleanup')

  call zeit_allflush(CPL_World,0,46)

end subroutine














