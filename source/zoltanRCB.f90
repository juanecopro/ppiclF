!!
!! Zoltan module for ppiclF
!! 
!! author: jcolmena@anl.gov
!!
  module zoltanRCB
  use zoltan

  implicit none
 
  !! particle data for RCB
  integer :: numGlobObjs, numLocObjs
  real(8), dimension(:,:), allocatable :: part_grid
  real(8) :: grid_dx
  integer(ZOLTAN_INT), dimension(:), allocatable :: gids, iwork, iprocp
  real(Zoltan_DOUBLE), dimension(3) :: locMin, locMax
  integer(Zoltan_INT), dimension(100) :: nbparts, nbprocs
  integer(Zoltan_int) :: numnbparts, numnbprocs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Zoltan data to store in module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(Zoltan_Struct), pointer :: zz_obj
  integer(ZOLTAN_INT) :: ierr
  real(ZOLTAN_FLOAT) :: version
  integer(ZOLTAN_INT) :: myrank, ndimpart
  LOGICAL :: changes 
  INTEGER(Zoltan_INT) :: numGidEntries, numLidEntries
  INTEGER(Zoltan_INT) :: numImport, numExport
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids, exportGlobalGids 
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids, exportLocalGids
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs, exportProcs
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart, exportToPart
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains
!-----------------------------------------------------------------
    subroutine zoltanRCB_cleanup
    use zoltan

    deallocate(gids)
    deallocate(iwork)
    deallocate(iprocp)
    deallocate(part_grid)
    call Zoltan_Destroy(zz_obj)


    end subroutine zoltanRCB_cleanup
!---------------------------------------------------------------
    subroutine partitionWithRCB(comm)
    use zoltan

    implicit none

    integer(4), intent(in) :: comm
    integer(4) :: i
    integer, save :: icalld = 0

    if(icalld.eq.0)then

      ierr = Zoltan_Initialize(version)
      nullify(zz_obj)
      !! This is Fortran 90 code! hope it works
      zz_obj => Zoltan_Create(comm)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! General Zoltan Parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Param(zz_obj, "LB_METHOD", "RCB")
!      ierr = Zoltan_Set_Param(zz_obj, "IMBALANCE_TOL", "4.0")
      ierr = Zoltan_Set_param(zz_obj, "KEEP_CUTS", "TRUE")
!      ierr = Zoltan_Set_param(zz_obj, "REMAP", "0")
      ierr = Zoltan_Set_param(zz_obj, "NUM_LOCAL_PARTS", "1")
      ierr = Zoltan_Set_param(zz_obj, "RCB_RECTILINEAR_BLOCKS", "TRUE")
      ierr = Zoltan_Set_param(zz_obj, "RETURN_LISTS", "NONE")
!      ierr = Zoltan_Set_param(zz_obj, "REDUCE_DIMENSIONS", "TRUE")
      ierr = Zoltan_Set_param(zz_obj, "DEBUG_LEVEL", "0")
!      ierr = Zoltan_Set_param(zz_obj, "AVERAGE_CUTS", "TRUE")
!      ierr = Zoltan_Set_param(zz_obj, "RCB_RECOMPUTE_BOX", "TRUE")
      ierr = Zoltan_Set_param(zz_obj, "RCB_OUTPUT_LEVEL", "0")
!      ierr = Zoltan_Set_param(zz_obj, "RCB_MAX_ASPECT_RATIO", "5")
      ierr = Zoltan_Set_param(zz_obj, "RCB_REUSE", "2")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! register query functions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zoltNumObjs)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zoltNumGeom)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, zoltGeom)
      
      icalld = 1
    endif

    !!-----------------------------
    !! Use Zoltan to partition the vertices in the simple mesh.
    !!
    !! Params:
    !!     zz_obj           -- input (all remaining fields are output)
    !!     changes          -- 1 if partition was changed, 0 otherwise 
    !!     numGidEntries    -- Number of integers used for a global ID 
    !!     numLidEntries    -- Number of integers used for a local ID 
    !!     numImport        -- Number of vertices to be sent to me 
    !!     importGlobalGids -- Global IDs of vertices to be sent to me 
    !!     importLocalGids  -- Local IDs of vertices to be sent to me 
    !!     importProcs      -- Process rank for source of each incoming vertex 
    !!     importToPart     -- New part for each incoming vertex 
    !!     numExport        -- Number of vertices I must send to other processes
    !!     exportGlobalGids -- Global IDs of the vertices I must send 
    !!     exportLocalGids  -- Local IDs of the vertices I must send 
    !!     exportProcs      -- Process to which I send each of the vertices 
    !!     exportToPart     -- Part to which each vertex will belong 
    !!-----------------------------
    ierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
                               numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
                               numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)

    ierr = Zoltan_RCB_Box(zz_obj,myrank,ndimpart,locMin(1),locMin(2),locMin(3), &
                          locMax(1),locMax(2),locMax(3))

!!    write(6,*) 'myrank', myrank
!!    write(6,'(2I4)') (exportProcs(i), exportToPart(i), i=1,numExport) 
!!    write(6,*) 'NUMOBJS',myrank, numExport
!!    write(6,'(A5,I5,2F10.3)') ('Boxes', myrank, locmin(i), locmax(i), i=1,ndimpart)

    ierr = Zoltan_LB_Free_Part(importGlobalGids, importLocalGids, importProcs, importToPart)
    ierr = Zoltan_LB_Free_Part(exportGlobalGids, exportLocalGids, exportProcs, exportToPart)

    end subroutine partitionWithRCB

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGeom(data, num_gid_entries, num_lid_entries, global_id, &
                    local_id, geom_vec, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
    integer(ZOLTAN_INT), intent(in) :: global_id
    integer(ZOLTAN_INT), intent(in) :: local_id
    real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
    integer(ZOLTAN_INT), intent(out) :: ierr
    integer :: i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,ndimpart
      geom_vec(i) =  part_grid(i,local_id)
    enddo

    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine zoltGeom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumObjs(data, ierr)
    use zoltan 
    implicit none

    ! Local declarations
    INTEGER(Zoltan_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zoltNumObjs = numLocObjs
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end function zoltNumObjs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGetObjs (data, num_gid_entries, num_lid_entries, global_ids, & 
                        local_ids, wgt_dim, obj_wgts, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
    integer(ZOLTAN_INT), intent(out) :: global_ids(*)
    integer(ZOLTAN_INT), intent(out) :: local_ids(*)
    integer(ZOLTAN_INT), intent(in) :: wgt_dim 
    real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    ! local declarations
    integer :: i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i= 1, numLocObjs
      global_ids(i)= gids(i)
      local_ids(i) = i
    end do
    
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumGeom(data, ierr)
    use zoltan 
    implicit none
    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zoltNumGeom = ndimpart
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end function zoltNumGeom
  end module
