!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
#else
      subroutine ppiclf_comm_InitMPI(comm,id,np)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 comm
      integer*4 id
      integer*4 np
!
      if (PPICLF_LINIT .or. PPICLF_LFILT .or. PPICLF_OVERLAP)
     >   call ppiclf_exittr('InitMPI must be called first$',0.0d0,0)

      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      call ppiclf_prints('   *Begin InitCrystal$')
         call ppiclf_comm_InitCrystal
      call ppiclf_prints('    End InitCrystal$')

      PPICLF_LCOMM = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitCrystal
!
      implicit none
!
      include "PPICLF"
!
      call pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateBin
#if PPICLF_ZOLTAN==1
      use zoltanRCB
#endif
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer*4 exit_1_array(3), exit_2_array(3), finished(3)
      integer*4 ppiclf_n_bins_tgt(3)
      common /tgtbins/ ppiclf_n_bins_tgt
      integer*4 ix, iy, iz, iperiodicx, iperiodicy, iperiodicz, 
     >          npt_total, j, i, idum, jdum, kdum, total_bin, 
     >          sum_value, count
      real*8 xmin, ymin, zmin, xmax, ymax, zmax, rduml, rdumr, rthresh,
     >       rmiddle, rdiff
      logical exit_1, exit_2
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      real*8 ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
      external ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
      integer icalld, istart
      save icalld
      data icalld /0/
!

! face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq. 3)
     >iz = 3

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
         
      ppiclf_d2chk(1) = max(ppiclf_d2chk(2),ppiclf_d2chk(3))

      ! binning requires > 1 global particle. This takes care of 
      ! single particle case
      npt_total = ppiclf_iglsum(ppiclf_npart,1)
c     if (npt_total .eq. 1) then
      if (.not. ppiclf_lproj .and. .not. ppiclf_lsubsubbin) 
     >   ppiclf_d2chk(1) = 1E-16

      ! compute binb
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      do i=1,ppiclf_npart
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(1)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(1)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         if (ppiclf_ndim .eq. 3) then
            rduml = ppiclf_y(iz,i) - ppiclf_d2chk(1)
            rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(1)
            if (rduml .lt. zmin) zmin = rduml
            if (rdumr .gt. zmax) zmax = rdumr
         endif
      enddo

      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = 0.0d0
      ppiclf_binb(6) = 0.0d0
      if(ppiclf_ndim .gt. 2) ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      if(ppiclf_ndim .gt. 2) ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      if (npt_total .gt. 0) then
      do i=1,ppiclf_ndim
         if (ppiclf_bins_balance(i) .eq. 1) then
            rmiddle = 0.0
            do j=1,ppiclf_npart
               rmiddle = rmiddle + ppiclf_y(i,j)
            enddo
            rmiddle = ppiclf_glsum(rmiddle,1)
            rmiddle = rmiddle/npt_total

            rdiff =  max(abs(rmiddle-ppiclf_binb(2*(i-1)+1)),
     >                   abs(ppiclf_binb(2*(i-1)+2)-rmiddle))
            ppiclf_binb(2*(i-1)+1) = rmiddle - rdiff
            ppiclf_binb(2*(i-1)+2) = rmiddle + rdiff
         endif
      enddo
      endif

      if (ppiclf_xdrange(2,1) .lt. ppiclf_binb(2) .or.
     >    ppiclf_xdrange(1,1) .gt. ppiclf_binb(1) .or. 
     >    iperiodicx .eq. 0) then
         ppiclf_binb(1) = ppiclf_xdrange(1,1)
         ppiclf_binb(2) = ppiclf_xdrange(2,1)
      endif

      if (ppiclf_xdrange(2,2) .lt. ppiclf_binb(4) .or.
     >    ppiclf_xdrange(1,2) .gt. ppiclf_binb(3) .or.
     >    iperiodicy .eq. 0) then
         ppiclf_binb(3) = ppiclf_xdrange(1,2)
         ppiclf_binb(4) = ppiclf_xdrange(2,2)
      endif

      if (ppiclf_ndim .gt. 2) then
      if (ppiclf_xdrange(2,3) .lt. ppiclf_binb(6) .or.
     >    ppiclf_xdrange(1,3) .gt. ppiclf_binb(5) .or. 
     >    iperiodicz .eq. 0) then
         ppiclf_binb(5) = ppiclf_xdrange(1,3)
         ppiclf_binb(6) = ppiclf_xdrange(2,3)
      endif
      endif

#if PPICLF_ZOLTAN==1
      if(icalld .eq.0) then
        allocate(gids(PPICLF_LPART))
        allocate(iprocp(ppiclf_np))
        allocate(iwork(ppiclf_np))
        allocate(part_grid(ppiclf_ndim,PPICLF_LPART))
        myrank = ppiclf_nid
        ndimpart = ppiclf_ndim
        icalld=1
      endif
#endif

      if (npt_total .lt. 1) return

      finished(1) = 0
      finished(2) = 0
      finished(3) = 0
      total_bin = 1 

      do i=1,ppiclf_ndim
         finished(i) = 0
         exit_1_array(i) = ppiclf_bins_set(i)
         exit_2_array(i) = 0
         if (ppiclf_bins_set(i) .ne. 1) then
            ppiclf_n_bins(i) = 1
         else
            ppiclf_n_bins(i) = ppiclf_n_bins_tgt(i)
         endif
         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                        ppiclf_binb(2*(i-1)+1)  ) / 
     >                       ppiclf_n_bins(i)
         ! Make sure exit_2 is not violated by user input
         if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
            do while (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1))
               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                              ppiclf_binb(2*(i-1)+1)  ) / 
     >                             ppiclf_n_bins(i)
            enddo
         endif
         total_bin = total_bin*ppiclf_n_bins(i)
      enddo

      ! Make sure exit_1 is not violated by user input
      count = 0
      do while (total_bin > ppiclf_np)
          count = count + 1;
          i = modulo((ppiclf_ndim-1)+count,ppiclf_ndim)+1
          ppiclf_n_bins(i) = max(ppiclf_n_bins(i)-1,1)
          ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                         ppiclf_binb(2*(i-1)+1)  ) / 
     >                        ppiclf_n_bins(i)
          total_bin = 1
          do j=1,ppiclf_ndim
             total_bin = total_bin*ppiclf_n_bins(j)
          enddo
          if (total_bin .le. ppiclf_np) exit
       enddo

       exit_1 = .false.
       exit_2 = .false.

       do while (.not. exit_1 .and. .not. exit_2)
          do i=1,ppiclf_ndim
             if (exit_1_array(i) .eq. 0) then
                ppiclf_n_bins(i) = ppiclf_n_bins(i) + 1
                ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                               ppiclf_binb(2*(i-1)+1)  ) / 
     >                              ppiclf_n_bins(i)

                ! Check conditions
                ! exit_1
                total_bin = 1
                do j=1,ppiclf_ndim
                   total_bin = total_bin*ppiclf_n_bins(j)
                enddo
                if (total_bin .gt. ppiclf_np) then
                   ! two exit arrays aren't necessary for now, but
                   ! to make sure exit_2 doesn't slip through, we
                   ! set both for now
                   exit_1_array(i) = 1
                   exit_2_array(i) = 1
                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
     >                                 ppiclf_n_bins(i)
                   exit
                endif
                
                ! exit_2
                if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
                   ! two exit arrays aren't necessary for now, but
                   ! to make sure exit_2 doesn't slip through, we
                   ! set both for now
                   exit_1_array(i) = 1
                   exit_2_array(i) = 1
                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
     >                                 ppiclf_n_bins(i)
                   exit
                endif
             endif
          enddo

          ! full exit_1
          sum_value = 0
          do i=1,ppiclf_ndim
             sum_value = sum_value + exit_1_array(i)
          enddo
          if (sum_value .eq. ppiclf_ndim) then
             exit_1 = .true.
          endif

          ! full exit_2
          sum_value = 0
          do i=1,ppiclf_ndim
             sum_value = sum_value + exit_2_array(i)
          enddo
          if (sum_value .eq. ppiclf_ndim) then
             exit_2 = .true.
          endif
       enddo

! -------------------------------------------------------
! SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! Check for too small bins 
      rthresh = 1E-12
      total_bin = 1
      do i=1,ppiclf_ndim
         total_bin = total_bin*ppiclf_n_bins(i)
         if (ppiclf_bins_dx(i) .lt. rthresh) ppiclf_bins_dx(i) = 1.0
      enddo

!     current box coordinates
      if (ppiclf_nid .le. total_bin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ppiclf_binx(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
         ppiclf_binx(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
         ppiclf_biny(1,1) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
         ppiclf_biny(2,1) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
         ppiclf_binz(1,1) = 0.0d0
         ppiclf_binz(2,1) = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            ppiclf_binz(1,1) = ppiclf_binb(5)+kdum    *ppiclf_bins_dx(3)
            ppiclf_binz(2,1) = ppiclf_binb(5)+(kdum+1)*ppiclf_bins_dx(3)
         endif
      endif

#if PPICLF_ZOLTAN==1
      numGlobObjs = npt_total
      numLocObjs = ppiclf_npart
      grid_dx = ppiclf_d2chk(1)

      call izero(gids,PPICLF_LPART)
      call izero(iprocp,ppiclf_np)
      call izero(iwork,ppiclf_np)
      do i=ppiclf_nid+1,ppiclf_np
        iprocp(i)=ppiclf_npart
      enddo
      call ppiclf_igop(iprocp,iwork,'+  ',ppiclf_np)
      
      istart = 0
      if(ppiclf_nid.gt.0) istart = iprocp(ppiclf_nid)
      do i=1,ppiclf_npart
        gids(i)=istart + i
        part_grid(ix,i)= ppiclf_y(ix,i)
        part_grid(iy,i)= ppiclf_y(iy,i)
        part_grid(iz,i)= ppiclf_y(iz,i)
!! MAP ONTO UNIFORM GRID
c        part_grid(ix,i)= floor((ppiclf_y(ix,i)-ppiclf_binb(1))
c     >                      /grid_dx)*grid_dx +
c     >                       ppiclf_binb(1)
c        part_grid(iy,i)= floor((ppiclf_y(iy,i)-ppiclf_binb(3))
c     >                      /grid_dx)*grid_dx +
c     >                       ppiclf_binb(3)
c        ! protect against 2D  case
c        part_grid(iz,i)= floor((ppiclf_y(iz,i)-ppiclf_binb(iz*2-1))
c     >                      /grid_dx)*grid_dx +
c     >                       ppiclf_binb(iz*2-1)
      enddo

      call partitionWithRCB(ppiclf_comm)
#endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin, idum, jdum, kdum, ndumx, ndumy, itmp, jtmp, ktmp,
     >          i, j, k
!

      nbin = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

c     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_bins_dx(1)/ppiclf_filter) + 1 + 1
         ppiclf_by = floor(ppiclf_bins_dx(2)/ppiclf_filter) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_bins_dx(3)/ppiclf_filter) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_ngrids
         ppiclf_by = ppiclf_by*ppiclf_ngrids
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_ngrids

         if (ppiclf_bx .gt. PPICLF_BX1)
     >      call ppiclf_exittr('Increase PPICLF_BX1$',0.,ppiclf_bx)
         if (ppiclf_by .gt. PPICLF_BY1)
     >      call ppiclf_exittr('Increase PPICLF_BY1$',0.,ppiclf_by)
         if (ppiclf_bz .gt. PPICLF_BZ1)
     >      call ppiclf_exittr('Increase PPICLF_BZ1$',0.,ppiclf_bz)

         ppiclf_rdx = ppiclf_bins_dx(1)/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_bins_dx(2)/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_rdz = ppiclf_bins_dx(3)/(ppiclf_bz-1)

         ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
         ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = sngl(ppiclf_binx(1,1) +
     >                                  (i-1)*ppiclf_rdx)
            ppiclf_grid_y(i,j,k) = sngl(ppiclf_biny(1,1) +
     >                                  (j-1)*ppiclf_rdy)
            ppiclf_grid_z(i,j,k) = sngl(ppiclf_binz(1,1) +
     >                                  (k-1)*ppiclf_rdz)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         enddo
         enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MapOverlapMesh
#if PPICLF_ZOLTAN==1
      use zoltanRCB
#endif
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Internal:
!
      integer*4 icalld
      save      icalld
      data      icalld /0/
      integer*4 nkey(2), i, j, k, ie, iee, ii, jj, kk, ndum, nrank,
     >          nl, nii, njj, nrr, ilow, jlow, klow, nxyz, il,
     >          ihigh, jhigh, khigh, ierr2
      real*8 rxval, ryval, rzval, rthresh
      logical partl
      real*8 ppiclf_vlmin, ppiclf_vlmax
      external ppiclf_vlmin, ppiclf_vlmax
      integer*4 ppiclf_ivlmin, ppiclf_ivlmax
      external ppiclf_ivlmin, ppiclf_ivlmax
#if PPICLF_ZOLTAN==1
      real(Zoltan_DOUBLE) coords(3)
#endif

      ! see which bins are in which elements
      ppiclf_neltb = 0

      rthresh = 1.e-12
c      if(abs(ppiclf_binb(2)-ppiclf_binb(1)).lt.rthresh
c     >   .and. abs(ppiclf_binb(4)-ppiclf_binb(3)).lt.rthresh
c     >   .and. abs(ppiclf_binb(6)-ppiclf_binb(5)).lt.rthresh )
c     >   return
      do ie=1,ppiclf_nee
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         rxval = ppiclf_xm1bs(i,j,k,1,ie)
         ryval = ppiclf_xm1bs(i,j,k,2,ie)
         rzval = 0.0d0
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1bs(i,j,k,3,ie)

         ! Check bounds of particle domain
         if (rxval .gt. ppiclf_binb(2)) goto 1233
         if (rxval .lt. ppiclf_binb(1)) goto 1233
         if (ryval .gt. ppiclf_binb(4)) goto 1233
         if (ryval .lt. ppiclf_binb(3)) goto 1233
         if (ppiclf_ndim.gt.2 .and. rzval .gt. ppiclf_binb(6)) 
     >      goto 1233
         if (ppiclf_ndim.gt.2 .and. rzval .lt. ppiclf_binb(5))
     >      goto 1233

         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim.lt.3) kk = 0
          if (ii .eq. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
          if (jj .eq. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
          if (kk .eq. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
          if (ii .eq. -1) ii = 0
          if (jj .eq. -1) jj = 0
          if (kk .eq. -1) kk = 0
         ! bin number
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
#if PPICLF_ZOLTAN==1
         coords(1) = rxval
         coords(2) = ryval
         coords(3) = rzval
         if(.not.associated(zz_obj)) goto 1233
         ierr = Zoltan_LB_Point_PP_Assign(zz_obj,coords,nrank,ndum)
#else
         nrank = ndum
#endif
         if (ii .lt. 0 .or. ii .gt. ppiclf_n_bins(1)-1) goto 1233
         if (jj .lt. 0 .or. jj .gt. ppiclf_n_bins(2)-1) goto 1233
         if (kk .lt. 0 .or. kk .gt. ppiclf_n_bins(3)-1) goto 1233

         ppiclf_neltb = ppiclf_neltb + 1
         if(ppiclf_neltb .gt. PPICLF_LEE) then
           call ppiclf_exittr('Increase PPICLF_LEE$',0.0d0,ppiclf_neltb)
         endif

         ppiclf_er_map(1,ppiclf_neltb) = ie
         ppiclf_er_map(2,ppiclf_neltb) = ppiclf_nid
         ppiclf_er_map(3,ppiclf_neltb) = ndum
         ppiclf_er_map(4,ppiclf_neltb) = nrank
         ppiclf_er_map(5,ppiclf_neltb) = nrank
         ppiclf_er_map(6,ppiclf_neltb) = nrank

         if (ppiclf_neltb .gt. 1) then
         do il=1,ppiclf_neltb-1
            if (ppiclf_er_map(1,il) .eq. ie) then
            if (ppiclf_er_map(4,il) .eq. nrank) then
               ppiclf_neltb = ppiclf_neltb - 1
               goto 1233
            endif
            endif
         enddo
         endif
 1233 continue
      enddo
      enddo
      enddo
      enddo

      ! copy send data to buffer
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
       iee = ppiclf_er_map(1,ie)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,1,ie)
     >                 ,ppiclf_xm1bs(1,1,1,1,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,2,ie)
     >                 ,ppiclf_xm1bs(1,1,1,2,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,3,ie)
     >                 ,ppiclf_xm1bs(1,1,1,3,iee),nxyz)
      enddo

      ppiclf_neltbb = ppiclf_neltb
      do ie=1,ppiclf_neltbb
         call ppiclf_icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie)
     >             ,PPICLF_LRMAX)
      enddo


      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*3
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltb
     >       ,PPICLF_LEE,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltb
     >       ,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,2)


      ! Map received element nodes to bins
      do ie=1,ppiclf_neltb
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         rxval = ppiclf_xm1b(i,j,k,1,ie)
         ryval = ppiclf_xm1b(i,j,k,2,ie)
         rzval = 0.0d0
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1b(i,j,k,3,ie)
         
         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim.eq.2) kk = 0
          if (ii .eq. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
          if (jj .eq. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
          if (kk .eq. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
          if (ii .eq. -1) ii = 0
          if (jj .eq. -1) jj = 0
          if (kk .eq. -1) kk = 0
          ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                 ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk

#if PPICLF_ZOLTAN==1
         coords(1) = rxval
         coords(2) = ryval
         coords(3) = rzval
         ierr = Zoltan_LB_Point_PP_Assign(zz_obj,coords,nrank,ndum)
#endif
         ppiclf_modgp(i,j,k,ie,1) = ii
         ppiclf_modgp(i,j,k,ie,2) = jj
         ppiclf_modgp(i,j,k,ie,3) = kk
         ppiclf_modgp(i,j,k,ie,4) = ndum ! bin/Zoltan partition
   
      enddo
      enddo
      enddo
      enddo

c      if(ppiclf_neltb.gt.0) 
c     >  write(6,*) 'nid, neltb', ppiclf_nid, ppiclf_neltb
      do ie=1,ppiclf_neltb
         ppiclf_xerange(1,1,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(2,1,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(1,2,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(2,2,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(1,3,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         ppiclf_xerange(2,3,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)

         ilow  = 
     >     floor((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         ihigh = 
     >     floor((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         jlow  = 
     >     floor((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         jhigh = 
     >     floor((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         klow  = 
     >     floor((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         khigh = 
     >     floor((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         if (ppiclf_ndim.lt.3) then
            klow = 0
            khigh = 0
         endif

#if PPICLF_ZOLTAN==1
         ierr = Zoltan_LB_Box_PP_Assign(zz_obj,
     >             ppiclf_xerange(1,1,ie),
     >             ppiclf_xerange(1,2,ie),
     >             ppiclf_xerange(1,3,ie),
     >             ppiclf_xerange(2,1,ie),
     >             ppiclf_xerange(2,2,ie),
     >             ppiclf_xerange(2,3,ie),
     >             nbprocs, numnbprocs, nbparts, numnbparts)
         ppiclf_el_map(1,ie) = ppiclf_ivlmin(nbprocs,numnbprocs)
         ppiclf_el_map(2,ie) = ppiclf_ivlmax(nbprocs,numnbprocs)
#else
         ppiclf_el_map(1,ie) = ilow  + ppiclf_n_bins(1)*jlow  
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*klow
         ppiclf_el_map(2,ie) = ihigh + ppiclf_n_bins(1)*jhigh 
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*khigh
#endif
         ppiclf_el_map(3,ie) = ilow
         ppiclf_el_map(4,ie) = ihigh
         ppiclf_el_map(5,ie) = jlow
         ppiclf_el_map(6,ie) = jhigh
         ppiclf_el_map(7,ie) = klow
         ppiclf_el_map(8,ie) = khigh
      enddo

      if (icalld .eq. 0) then 

         icalld = icalld + 1

         call ppiclf_prints('   *Begin mpi_comm_split$')
            call mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr2)
         call ppiclf_prints('    End mpi_comm_split$')

         call ppiclf_prints('   *Begin InitSolve$')
            call ppiclf_solve_InitSolve
         call ppiclf_prints('    End InitSolve$')

         call ppiclf_io_OutputDiagGrid
      endif

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
     > bind(C, name="ppiclc_comm_InitOverlapMesh")
#else
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8    xgrid(*)
      real*8    ygrid(*)
      real*8    zgrid(*)
!
! External:
!
      integer*4 nxyz, i, j, ie
!
      ppiclf_overlap = .true.

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      if (ncell .gt. PPICLF_LEE .or. ncell .lt. 0) 
     >   call ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      if (lx1 .ne. PPICLF_LEX) 
     >   call ppiclf_exittr('LX1 != LEX in InitOverlap$',0.0d0,ncell)
      if (ly1 .ne. PPICLF_LEY)
     >   call ppiclf_exittr('LY1 != LEY in InitOverlap$',0.0d0,ncell)
      if (lz1 .ne. PPICLF_LEZ)
     >   call ppiclf_exittr('LZ1 != LEZ in InitOverlap$',0.0d0,ncell)

      ppiclf_nee = ncell
      nxyz       = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      call ppiclf_prints("-Begin InitOverlapMesh$")
      do ie=1,ppiclf_nee
      do i=1,nxyz
         j = (ie-1)*nxyz + i
         ppiclf_xm1bs(i,1,1,1,ie) = xgrid(j)
         ppiclf_xm1bs(i,1,1,2,ie) = ygrid(j)
         ppiclf_xm1bs(i,1,1,3,ie) = zgrid(j)
      enddo
      enddo

      call ppiclf_comm_MapOverlapMesh
      
      call ppiclf_prints("-End InitOverlapMesh$")

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_FindParticle
#if PPICLF_ZOLTAN==1
      use zoltanRCB
#endif 
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 ix, iy, iz, i, ii, jj, kk, ndum, nrank
#if PPICLF_ZOLTAN==1
      real(Zoltan_DOUBLE) coords(3)
      integer(Zoltan_INT) nzrank, nzdum
#endif 
!
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq.3)
     >iz = 3

      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ii  = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj  = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim .lt. 3) kk = 0
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
#if PPICLF_ZOLTAN==1

         coords(1) = ppiclf_y(ix,i)
         coords(2) = ppiclf_y(iy,i)
         coords(3) = ppiclf_y(iz,i)
c         if(iz .eq. 1) coords(3) = 0.0
         ! rank from zoltan
         ierr = Zoltan_LB_Point_PP_Assign(zz_obj, coords, nzrank, nzdum)
         nrank = nzrank
         ndum = nzdum
         write(6,*) 'PART',ppiclf_nid, i, nrank, ndum
#else
         ! rank from binning
         nrank = ndum
#endif
         ppiclf_iprop(8 ,i)  = ii
         ppiclf_iprop(9 ,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(3,i)  = nrank ! where particle is actually moved
         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
#if PPICLF_ZOLTAN==1
         ppiclf_iprop(5,i)  = nrank ! CHANGED TAG
         ppiclf_iprop(6,i)  = gids(i) ! CHANGED TAG
         ppiclf_iprop(7,i)  = ppiclf_nid ! CHANGED TAG
#endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveParticle
!
#if PPICLF_ZOLTAN==1
      use zoltanRCB
#endif
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl    
      integer*4 lrf
      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2)
      real*8 rwork(lrf,PPICLF_LPART)
      integer*4 i, ic, j0
!

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
      enddo

      j0 = 4
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart,PPICLF_LPART
     >                                  ,ppiclf_iprop,PPICLF_LIP
     >                                  ,partl,0
     >                                  ,rwork,lrf
     >                                  ,j0)

      if (ppiclf_npart .gt. PPICLF_LPART .or. ppiclf_npart .lt. 0)
     >   call ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)

      write(6,*) 'NumObjs', myrank, ppiclf_npart
      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
      enddo
        
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateGhost
!
#if PPICLF_ZOLTAN==1
      use zoltan
      use zoltanRCB
#endif
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0d0

      do ip=1,ppiclf_npart

         ! preparing particles for projection
         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         ! particle coordinates
         rxval = ppiclf_cp_map(1,ip)
         ryval = ppiclf_cp_map(2,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Particle proc based on binning !!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         iip    = ppiclf_iprop(8,ip)
         jjp    = ppiclf_iprop(9,ip)
         kkp    = ppiclf_iprop(10,ip)

#if PPICLF_ZOLTAN==1
         !! Evaluate bounding box for particle
         distchk = (rfac*grid_dx)
         if(  abs(rxval-locMin(1)).lt.distchk 
     >   .or. abs(rxval-locMax(1)).lt.distchk
     >   .or. abs(ryval-locMin(2)).lt.distchk
     >   .or. abs(ryval-locMax(2)).lt.distchk
     >   .or. abs(rzval-locMin(3)).lt.distchk
     >   .or. abs(rzval-locMax(3)).lt.distchk) then
         !! if near partition boundaries,
         !! check which partitions are within box
!! This would be better if we could check with a sphere
!! instead of a box, but this is good enough for now.
           ierr =  Zoltan_LB_Box_PP_Assign(zz_obj
     >            ,rxval-distchk
     >            ,ryval-distchk
     >            ,rzval-distchk
     >            ,rxval+distchk
     >            ,ryval+distchk
     >            ,rzval+distchk
     >            ,nbprocs
     >            ,numnbprocs
     >            ,nbparts
     >            ,numnbparts)
           do j=1,numnbprocs

             if(nbprocs(j).eq.ppiclf_nid) cycle

             nrank = nbprocs(j)
             ndumn = nrank
             rxnew(1) = rxval
             rxnew(2) = ryval
             rxnew(3) = rzval
       
             ppiclf_npart_gp = ppiclf_npart_gp + 1
             ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
             ppiclf_iprop_gp(2,ppiclf_npart_gp) = iip
             ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjp
             ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkp
             ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn ! destination

             ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
             ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
             ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
             do k=4,PPICLF_LRP_GP
                ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
             enddo

           enddo
         endif
#else
         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip
         rxr = rxl + ppiclf_bins_dx(1)
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
         ryr = ryl + ppiclf_bins_dx(2)
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
            rzr = rzl + ppiclf_bins_dx(3)
         endif

         isave = 0

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Determine where ghost particles will go based on !!
         !! faces, edges and corners of cartesian bins       !!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle !not periodic == 1
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo
#endif

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rxdrng(3)
      integer*4 iadd(3)
!
! Input/Output:
!
      real*8 rxnew(3)
!
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .ge. ppiclf_n_bins(1)) then
         rxnew(1) = rxnew(1) - rxdrng(1)
         goto 123
      endif
      endif
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .lt. 0) then
         rxnew(1) = rxnew(1) + rxdrng(1)
         goto 123
      endif
      endif

  123 continue    
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .ge. ppiclf_n_bins(2)) then
         rxnew(2) = rxnew(2) - rxdrng(2)
         goto 124
      endif
      endif
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .lt. 0) then
         rxnew(2) = rxnew(2) + rxdrng(2)
         goto 124
      endif
      endif
  124 continue

      if (ppiclf_ndim .gt. 2) then
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .ge. ppiclf_n_bins(3)) then
            rxnew(3) = rxnew(3) - rxdrng(3)
            goto 125
         endif
         endif
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .lt. 0) then
            rxnew(3) = rxnew(3) + rxdrng(3)
            goto 125
         endif
         endif
      endif
  125 continue

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl         
!
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
