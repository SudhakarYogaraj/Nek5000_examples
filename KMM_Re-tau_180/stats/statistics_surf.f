!====================================================================== 
!    Name           : statistics_surf.f
!    Author         : Prabal S. Negi
!    Last Modified  : June 10, 2016 
!    Description    : Statistics for defined surfaces 
!    Notes          : Requires the main statistics toolbox developed by
!                   : Prabal S. Negi and Adam Peplinski
!======================================================================  
!---------------------------------------------------------------------- 
!     read parameters Surface statistics 
      subroutine SURF_PARAM_IN(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'SURF_STATS'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SURF/ NSURFS,SURF_DEF,SURF_COMP,SURF_AVGITER 

!     default values
      NSURFS        = 1             ! No. of surfaces
      SURF_DEF      = 'W  '         ! Surface definition
      SURF_COMP     = 5             ! compute interval 
      SURF_AVGITER  = 1             ! saving interval
     
     
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=SURF,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading SURF STAT parameters.$')

!     broadcast data
      call bcast(NSURFS,          ISIZE)
      call bcast(SURF_DEF,        NSURFS*3*CSIZE)
      call bcast(SURF_COMP,       ISIZE)
      call bcast(SURF_AVGITER,    ISIZE)  

      return
      end
!-----------------------------------------------------------------------
!     write parameters relaxation term filtering 
      subroutine SURF_PARAM_OUT(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'SURF_STATS'

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SURF/ NSURFS,SURF_DEF,SURF_COMP,SURF_AVGITER 

!     read the file
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=SURF,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing SURF STAT parameters.$')

      return
      end

!---------------------------------------------------------------------- 

      subroutine SURF_STATS_AVG

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SURF_STATS'
      include 'STATS'

      integer COMPS
      save COMPS
      data COMPS /0/

      integer icalld
      save icalld
      data icalld /0/

      real dtime,talpha,tbeta

      SURF_IDIR = STAT_IDIR

      if (icalld.eq.0) then
        icalld = icalld+1
        call SURF_STAT_INIT
        return
      endif

      if (mod(ISTEP,SURF_COMP).LT.SURF_AVGITER) then

          if (COMPS.eq.0) then
              SURF_TSTART=TIME-DT
              dtime=TIME-SURF_TSTART
              SURF_ATIME=TIME-SURF_TSTART
              tbeta=dtime/SURF_ATIME
              talpha=1.0 - tbeta
          else
              dtime=TIME-SURF_ATIME-SURF_TSTART
              SURF_ATIME=TIME-SURF_TSTART
              tbeta=dtime/SURF_ATIME
              talpha=1.0 - tbeta
          endif
 
          COMPS=COMPS+1
          call SURF_COMPUTE(talpha,tbeta,comps)

          if (COMPS.eq.SURF_AVGITER) then
              call SURF_STATS_SAVE
              COMPS=0
          endif
      endif

      if (SURF_NSAVES.eq.SURF_MAXTSAVES) then
          call SURF_COLLATE
          call SURF_STATS_OUT
!          call exitt                ! prabal  
          continue
      endif

      return
      end subroutine SURF_STATS_AVG

!---------------------------------------------------------------------- 
!    Initialization
      subroutine SURF_STAT_INIT

      implicit none 

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'

      integer lt

      call SURF_DEFINE

!    Assumes mapping for the main statistics has already been done
      call SURF_MAPPING

      call SURF_MASS1D

      call SURF_COMM_MAP

!     Initialize compute arrays
      lt=LX1*LELT*SURF_LVAR*MAXSURFS
      call rzero(SURF_RAVG,lt)

      lt=3*SURF_LVAR*MAXSURFS
      call rzero(SURF_INTV,lt)


!     prabal      
!      call SURF_DEBUG_COMM
!      call exitt        

      return
      end subroutine SURF_STAT_INIT
!----------------------------------------------------------------------

      subroutine SURF_DEFINE

      implicit none 

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'SURF_STATS'

      character cb*3

      integer iel,iface,nfaces,ifld,isf

      ifld=1
      nfaces=2*ndim

      call izero(SURF_MEMCOUNT,NSURFS)
      call izero(SURF_OBJ,2*LELT*6*MAXSURFS)
      if (nid.eq.0) write(6,*) 'Surface Stats: Defining Surface(s)'

      do iel=1,nelt
        do iface = 1,nfaces
          cb = CBC(iface,iel,ifld)
          do isf = 1,nsurfs
            if (cb.EQ.SURF_DEF(isf)) then
               SURF_MEMCOUNT(isf)=SURF_MEMCOUNT(isf)+1 
               SURF_OBJ(1,SURF_MEMCOUNT(isf),isf) = iel
               SURF_OBJ(2,SURF_MEMCOUNT(isf),isf) = iface
            endif
          enddo
        enddo
      enddo

      return
      end subroutine SURF_DEFINE

!----------------------------------------------------------------------  

      subroutine SURF_MAPPING

!    Assumes the main statistics module has been initialized.

      implicit none
     
      include 'SIZE_DEF'
      include 'SIZE'
      include 'STATS'
      include 'SURF_STATS'

      integer isf,sel,sel2,iel,st_el
      integer cntr
      integer st_lmap,st_gmap
      integer tsum

      integer iglsum     ! function: global sum

      integer tmp_lmap(lelt),tmp_gmap(lelt)
    
!    Build local mapping for 2D=>1D
      call izero(SURF_LMAP,LELT*MAXSURFS)
      call ifill(SURF_GMAP,-1,LELT*MAXSURFS)
      call ifill(SURF_OWN,-1,LELT*MAXSURFS)

      call izero(SURF_UNQ,NSURFS)
      call izero(SURF_TUNQ,NSURFS)

      do isf=1,NSURFS
           call izero(tmp_lmap,nelt)

           do sel=1,SURF_MEMCOUNT(isf)
               iel = SURF_OBJ(1,sel,isf)               ! Element no
               st_lmap=STAT_LMAP(iel)                  ! Local mapping from STATS
               SURF_STAT_LMAP(sel,isf)=st_lmap         ! Save this for surface stats
               tmp_lmap(sel) = st_lmap
               tmp_gmap(sel) = STAT_GMAP(st_lmap)
           enddo

           cntr=0
           do sel=1,SURF_MEMCOUNT(isf)
               if (tmp_lmap(sel).ne.-2) then           ! -2 ==> has already been mapped.
                 
!                  if (.neq.-2) ==> Do mapping
!                  Local mapping 2D=>1D

                   cntr=cntr+1
                   st_lmap=SURF_STAT_LMAP(sel,isf)
                   st_gmap=tmp_gmap(sel)
                   do sel2=sel,SURF_MEMCOUNT(isf)
                      if (tmp_gmap(sel2).eq.st_gmap) then
                         SURF_LMAP(sel2,isf)=cntr           ! New mapping for surface stats
                         tmp_lmap(sel2) = -2
                      endif
                   enddo
!                  Global mapping
!                  and Ownership
!                   do st_el=1,STAT_LNUM
!                        if (st_lmap.eq.STAT_GMAP(st_el)) then
!                             SURF_GMAP(cntr,isf)=STAT_GMAP(st_el)
!                             SURF_OWN(cntr,isf)=STAT_OWN(st_el)
!                             exit
!                        endif
!                   enddo
                   SURF_GMAP(cntr,isf)=STAT_GMAP(st_lmap)
                   SURF_OWN(cntr,isf)=STAT_OWN(st_lmap)
                   if (SURF_OWN(cntr,isf).eq.nid) then
                      SURF_UNQ(isf) = SURF_UNQ(isf)+1
                   endif     

                endif
          enddo               ! sel

          SURF_LNUM(isf)=cntr          ! No. of unique surface statistics elements
      enddo                   ! isf

      do isf=1,NSURFS
            tsum = SURF_UNQ(isf)
            SURF_TUNQ(isf) = iglsum(tsum,1)
      enddo
!      write(6,*) nid, 'UNQ: ', SURF_UNQ,SURF_TUNQ,SURF_LNUM     ! prabal
 

      if (nid.eq.0) write(6,*) 'Surface Stats: Element mapping done'

      return
      end subroutine SURF_MAPPING

!---------------------------------------------------------------------- 

      subroutine SURF_MASS1D

      implicit none 

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SURF_STATS'
      include 'MASS_DEF'
      include 'MASS'
      include 'DXYZ_DEF'
      include 'DXYZ'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'WZ_DEF'
      include 'WZ'


      character cb*3

      integer iel,iface,isf
      integer ix,iy,iz
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer cntr

      integer i, j, k, e
      integer el,ie,selno
      integer elpos

      real lwm1(LX1)
      real tmp_sv(LX1,LELT)

!     scratch space
      real lxyzd(LX1,LY1,LZ1,LELT,3)
      common /SCRSF/ lxyzd      ! coordinate derivatives

      real tmp_XM1(LX1,LELT)
      real tmp_YM1(LX1,LELT)

!     copy wieghts depending on the uniform direction
      if (SURF_IDIR.eq.1) then
         call copy(lwm1,WXM1,NX1)
!     get coordinates derivatives d[XYZ]/dr
         i = NY1*NZ1
         do e = 1, NELT
              call mxm(DXM1,NX1,XM1(1,1,1,e),NX1,lxyzd(1,1,1,e,1),i)
              call mxm(DXM1,NX1,YM1(1,1,1,e),NX1,lxyzd(1,1,1,e,2),i)
              call mxm(DXM1,NX1,ZM1(1,1,1,e),NX1,lxyzd(1,1,1,e,3),i)
         enddo
      elseif (SURF_IDIR.eq.2) then
         call copy(lwm1,WYM1,NY1)
!     get coordinates derivatives d[XYZ]/ds
         do e = 1, NELT
              do i=1, NZ1
                 call mxm(XM1(1,1,i,e),NX1,DYTM1,NY1,
     $                lxyzd(1,1,i,e,1),NY1)
                 call mxm(YM1(1,1,i,e),NX1,DYTM1,NY1,
     $                lxyzd(1,1,i,e,2),NY1)
                 call mxm(ZM1(1,1,i,e),NX1,DYTM1,NY1,
     $                lxyzd(1,1,i,e,3),NY1)
              enddo
         enddo
      else
         if (if3d) then            ! #2D 
              call copy(lwm1,WZM1,NZ1)
         else
              call rone(lwm1,NX1)
         endif
!     get coordinates derivatives d[XYZ]/dt
         i = NX1*NY1
         do e = 1, NELT
              call mxm(XM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,1),NZ1)
              call mxm(YM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,2),NZ1)
              call mxm(ZM1(1,1,1,e),i,DZTM1,NZ1,lxyzd(1,1,1,e,3),NZ1)
         enddo

      endif

      i = NX1*NY1*NZ1
!     get arc length
      do e = 1, NELT
           call vsq(lxyzd(1,1,1,e,1),i)
           call vsq(lxyzd(1,1,1,e,2),i)
           call vsq(lxyzd(1,1,1,e,3),i)
      
           call add2(lxyzd(1,1,1,e,1),lxyzd(1,1,1,e,2),i)
           call add2(lxyzd(1,1,1,e,1),lxyzd(1,1,1,e,3),i)

           call vsqrt(lxyzd(1,1,1,e,1),i)
      enddo


!     multiply by appropriate wieghts
      if (SURF_IDIR.EQ.1) then
           do e=1, NELT
                do k=1,NZ1
                do j=1,NY1
                do i=1,NX1
                    SURF_BM1D(i,j,k,e) = lwm1(i)*lxyzd(i,j,k,e,1)
                enddo
                enddo
                enddo
           enddo
      elseif (SURF_IDIR.EQ.2) then
           do e=1, NELT
                do k=1,NZ1
                do j=1,NY1
                do i=1,NX1
                    SURF_BM1D(i,j,k,e) = lwm1(j)*lxyzd(i,j,k,e,1)
                enddo
                enddo
                enddo
           enddo
      elseif (SURF_IDIR.EQ.3) then
           if (if3d) then
                do e=1, NELT
                     do k=1,NZ1
                     do j=1,NY1
                     do i=1,NX1
                         SURF_BM1D(i,j,k,e)=lwm1(k)*lxyzd(i,j,k,e,1)
                     enddo
                     enddo
                     enddo
                enddo
           else
                call rone(SURF_BM1D,lx1*ly1*lz1*nelt)
           endif
      endif

!     get total line length
!     sum contributions within different 3D elements to get 
!     elementwise arc length

      i = LX1*LELT
      call rzero(tmp_sv,i)

      selno=0
      do isf=1,NSURFS
          do ie =1,SURF_MEMCOUNT(isf)
          selno=selno+1
          iel = SURF_OBJ(1,ie,isf)
          iface = SURF_OBJ(2,ie,isf)
                   
          CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $                  ,iface)

!         integration within the element.
          cntr=0
          if (SURF_IDIR.eq.3) then
               do ix=kx1,kx2
                  do iy=ky1,ky2
                     cntr=cntr+1
                     tmp_XM1(cntr,selno)=XM1(ix,iy,kz1,iel)
                     tmp_YM1(cntr,selno)=YM1(ix,iy,kz1,iel)
                     do iz=kz1,kz2
                        if (if3d) then 
                              tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)
                        else
                              tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)
                        endif
                     enddo    ! iz
                  enddo       ! iy
               enddo          ! ix
          elseif (SURF_IDIR.eq.2) then
               do ix=kx1,kx2
                  do iz=kz1,kz2
                     cntr=cntr+1
                     tmp_XM1(cntr,selno)=XM1(ix,ky1,iz,iel)
                     tmp_YM1(cntr,selno)=ZM1(ix,ky1,iz,iel)
                     do iy=ky1,ky2
                         tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)
                     enddo    ! iy
                  enddo       ! iz
               enddo          ! ix
          elseif (SURF_IDIR.eq.1) then
               do iz=kz1,kz2
                  do iy=ky1,ky2
                     cntr=cntr+1
                     tmp_XM1(cntr,selno)=ZM1(kx1,iy,iz,iel)
                     tmp_YM1(cntr,selno)=YM1(kx1,iy,iz,iel)
                     do ix=kx1,kx2
                         tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)
                     enddo    ! ix
                  enddo       ! iy
               enddo          ! iz
          endif

          enddo               ! ie
      enddo                   ! isf

!     Collapse Elements
      call rzero(SURF_ABM1,LX1*LELT*MAXSURFS)
      selno=0
      do isf=1,NSURFS
         do e=1,SURF_MEMCOUNT(isf)
            selno=selno+1
            elpos = SURF_LMAP(e,isf)
            do ix=1,lx1
                 SURF_ABM1(ix,elpos,isf)=SURF_ABM1(ix,elpos,isf)+
     $          tmp_sv(ix,selno)
                 SURF_XM1(ix,elpos,isf) = tmp_XM1(ix,selno)
                 SURF_YM1(ix,elpos,isf) = tmp_YM1(ix,selno)
            enddo     ! ix
         enddo         ! e
      enddo         ! isf

      return
      end subroutine SURF_MASS1D

!----------------------------------------------------------------------  

      subroutine SURF_COMM_MAP

      implicit none
     
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'SURF_STATS'

      integer e,e2,ii,arr_size

      integer collsteps,poffset,dstid,recno,tosend
      integer nt,csteps,lsend,nprocs,tmp_max,pos_cnt,ivlmax

      parameter (arr_size=lelt)
      
      integer logi_procpos(arr_size),tmp_procpos(arr_size)
      integer tmp2_procpos(arr_size)

!------------------------------

      call ifill(SURF_PROCID,-1,lelt)
      call izero(SURF_PROCPOS,lelt)

      collsteps = np/arr_size
      if (mod(np,arr_size).gt.0) then
            collsteps=collsteps+1
      endif

      nprocs=arr_size         ! no of processes to check in each cycle
      SURF_SND_CNT=0
      pos_cnt=0
      recno=0
      SURF_MAXREC=0

      do ii=0,collsteps-1

      call izero(tmp_procpos,arr_size)
      call izero(logi_procpos,arr_size)
      call izero(tmp2_procpos,arr_size)

      if (ii.eq.collsteps-1) then
            if (mod(np,arr_size).gt.0) then
               nprocs=mod(np,arr_size)
            endif
      endif

      poffset=ii*arr_size

      do e=1,nprocs
            dstid=poffset+e-1
            
            if (dstid.eq.nid) then
                  tosend=0
            else
                  call get_ifsurf_send(tosend,dstid)
            endif
            
            if (tosend.gt.0) then
                 SURF_SND_CNT=SURF_SND_CNT+1
                 SURF_PROCID(SURF_SND_CNT)=dstid
            endif

            tmp_procpos(e)=tosend
            logi_procpos(e)=tosend

      enddo       ! e=1,nprocs
      
      call ivgl_running_sum(tmp2_procpos,tmp_procpos,nprocs)          ! running sum across all processes
      call icopy(tmp_procpos,tmp2_procpos,nprocs)
     
      call ibcastn(tmp2_procpos,nprocs,np-1)           ! broadcasts the max receives for each process
                                                       ! for: poffset <= nid < pffset+nprocs       

      tmp_max=ivlmax(tmp2_procpos,nprocs)
      if (tmp_max.gt.SURF_MAXREC) then
            SURF_MAXREC=tmp_max           ! update max number of receives
      endif                               ! across all processes.

      if ((nid.ge.poffset).and.(nid.lt.poffset+nprocs)) then
            SURF_RECNO=tmp2_procpos(nid+1-poffset)    ! get how many NID receives
      endif

      do e2=1,nprocs
            tmp2_procpos(e2)=tmp_procpos(e2)*logi_procpos(e2)
            if (logi_procpos(e2).gt.0) then
                  pos_cnt=pos_cnt+1
                  SURF_PROCPOS(pos_cnt)=tmp2_procpos(e2)
            endif
      enddo       ! end e2=1,nprocs

      enddo       ! end ii=1:collsteps-1

      if (nid.eq.0) write(6,*) 'Surface Stats: Communication map done'
!      write(6,*) nid, 'REC No: ', SURF_RECNO,SURF_MAXREC    ! prabal

      return
      end subroutine SURF_COMM_MAP

!---------------------------------------------------------------------- 

      subroutine SURF_COMPUTE(alpha,beta,comps)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'SURF_STATS'
      include 'PARALLEL_DEF'  ! WDSIZE
      include 'PARALLEL'
      include 'GEOM_DEF'
      include 'GEOM'

      integer surf_npos,surf_lnvar
      integer comps

      real alpha, beta

!     work arrays
      real slvel(LX1,LY1,LZ1,LELT,3), slp(LX1,LY1,LZ1,LELT)
      common /SCRMG/ slvel, slp
      real tmpvel(LX1,LY1,LZ1,LELT,3), tmppr(LX1,LY1,LZ1,LELT)
      common /SCRUZ/ tmpvel, tmppr

      real dudx(LX1,LY1,LZ1,LELT,3) ! velocity derivatives; U
      real dvdx(LX1,LY1,LZ1,LELT,3) ! V
      real dwdx(LX1,LY1,LZ1,LELT,3) ! W
      common /SCRNS/ dudx, dvdx
      common /SCRSF/ dwdx

      character(24) outfmt

      integer ix,iy,iz,e
      integer rix,riy,riz,rie          ! indicies of reference point.
      integer refp                     ! proc owning reference point.
      real x0,y0,z0                    ! reference location.
      real cpref                       ! reference pressure.

      real comp_timing             ! timing
      real dnekclock               ! function
      
      integer len

!-------------------- 
      integer lt
      integer nij
      real sij(lx1,ly1,lz1,3*ldim-3,lelt)

      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      real ur(nxyz),us(nxyz),ut(nxyz)
     $   , vr(nxyz),vs(nxyz),vt(nxyz)
     $   , wr(nxyz),ws(nxyz),wt(nxyz)


      comp_timing = dnekclock()

!    Find point to get reference pressure values.
!    Slight offset to avoid being on element edges
      x0=-4.999
      y0=0.0001
      z0=0.0001
      call get_ref_ind(rix,riy,riz,rie,refp,x0,y0,z0)

      lt=nxyz*nelt

!    prabal
!      call gradm1(dudx(1,1,1,1,1),dudx(1,1,1,1,2),dudx(1,1,1,1,3),VX)
!      call gradm1(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),dvdx(1,1,1,1,3),VY)
!      call gradm1(dwdx(1,1,1,1,1),dwdx(1,1,1,1,2),dwdx(1,1,1,1,3),VZ)

      call mappr(tmppr,PR,tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,3))

      if (refp.eq.nid) then
          cpref = -tmppr(rix,riy,riz,rie)        ! needs to be subrtacted
      endif

      len=WDSIZE*1
      call rbcastn(cpref,len,refp)
!
      call cadd(tmppr,cpref,lt)
      call cmult(tmppr,2.,lt)

      surf_lnvar = 0
      surf_npos = 0

!-------------------------------------------------- 
!     XM1 
      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(XM1,surf_npos,alpha,beta)

!    surf_lnvar = 1
!-------------------------------------------------- 
!     XM1 
      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(YM1,surf_npos,alpha,beta)

!    surf_lnvar = 2
!-------------------------------------------------- 
!     Surface Pressure 
      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(tmppr,surf_npos,alpha,beta)

!    surf_lnvar = 3
!-------------------------------------------------- 
!     ! testing 2nd variable 
      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call rone(tmppr,lt)
!      call vsq(tmppr,lt)

      call SURF_COMPUTE_1DAV1(tmppr,surf_npos,alpha,beta)

!    surf_lnvar = 4
!-------------------------------------------------- 
!     Testing viscous stresses

!         s11 = sij(j1,j2,1,1,e)
!         s21 = sij(j1,j2,1,4,e)
!         s31 = sij(j1,j2,1,6,e)
!
!         s12 = sij(j1,j2,1,4,e)
!         s22 = sij(j1,j2,1,2,e)
!         s32 = sij(j1,j2,1,5,e)
!
!         s13 = sij(j1,j2,1,6,e)
!         s23 = sij(j1,j2,1,5,e)
!         s33 = sij(j1,j2,1,3,e)
!
!         dg(1,1) = pm1(j1,j2,1,e)*n1     ! pressure drag
!         dg(2,1) = pm1(j1,j2,1,e)*n2
!         dg(3,1) = pm1(j1,j2,1,e)*n3
!
!         dg(1,2) = -v*(s11*n1 + s12*n2 + s13*n3) ! viscous drag
!         dg(2,2) = -v*(s21*n1 + s22*n2 + s23*n3)
!         dg(3,2) = -v*(s31*n1 + s32*n2 + s33*n3)

      nij = 3
      if (if3d.or.ifaxis) nij=6
      call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

!     x-component of normal vector?
!     sij needs reshuffling.       
      do e=1,nelt
        call copy(slvel(1,1,1,e,1),sij(1,1,1,1,e),nxyz)
        call copy(slvel(1,1,1,e,2),sij(1,1,1,4,e),nxyz)
        call copy(slvel(1,1,1,e,3),sij(1,1,1,6,e),nxyz)
      enddo
  
      call get_normal_comp(tmpvel(1,1,1,1,1),slvel(1,1,1,1,1),1)
      call get_normal_comp(tmpvel(1,1,1,1,2),slvel(1,1,1,1,2),2)
      call get_normal_comp(tmpvel(1,1,1,1,3),slvel(1,1,1,1,3),3)

      call add3(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,3),
     $      lt)
      call col2(tmpvel,vdiff,lt)

      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(tmpvel(1,1,1,1,1),surf_npos,alpha,beta)

!    surf_lnvar = 5
!-------------------------------------------------- 
!     y-component of normal vector?
      do e=1,nelt
        call copy(slvel(1,1,1,e,1),sij(1,1,1,4,e),nxyz)
        call copy(slvel(1,1,1,e,2),sij(1,1,1,2,e),nxyz)
        call copy(slvel(1,1,1,e,3),sij(1,1,1,5,e),nxyz)
      enddo
  
      call get_normal_comp(tmpvel(1,1,1,1,1),slvel(1,1,1,1,1),1)
      call get_normal_comp(tmpvel(1,1,1,1,2),slvel(1,1,1,1,2),2)
      call get_normal_comp(tmpvel(1,1,1,1,3),slvel(1,1,1,1,3),3)

      call add3(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,3),
     $      lt)
      call col2(tmpvel,vdiff,lt)

      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(tmpvel(1,1,1,1,1),surf_npos,alpha,beta)

!    surf_lnvar = 6
!-------------------------------------------------- 
!     z-component of normal vector?
      do e=1,nelt
        call copy(slvel(1,1,1,e,1),sij(1,1,1,6,e),nxyz)
        call copy(slvel(1,1,1,e,2),sij(1,1,1,5,e),nxyz)
        call copy(slvel(1,1,1,e,3),sij(1,1,1,3,e),nxyz)
      enddo
  
      call get_normal_comp(tmpvel(1,1,1,1,1),slvel(1,1,1,1,1),1)
      call get_normal_comp(tmpvel(1,1,1,1,2),slvel(1,1,1,1,2),2)
      call get_normal_comp(tmpvel(1,1,1,1,3),slvel(1,1,1,1,3),3)

      call add3(tmpvel(1,1,1,1,1),tmpvel(1,1,1,1,2),tmpvel(1,1,1,1,3),
     $      lt)
      call col2(tmpvel,vdiff,lt)

      surf_lnvar = surf_lnvar + 1
      surf_npos = surf_lnvar

      call SURF_COMPUTE_1DAV1(tmpvel(1,1,1,1,1),surf_npos,alpha,beta)

!    surf_lnvar = 7
!--------------------------------------------------

      if (surf_lnvar.gt.SURF_LVAR) then
        if (nid.eq.0) then
          write(6,*) 'Inadequate SURF_LVAR. Modify in SURF_STATS
     $and recompile'
        endif    
        call exitt
      endif 

      comp_timing = dnekclock() - comp_timing

 9005 format(A24,2X,I5,2X,A7,2X,E15.8E2)
      if (nid.eq.0) write(6,9005) 'Surface Stats: Compute:',comps,
     $         'Timing:', comp_timing

      return
      end subroutine SURF_COMPUTE 

!---------------------------------------------------------------------- 

      subroutine SURF_COMPUTE_1DAV1(lvar,npos,alpha,beta) 

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'         ! 2D statistics speciffic variables
      include 'INPUT_DEF'
      include 'INPUT'              ! IF3D
      include 'MASS_DEF' 
      include 'MASS'               ! BM1   


!     argument list
      real lvar(LX1,LY1,LZ1,LELT) ! integrated variable
      integer npos              ! Variable position in SURF_AVG_VAR
      real alpha, beta          ! time averaging parameters

!     local variables
      integer i, j, k, e        ! loop index
      integer iel               ! index of 2D element

!     New variables for surface integeration 
      integer isf,iface
      integer cntr,selno        ! counters
      integer avgpts
      integer elpos

      integer ix,iy,iz,kx1,kx2,ky1,ky2,kz1,kz2

      real tmp_sv(lx1,lelt)
      real TMP_AVG_VAR(lx1,lelt)
      real TMP_SINTV(3,MAXSURFS)

!-------------------------------------------------- 

!     zero work array
      i=LX1*LELT
      call rzero(tmp_sv,i)
      i=3*NSURFS
      call rzero(TMP_SINTV,i)

!     perform 1D integral
      selno=0
      do isf=1,NSURFS
      do e = 1,SURF_MEMCOUNT(isf)
          selno=selno+1
          iel = SURF_OBJ(1,e,isf)                       ! Element no.
          iface = SURF_OBJ(2,e,isf)                     ! Face no 
!                   Get appropriate indicies

          call CALC_SURF_INTEGRAL(TMP_SINTV(1,isf),lvar,iface,iel)   ! calculate lvar*area (3 components) 

          CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $                  ,iface)
    
!         integration within the element.
          cntr=0
          if (SURF_IDIR.eq.3) then
               do ix=kx1,kx2
                  do iy=ky1,ky2
                     cntr=cntr+1
                     do iz=kz1,kz2
                        if (if3d) then 
                              tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)*lvar(ix,iy,iz,iel)
                        else
                              tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    lvar(ix,iy,iz,iel)
                        endif
                     enddo    ! iz
                  enddo       ! iy
               enddo          ! ix
          elseif (SURF_IDIR.eq.2) then
               do ix=kx1,kx2
                  do iz=kz1,kz2
                     cntr=cntr+1
                     do iy=ky1,ky2
                         tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)*lvar(ix,iy,iz,iel)
                     enddo    ! iy
                  enddo       ! iz
               enddo          ! ix
          elseif (SURF_IDIR.eq.1) then
               do iz=kz1,kz2
                  do iy=ky1,ky2
                     cntr=cntr+1
                     do ix=kx1,kx2
                         tmp_sv(cntr,selno)=tmp_sv(cntr,selno)+
     $    SURF_BM1D(ix,iy,iz,iel)*lvar(ix,iy,iz,iel)
                     enddo    ! ix
                  enddo       ! iy
               enddo          ! iz
          endif
      enddo         ! e=1,SURF_MEMCOUNT(isf)
      enddo         ! isf=1,NSURFS

!     Collapse Elements
      selno=0
      do isf=1,NSURFS
         call rzero(TMP_AVG_VAR,LX1*SURF_LNUM(isf))
         do e=1,SURF_MEMCOUNT(isf)
             selno=selno+1
             elpos = SURF_LMAP(e,isf)
             call add2(TMP_AVG_VAR(1,elpos),tmp_sv(1,selno),LX1)
         enddo         ! e

!        Time average. 1D integeral
         avgpts=lx1*SURF_LNUM(isf)
         call add2sxy(SURF_RAVG(1,1,npos,isf),alpha,
     $         TMP_AVG_VAR(1,1),beta,avgpts)

!        Time average. Surf integeral
         avgpts=3
         call add2sxy(SURF_INTV(1,npos,isf),alpha,
     $         TMP_SINTV(1,isf),beta,avgpts)

      enddo         ! isf

      return
      end subroutine SURF_COMPUTE_1DAV1

!----------------------------------------------------------------------

      subroutine SURF_STATS_SAVE

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SURF_STATS'

      integer isf,svar
      integer i,lt

      SURF_NSAVES=SURF_NSAVES+1
      SURF_TSTAMP(SURF_NSAVES)=TIME

      if (nid.eq.0) write(6,*) 'Saving Surface history'

      do isf=1,NSURFS
         do svar=1,SURF_LVAR

!           Save 1D integral
             lt=SURF_LNUM(isf)*lx1
             call copy(SURF_AVG_HIST(1,1,SURF_NSAVES,svar,isf),
     $            SURF_RAVG(1,1,svar,isf),lt)
             call rzero(SURF_RAVG(1,1,svar,isf),lt)       ! zero out parts already saved

!           Save surface integral components
            do i=1,ndim
              call copy(SURF_INTV_HIST(SURF_NSAVES,i,svar,isf),
     $           SURF_INTV(i,svar,isf),1)
              SURF_INTV(i,svar,isf) = 0.
            enddo

         enddo
      enddo

      return
      end subroutine SURF_STATS_SAVE

!---------------------------------------------------------------------- 

      subroutine SURF_COLLATE

      implicit none

      include 'mpif.h'
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'SURF_STATS'
      include 'INPUT_DEF'
      include 'INPUT'                   !if3d

      integer e2,cnt2,cnt3,len
      integer isf,it
      integer c_step                        ! communication step no.
      integer recsize
      integer recmap(lelt)
      integer rec_dist(MAXSURFS)

      integer status(mpi_status_size)
      integer msg_id1,msg_id2,msg_id3,msg_id4          ! request handles for irecv
      integer msgtag
      integer ierr
      integer dstid
      integer maxpos
     
      integer ivlmax          ! function     integer vector max
      integer ivlsum          ! function     integer vector sum

      integer snds(lelt)
      integer irecv

      real recdata(LX1,LELT*SURF_LVAR*SURF_MAXTSAVES*MAXSURFS)
      real arclen(LX1,LELT*MAXSURFS)

      integer tot_surf_lnum         ! Total unique surf elements  
      integer tot_surf_rec          ! Total surface elements received (in each step)  

      integer pos

      real wk(3*SURF_LVAR*SURF_MAXTSAVES*MAXSURFS)

!     Collate Surface integral values
      call gop(SURF_INTV_HIST,wk,'+  ',
     $       3*SURF_LVAR*SURF_MAXTSAVES*MAXSURFS)

      maxpos=ivlmax(SURF_PROCPOS,SURF_MAXREC)     ! local comm. steps 
      tot_surf_lnum=ivlsum(SURF_LNUM,NSURFS)

      do c_step=1,SURF_MAXREC
!       Receive distribution            
        if (c_step.le.SURF_RECNO) then              ! How many elements for each surface 
          len=(NSURFS)*ISIZE
          msgtag=0 
          msg_id1 = irecv(msgtag,rec_dist,len)       ! rec surface distribution
        endif

!       Send Distribution    
        cnt2=0                  ! counter for processes to send to in this step
        if (c_step.le.maxpos) then
          do e2=1,SURF_SND_CNT                            ! could send to multiple
                                                          ! processes in one step
            if (SURF_PROCPOS(e2).eq.c_step) then
              dstid=SURF_PROCID(e2)
              cnt2=cnt2+1
              snds(cnt2)=dstid
              
              call SENDDIST(dstid)                ! send surface distribution  
            endif
          enddo
        endif

        if (c_step.le.SURF_RECNO) then
          call msgwait(msg_id1)                     ! wait to receive distribution 
          tot_surf_rec = ivlsum(rec_dist,NSURFS)
        endif

!       Receive mapping
        if (c_step.le.SURF_RECNO) then
          len=(tot_surf_lnum)*isize
          msgtag=1 
          msg_id2 = irecv(msgtag,recmap,len)             ! rec surface mapping 
        endif

!       Send mapping 
        if (c_step.le.maxpos) then
          do e2=1,SURF_SND_CNT                       ! could send to multiple
                                                     ! processes in one step
            if (SURF_PROCPOS(e2).eq.c_step) then
              dstid=SURF_PROCID(e2)
              call SENDSURFMAP(dstid)            ! send mapping 
            endif
          enddo
        endif

!       Receive Data/Arc length 
        if (c_step.le.SURF_RECNO) then
          call msgwait(msg_id2)                ! wait to receive mapping

          call TRNSFM_SURF_RECMAP(recmap,rec_dist) ! transform global to local mapping

          len=WDSIZE*tot_surf_rec*LX1*SURF_LVAR*SURF_MAXTSAVES
          msgtag=2
          msg_id3 = irecv(msgtag,recdata,len)      ! rec data

          len=wdsize*tot_surf_rec*LX1
          msgtag=3
          msg_id4 = irecv(msgtag,arclen,len)       ! rec arclength
        endif
        
!       Send data
        do e2=1,cnt2
          dstid=snds(e2)
          call SENDSURFDATA(dstid)
        enddo

!       Waiting to receive data and arclength
        if (c_step.le.SURF_RECNO) then
          call msgwait(msg_id3)               ! wait to receive data
          call msgwait(msg_id4)               ! wait to receive arclength

          call COLLSURFDATA(recmap,rec_dist,recdata,arclen)
        endif

        call nekgsync            
      enddo             ! end e=1,surf_maxrec

      SURF_INI=.true.

!     Divide data by arc length
      if (if3d) then
        do isf=1,NSURFS
          do e2=1,SURF_LVAR
             do it=1,SURF_MAXTSAVES
               len=LX1*SURF_LNUM(isf)
               call invcol2(SURF_AVG_HIST(1,1,it,e2,isf),
     $             SURF_ABM1(1,1,isf),len)
             enddo
          enddo
        enddo
      endif

      SURF_NSAVES=0
      if (nid.eq.0) write(6,*) 'Surface Stats: Collation done'
      
      return
      end subroutine SURF_COLLATE 

!----------------------------------------------------------------------
 
      subroutine GET_IFSURF_SEND(tosend,dstid)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'

      integer dstid,tosend
      integer isf,sel 

      tosend=0

      do isf=1,NSURFS
        do sel=1,SURF_LNUM(isf)
          if (SURF_OWN(sel,isf).eq.dstid) then
            tosend=1
            return
          endif
        enddo
      enddo

      return
      end subroutine GET_IFSURF_SEND

!----------------------------------------------------------------------

      subroutine SENDDIST(dstid)
!     Number of elements that are being communicated for each surface object
!     Send: (5 8 0) ==> 5 elemts for surface 1, 8 elements for surface 2, 0 elements for surface 3
!     will be sent to nid==dstid

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'SURF_STATS'

      integer dstid
      integer len
      integer isf,e,s_els

      integer msgtag

      integer snddist(MAXSURFS)
      
      len=NSURFS*ISIZE

      do isf=1,NSURFS
        s_els=0                                ! elements to send for surface=isf
        do e=1,SURF_LNUM(isf)
          if (SURF_OWN(e,isf).eq.dstid) s_els=s_els+1
        enddo
        snddist(isf)=s_els
      enddo
      
      msgtag = 0
      call csend(msgtag,snddist,len,dstid,0)

      return
      end subroutine SENDDIST

!---------------------------------------------------------------------- 

      subroutine SENDSURFMAP(dstid)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer dstid
      integer len
      integer isf,e,cntr
      integer msgtag

      integer sndmap(lelt)

      cntr=0
      do isf=1,NSURFS
        do e=1,SURF_LNUM(isf)
          if (SURF_OWN(e,isf).eq.dstid) then
            cntr=cntr+1
            sndmap(cntr)=SURF_GMAP(e,isf)
          endif
        enddo
      enddo
      
      len=cntr*isize
      msgtag=1
      call csend(msgtag,sndmap,len,dstid,0)

      end subroutine SENDSURFMAP

!----------------------------------------------------------------------

      subroutine TRNSFM_SURF_RECMAP(recmap,rec_dist)

      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'

      integer rec_dist(MAXSURFS)
      integer recmap(lelt)

      integer i,e,cntr
      integer isf
      integer gl_no,scount

      
      cntr=0
      do isf=1,NSURFS
        scount=rec_dist(isf)
        do e=1,scount
          cntr=cntr+1
          gl_no=recmap(cntr)
          do i=1,SURF_LNUM(isf)
            if (gl_no.eq.SURF_GMAP(i,isf)) then
              recmap(cntr)=i
              exit
            endif
          enddo
        enddo
      enddo

      return
      end subroutine TRNSFM_SURF_RECMAP

!----------------------------------------------------------------------

      subroutine COLLSURFDATA(recmap,rec_dist,recdata,arclen)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'
      
      integer e,it,svar,isf
      integer len
      integer recmap(LELT)         ! Mapping
      integer pos                  ! Mapped position in array
      integer cntr                 ! Counter
      integer rec_dist(MAXSURFS)   ! No. elements for each surface
      integer scount
      integer ofst                 ! offset
      integer ivlsum               ! function 

      real recdata(LX1,LELT*SURF_MAXTSAVES*SURF_LVAR*MAXSURFS)
      real arclen(LX1,LELT*MAXSURFS)

      cntr=0
      ofst=0
      do isf=1,NSURFS
        scount=rec_dist(isf)
        do svar=1,SURF_LVAR
          do it=1,SURF_MAXTSAVES
            ofst = ivlsum(rec_dist,isf-1)
            do e=1,scount
              pos=recmap(ofst+e)
              len=LX1
              cntr=cntr+1
              call add2(SURF_AVG_HIST(1,pos,it,svar,isf),
     $               recdata(1,cntr),len)
            enddo
          enddo
        enddo
      enddo

!     Sum up arc lengths
      if (SURF_INI.eqv..false.) then
        cntr=0 
        do isf=1,NSURFS
          scount=rec_dist(isf)
          do e=1,scount
            cntr=cntr+1
            pos=recmap(cntr)
            len=LX1
            call add2(SURF_ABM1(1,pos,isf),arclen(1,cntr),len)
          enddo
        enddo
      endif

      return
      end subroutine COLLSURFDATA 

!---------------------------------------------------------------------- 

      subroutine SENDSURFDATA(dstid)
      
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SURF_STATS'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer dstid
      integer len
      integer isf,e,svar,sel,cntr
      integer it

      integer msgtag

      real snddata(LX1,LELT*SURF_MAXTSAVES*SURF_LVAR*MAXSURFS)
      real arclen(LX1,LELT*MAXSURFS)

      cntr=0
      do isf=1,NSURFS
        do svar=1,SURF_LVAR
          do it=1,SURF_MAXTSAVES
            do sel=1,SURF_LNUM(isf)
              if (SURF_OWN(sel,isf).eq.dstid) then
                cntr=cntr+1
                len=LX1
                call copy(snddata(1,cntr),
     $              SURF_AVG_HIST(1,sel,it,svar,isf),len)
              endif 
            enddo
          enddo
        enddo
      enddo
      
      len=cntr*LX1*wdsize
      msgtag=2
      call csend(msgtag,snddata,len,dstid,0)      ! send data

      cntr=0
      do isf=1,NSURFS
        do sel=1,SURF_LNUM(isf)
          if (SURF_OWN(sel,isf).eq.dstid) then
            cntr=cntr+1
            len=LX1
            call copy(arclen(1,cntr),SURF_ABM1(1,sel,isf),len)
          endif
        enddo
      enddo
      
      len=cntr*LX1*wdsize
      msgtag=3
      call csend(msgtag,arclen,len,dstid,0)       ! send arclength

      return
      end subroutine SENDSURFDATA

!---------------------------------------------------------------------- 
      subroutine SURF_STATS_OUT

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'TSTEP_DEF'
      include 'TSTEP'                     ! Time
      include 'SURF_STATS'      
!      include 'RESTART'                   ! nek IO
      include 'mpif.h'

      integer iel,counter,cntr

      real wallxy(LX1,LELT*MAXSURFS,2)
      real wallvar(LX1,LELT*SURF_LVAR*SURF_MAXTSAVES*MAXSURFS)

      integer i,j,ix,it
      integer isf,sel,selno,selno2,svar,scomp

      integer cnt_rem(MAXSURFS)
      integer totpts

      integer msg_id1,msg_id2,msg_id3
      integer msgtag

      integer ivlsum,iglsum,irecv                ! function
      integer funits(MAXSURFS)
      integer tsum

      character(39) outfmt
      character(51) outfmt1
      character(25) outfmt2
      character(16) realfmt1
      character(19) realfmt2
      character(17) realfmt3
      character(30) fname

      character(132) shdr

      integer ierr,length,ip

      logical ifbinary                          ! binary or ASCII output
      character(4) fext

      integer s_wdsizeio,s_isizeio,s_csizeio    ! modifiers for 32/64 bit  
      integer recpos                            ! rec position
      integer reclen                            ! rec length
      integer hdrofst                           ! offset for header
      integer tstampofst                        ! timestamp offset
      integer intgofst                          ! integral vals offset
      integer cnt_rsum(MAXSURFS)                ! running sum of elements

      real sout_timing                          ! timing for I/O
      real dnekclock                            ! function 


      sout_timing = dnekclock()

      s_wdsizeio=1                              ! double precision
      s_isizeio=1
      s_csizeio=1
      ifbinary=.true.

      if (ifbinary) then
        fext='.bin'
      else
        fext='.dat'
      endif

! 1001 format (I5.5,2X,E14.6E3,2X,E14.6E3,2X,E14.6E3)    ! nid,x,y,p
      outfmt2 = "(A5,2X,A14,2X,A14,2X,A14)"    ! nid,x,y,p
      write(realfmt1,"(A1,I2.2,A13)") "(", LX1, "(E14.6E3,2X))"
      write(realfmt2,"(A3,I2.2,A14)") "(2(", LX1, "(E14.6E3,2X)))"
      write(realfmt3,"(A1,I3.3,A13)") "(", SURF_MAXTSAVES,
     $              "(E14.6E3,2X))"

      call blank(shdr,132)
      recpos=1
      reclen=S_WDSIZEIO       ! just initialization
      hdrofst=9
      tstampofst=SURF_MAXTSAVES*1
      
      totpts=0
      if (nid.eq.0) then

         if (ifbinary) then 
           write(6,*) 'Output Surface history: Binary file'
         else 
           write(6,*) 'Output Surface history: Ascii file.'
           write(6,*) 'Ascii output not implemented yet.'
           write(6,*) 'Switching to binary Output.'
           ifbinary=.true.
           fext='.bin'
         endif

!     File I/O
         do i = 1,NSURFS
           funits(i)=73+i
           write(fname,9002) 'surf_data',i,'_',TIME,fext
 9002      format (A9,I2.2,A1,E14.8E2,A4)

           open(unit=funits(i), FILE=trim(fname),form="unformatted",
     $         ACCESS="DIRECT",recl=8)

!          call mbyte_open(trim(fname),fid0,ierr)           ! nek IO. prabal

!     Write hdr
           write(shdr,9003) '#',s_wdsizeio,'Surf No:',i,LX1, 
     $     SURF_TUNQ(i),'XS',SURF_LVAR,SURF_MAXTSAVES,TIME

 9003      format (A1,2X,I1,2X,A8,2X,I3.3,2X,I2.2,2X,I5,2X,A2,I2.2,2X,
     $         I3.3,2X,E14.6E3)

           if (ifbinary) then
!             write(funits(i),rec=recpos) shdr                 ! header
!             write hdr
!             call byte_write(shdr,iHeaderSize/4,ierr)         ! nek IO. prabal
!             call byte_write(test_pattern,1,ierr)

!           prabal
!           nek IO
!             ! if we want to switch the bytes for output
!             ! switch it again because the hdr is in ASCII
!             call get_bytesw_write(ibsw_out)
!             if (ibsw_out.ne.0) call set_bytesw_write(ibsw_out)
!             if (ibsw_out.ne.0) call set_bytesw_write(0)  
!           nek IO

             recpos=1
!             write(funits(i),rec=recpos) '##'
!             recpos=recpos+1
             write(funits(i),rec=recpos) s_wdsizeio+0.
             recpos=recpos+1
             write(funits(i),rec=recpos) i+0.
             recpos=recpos+1
             write(funits(i),rec=recpos) LX1+0.
             recpos=recpos+1
             write(funits(i),rec=recpos) SURF_TUNQ(i)+0.
!             recpos=recpos+1
!             write(funits(i),rec=recpos) 'XS'
             recpos=recpos+1
             write(funits(i),rec=recpos) SURF_LVAR+0.
             recpos=recpos+1
             write(funits(i),rec=recpos) SURF_MAXTSAVES+0.
             recpos=recpos+1
             write(funits(i),rec=recpos) TIME
!            End of header

             hdrofst=recpos 

!            write Time stamps
             do ix=1,SURF_MAXTSAVES
               recpos=recpos+1
               write(funits(i),rec=recpos) SURF_TSTAMP(ix)
             enddo
             tstampofst=SURF_MAXTSAVES 

!            write integral values
             do sel=1,SURF_LVAR
               do scomp=1,NDIM 
                 do ix=1,SURF_MAXTSAVES
                   recpos=recpos+1
                   write(funits(i),rec=recpos) 
     $                   SURF_INTV_HIST(ix,scomp,sel,i)
                 enddo
               enddo
             enddo   
             intgofst=NDIM*SURF_LVAR*SURF_MAXTSAVES 
           endif
         enddo 

         counter=0 
         do isf = 1,NSURFS
           do sel = 1,SURF_LNUM(isf)
              if (SURF_OWN(sel,isf).eq.nid) then   
!     XM1/YM1
                if (ifbinary) then

!                  write(6,*) 'recpos', recpos         ! prabal
                  do ix=1,LX1
                    recpos=hdrofst+tstampofst+intgofst+2*counter + ix
!                    write(6,*) recpos                  ! prabal

!                    call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                    call byte_write(lglel,nelt,ierr)

                    write(funits(isf),rec=recpos) SURF_XM1(ix,sel,isf)
                  enddo
                  do ix=1,LX1 
                    recpos=hdrofst+tstampofst+intgofst+2*counter+LX1+ix
!                    write(6,*) recpos                  ! prabal

!                    nek IO. prabal   
!                    call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                    call byte_write(lglel,nelt,ierr)
!                    nek IO. prabal
   
                    write(funits(isf),rec=recpos) SURF_YM1(ix,sel,isf)
                  enddo 

                endif
                counter=counter+LX1
              endif
           enddo        ! sel
!     data
           do svar=1,SURF_LVAR 
             do it=1,SURF_MAXTSAVES
               cntr=0
               do sel=1,SURF_LNUM(isf)
                 if (SURF_OWN(sel,isf).eq.nid) then
                   cntr=cntr+1   
                   if (ifbinary) then
                     do ix=1,LX1
                       recpos = hdrofst+tstampofst+intgofst
     $                 + ( 2*SURF_TUNQ(isf)*LX1
     $                 + SURF_TUNQ(isf)*(svar-1)*LX1*SURF_MAXTSAVES
     $                 + SURF_TUNQ(isf)*LX1*(it-1)
     $                 + (cntr-1)*LX1 ) + ix

!                      nek IO. prabal   
!                      call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                      call byte_write(lglel,nelt,ierr)
!                      nek IO. prabal

                       write(funits(isf),rec=recpos) SURF_AVG_HIST
     $                  (ix,sel,it,svar,isf)
                     enddo 
                   endif
                 endif  
               enddo    ! sel
             enddo      ! it
           enddo        ! svar 
           cnt_rsum(isf) = counter/LX1
         enddo          ! isf

         totpts=counter
         
         do ip=1,np-1
!     hand sahking
           length = 1*isize
           ierr = ip
           call csend(ip,ierr,length,ip,0)
!     count number
           length=NSURFS*isize
           call crecv(ip,cnt_rem,length)
!     get corrdinates 
           counter=ivlsum(cnt_rem,NSURFS)
           totpts=totpts+counter*LX1
           if (counter.gt.0) then
             length = wdsize*counter*LX1
             msgtag=0
             msg_id1 = irecv(msgtag,wallxy(1,1,1),length)         ! XM1
             msgtag=1
             msg_id2 = irecv(msgtag,wallxy(1,1,2),length)         ! YM1

             length=wdsize*counter*LX1*SURF_LVAR*SURF_MAXTSAVES  
             msgtag=2
             msg_id3 = irecv(msgtag,wallvar(1,1),length)         ! data

             call msgwait(msg_id1)
             call msgwait(msg_id2)
             call msgwait(msg_id3)

!     write XM1/YM1/Data
             selno=0
             selno2=0  
             do isf = 1,NSURFS
               if (ifbinary) then
                 do sel = 1,cnt_rem(isf)
                   selno=selno+1
!     XM1
                   do ix=1,LX1
                     recpos=hdrofst + tstampofst + intgofst
     $            + 2*cnt_rsum(isf)*LX1 + 2*(sel-1)*LX1 + ix

!                      nek IO. prabal   
!                      call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                      call byte_write(lglel,nelt,ierr)
!                      nek IO. prabal


                     write(funits(isf),rec=recpos) wallxy(ix,selno,1)
                   enddo 
!    YM1
                   do ix=1,LX1
                     recpos=hdrofst + tstampofst + intgofst 
     $               + 2*cnt_rsum(isf)*LX1 + 2*(sel-1)*LX1 + LX1 + ix

!                      nek IO. prabal   
!                      call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                      call byte_write(lglel,nelt,ierr)
!                      nek IO. prabal

                     write(funits(isf),rec=recpos) wallxy(ix,selno,2)
                   enddo 

                 enddo 
!     data
                 do svar=1,SURF_LVAR 
                   do it=1,SURF_MAXTSAVES
                     do sel=1,cnt_rem(isf)
                       selno2=selno2+1
                       do ix=1,LX1  
                         recpos = hdrofst + tstampofst + intgofst
     $                   + ( 2*SURF_TUNQ(isf)*LX1
     $                   + SURF_TUNQ(isf)*(svar-1)*LX1*SURF_MAXTSAVES
     $                   + SURF_TUNQ(isf)*LX1*(it-1)
     $                   + cnt_rsum(isf)*LX1 
     $                   + (sel-1)*LX1 ) + ix

!                        nek IO. prabal   
!                        call byte_set_view (ioff,ifh_mbyte)    ! nek IO. prabal
!                        call byte_write(lglel,nelt,ierr)
!                        nek IO. prabal

                         write(funits(isf),rec=recpos) wallvar
     $                                       (ix,selno2)
                       enddo  
                     enddo    ! sel
                   enddo      ! it
                 enddo        ! svar 
               endif          ! ifbinary
               cnt_rsum(isf) = cnt_rsum(isf) + cnt_rem(isf)
             enddo            ! isf
           endif              ! counter.gt.0
         enddo                ! ip 
      else
!    (nid.ne.0)

!     Coordinates 
!     hand shaking
        length = 1*isize
        call crecv(nid,ierr,length)

        length = NSURFS*isize
        call csend(nid,SURF_UNQ,length,0,0)

        counter=ivlsum(SURF_UNQ,NSURFS)
!       send data
        if(counter.gt.0) then
!         XM1/YM1
          cntr=0   
          do isf=1,NSURFS
            do sel=1,SURF_LNUM(isf)
              if (SURF_OWN(sel,isf).eq.nid) then
                cntr=cntr+1   
                length=LX1
!               XM1   
                call copy(wallxy(1,cntr,1),
     $              SURF_XM1(1,sel,isf),length)
!               YM1
                call copy(wallxy(1,cntr,2),
     $              SURF_YM1(1,sel,isf),length)
              endif
            enddo
          enddo
          length = wdsize*cntr*LX1
          msgtag=0
          call csend(msgtag,wallxy(1,1,1),length,0,0)
          msgtag=1
          call csend(msgtag,wallxy(1,1,2),length,0,0)

!!        Data
          cntr=0 
          do isf=1,NSURFS
            do svar=1,SURF_LVAR
              do it=1,SURF_MAXTSAVES   
!             copy data
                do sel=1,SURF_LNUM(isf)
                  if (SURF_OWN(sel,isf).eq.nid) then
                      cntr=cntr+1 
                      length=LX1
                      call copy(wallvar(1,cntr),
     $          SURF_AVG_HIST(1,sel,it,svar,isf),length)
                  endif
                enddo             ! sel   
              enddo               ! it
            enddo                 ! svar
          enddo                   ! isf  
!
!         send data
          length = WDSIZE*cntr*LX1
          msgtag=2
          call csend(msgtag,wallvar(1,1),length,0,0)
        endif

      endif         ! nid.eq.0

!     End of writing 
!---------------------------------------- 
!     Close all files
      do i=1,NSURFS
         close(funits(i))
      enddo
!
      call nekgsync()

      sout_timing = dnekclock() - sout_timing

      if(nid.eq.0)then
         write(6,9004)  'Total PTS:', totpts, 'I/O Time:', sout_timing
      end if

 9004 format (A10,2X,I5,2X,A9,2X,E14.6E2)

      return
      end subroutine SURF_STATS_OUT

!---------------------------------------------------------------------- 
!
      subroutine SURF_DEBUG_COMM

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'SURF_STATS'      
      include 'STATS'
      include 'mpif.h'

      integer iel,counter

      real wall(lx1,lelt*MAXSURFS),wall_tmp(lx1,lelt*MAXSURFS)
      integer iwall(lx1*lelt*MAXSURFS),iwall_tmp(lx1*lelt*MAXSURFS)

      integer i,j,ix,iy,iz

      integer isf,sel,selno,smap,lmap
      integer cnt_rem
      integer totpts

      integer msg_id1,msg_id2

      integer ivlsum,irecv                ! function

      character(8) outfmt

      integer ierr,length,ip
      integer msgtag

!      write (outfmt, "(A1,I1,A6)") '(', 1, 'F10.6)'
!      if (nid.eq.0) write(*,*) 'outfmt', outfmt

      totpts=0
      if (nid.eq.0) then

!         LMAP
          write(6,'(5(A10,1x))') 'nid', 'SURF_LMAP','STAT_LMAP'
     $      ,'SURF_GMAP','SURF_OWN'     
         do isf=1,NSURFS
         do sel=1,SURF_MEMCOUNT(isf)
               write(6,'(5I10)') nid,SURF_LMAP(sel,isf),
     $   SURF_STAT_LMAP(sel,isf),SURF_GMAP(sel,isf),SURF_OWN(sel,isf)
         enddo
         enddo      

!         Arc Length            
!         do isf=1,NSURFS
!         do iel=1,SURF_LNUM(isf)
!               write(*,*) nid, ':', (SURF_ABM1(ix,iel,isf), ix=1,LX1)
!         enddo
!         enddo      

!        Integration Weights            
!         do isf=1,NSURFS
!         do sel=1,SURF_MEMCOUNT(isf)
!         do ix=1,1
!               iel = SURF_OBJ(1,sel,isf)
!               write(*,*) nid, ':', (SURF_BM1D(ix,1,iz,iel),
!     $                   iz=1,LZ1)
!         enddo
!         enddo   
!         enddo   

!         OWNERSHIP SURF/STATS
!         do isf=1,NSURFS
!         do sel=1,SURF_MEMCOUNT(isf)
!               iel = SURF_OBJ(1,sel,isf)
!               smap = SURF_LMAP(sel,isf)
!               lmap = STAT_LMAP(iel)    
!               write(*,*) nid,SURF_OWN(smap,isf),STAT_OWN(lmap)
!         enddo
!         enddo 

!         PROCID          
!         do sel=1,SURF_SND_CNT
!               write(*,*) nid,SURF_PROCID(sel)
!         enddo

!         SURF_LNUM          
!         write(*,*) nid,'SURF_LNUM:',SURF_LNUM,'MEMCOUNT:',
!     $                  SURF_MEMCOUNT

         open(unit=75, FILE='surf_data.debug')
         write(75,*) 'Proc/PROCID',nid,(SURF_PROCID(i),i=1,SURF_SND_CNT)
         
         do ip=1,np-1
!     hand sahking
            length = 1*isize
            ierr = ip
            call csend(ip,ierr,length,ip,0)
!     count number
            call crecv(ip,cnt_rem,length)
            write(75,*) 'Test #',ip,cnt_rem
!     get data
            if (cnt_rem.gt.0) then
               length = isize*cnt_rem
               call crecv(0,iwall,length)         ! PROCID 

               write(75,*) 'Proc/PROCID',ip, (iwall(i),i=1,cnt_rem)

            endif
!            counter = counter+cnt_rem
         enddo
         close(75)

      else
!     hand shaking
         length = 1*isize
         call crecv(nid,ierr,length)

!        LMAP 
        if (nid.eq.8) then    
         do isf=1,NSURFS
         do sel=1,SURF_MEMCOUNT(isf)   
               write(*,*) nid,SURF_LMAP(sel,isf),
     $   SURF_STAT_LMAP(sel,isf),SURF_GMAP(sel,isf),SURF_OWN(sel,isf)
         enddo
         enddo      
       endif
!         Arc Length            
!         do isf=1,NSURFS
!         do iel=1,SURF_LNUM(isf)
!               write(*,*) nid, ':', (SURF_ABM1(ix,iel,isf), ix=1,LX1)
!         enddo
!         enddo      

!        Integration Weights            
!         do isf=1,NSURFS
!         do sel=1,SURF_MEMCOUNT(isf)
!         do ix = 1,1
!               iel = SURF_OBJ(1,sel,isf)
!               write(*,*) nid, ':', (SURF_BM1D(ix,1,iz,iel),
!     $                   iz=1,LZ1)
!         enddo
!         enddo
!         enddo      
       
!         OWNERSHIP SURF/STATS          
!         do isf=1,NSURFS
!         do sel=1,SURF_MEMCOUNT(isf)
!               iel = SURF_OBJ(1,sel,isf)
!               smap = SURF_LMAP(sel,isf)
!               lmap = STAT_LMAP(iel)    
!               write(*,*) nid,SURF_OWN(smap,isf),STAT_OWN(lmap)
!         enddo
!         enddo 


!         PROCID          
!         do sel=1,SURF_SND_CNT
!               write(*,*) nid,SURF_PROCID(sel)
!         enddo

!         SURF_LNUM          
!         write(*,*) nid,'SURF_LNUM:',SURF_LNUM,'MEMCOUNT:',
!     $               SURF_MEMCOUNT

!     number of points
         call csend(nid,SURF_SND_CNT,length,0,0)

!     send data
         if(SURF_SND_CNT.gt.0) then

            length = isize*SURF_SND_CNT
            msgtag = 0
            call csend(msgtag,SURF_PROCID,length,0,0)

         endif

      endif
      call nekgsync()

      return
      end subroutine SURF_DEBUG_COMM            
!---------------------------------------------------------------------- 

      subroutine get_ref_ind(refx,refy,refz,refe,refp,x0,y0,z0)

!     indicies and proc id of the reference point defined for Cp

      implicit none
      
      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer ii,jj,kk,ee

      real x0,y0,z0,dmin,dp
      integer refx,refy,refz,refe,refp
      integer rix,riy,riz,rie,rip 
      save rix,riy,riz,rie,rip

      real tmp_tol
      parameter (tmp_tol = 1E-08)

      integer icalld
      save icalld
      data icalld /0/

      integer msgtag,len

      real glamin       ! function global absolute minimum

!---------------------------------------------------------------------- 
!     indicies of reference location
!     Only need to be calculated once.
!     That might change if meshes are moving! 
            
      if (icalld.eq.0) then
            icalld = 1

            dmin=1.0E+5
            do 20 ee=1,nelv
            do 20 kk=1,lz1
            do 20 jj=1,ly1
            do 20 ii=1,lx1
            dp=sqrt((XM1(ii,jj,kk,ee)-x0)**2+
     $     (YM1(ii,jj,kk,ee)-y0)**2 + (ZM1(ii,jj,kk,ee)-z0)**2)
            if (dp.lt.dmin) then
                  dmin=dp
                  rix=ii
                  riy=jj
                  riz=kk
                  rie=ee
            endif
   20       continue
             
            dp = glamin(dmin,1)
            rip = -1
   
            if (abs(dmin-dp).lt.tmp_tol) then
                  rip = nid               ! which process
!                  write(*,*) 'Query loc:', x0,y0
!                  write(*,*) 'Found point:',nid,xm1(tix,tiy,tiz,tie),
!     $                       YM1(tix,tiy,tiz,tie),dmin,dmin-dp
            endif

            len=1*ISIZE
            msgtag=0
            if (nid.eq.0) then
               if (rip.eq.-1) then
                    call crecv(msgtag,rip,len)         ! receive which proc has ref point.
               endif
            elseif (rip.ne.-1) then
               call csend(msgtag,rip,len,0,0)          ! send owner pid to 0
            endif

            call bcast(rip,len)                        ! broadcast owner pid                 
            if (nid.eq.0) then
               write(6,*) 'Process with reference point:,',rip
            endif
      endif

      refx=rix
      refy=riy
      refz=riz
      refe=rie
      refp=rip


      return
      end subroutine get_ref_ind

!----------------------------------------------------------------------

      subroutine get_normal_comp(svn1,svar,ndir)

      implicit none
      
      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'TOPOL_DEF'
      include 'TOPOL'         ! skpdat,eface1
      include 'SURF_STATS'

      real svn1(lx1,ly1,lz1,lelt),svar(lx1,ly1,lz1,lelt)
      integer isf,iel,f,sel
      integer ndir
      real sn
      integer pf,j1,j2,js1,jf1,jskip1,js2,jf2,jskip2
      integer cntr


      call copy(svn1,svar,lx1*ly1*lz1*nelt)

      do isf=1,NSURFS
        do sel=1,SURF_MEMCOUNT(isf)
          f   = SURF_OBJ(2,sel,isf)
          iel = SURF_OBJ(1,sel,isf)    

          call dsset(nx1,ny1,nz1)    ! set up counters
          pf     = eface1(f)     ! convert from preproc. notation
          js1    = skpdat(1,pf)
          jf1    = skpdat(2,pf)
          jskip1 = skpdat(3,pf)
          js2    = skpdat(4,pf)
          jf2    = skpdat(5,pf)
          jskip2 = skpdat(6,pf)
 
          cntr=0 
          do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
              cntr=cntr+1
              if (ndir.eq.1) then 
                sn = unx(cntr,1,f,iel)
              elseif (ndir.eq.2) then  
                sn = uny(cntr,1,f,iel)
              elseif (ndir.eq.3) then  
                sn = unz(cntr,1,f,iel)
              else
                sn = 1.
              endif

              svn1(j1,j2,1,iel) = sn*svar(j1,j2,1,iel) 

            enddo      ! j1
          enddo        ! j2
        enddo          ! sel 
      enddo            ! isf 

      return
      end subroutine get_normal_comp
!---------------------------------------------------------------------- 
      subroutine CALC_SURF_INTEGRAL(intgvar,lvar,f,e)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'TOPOL_DEF'
      include 'TOPOL'
      include 'SURF_STATS'

      real lvar(lx1,ly1,lz1,lelt),intgvar(3)
      integer i,e,f
      integer pf,j1,j2,js1,jf1,jskip1,js2,jf2,jskip2
     
      call dsset(nx1,ny1,nz1)    ! set up counters
      pf     = eface1(f)     ! convert from preproc. notation
      js1    = skpdat(1,pf)
      jf1    = skpdat(2,pf)
      jskip1 = skpdat(3,pf)
      js2    = skpdat(4,pf)
      jf2    = skpdat(5,pf)
      jskip2 = skpdat(6,pf)

      i=0 
      do j2=js2,jf2,jskip2
        do j1=js1,jf1,jskip1
          i=i+1
          intgvar(1) = intgvar(1) + lvar(j1,j2,1,e)*unx(i,1,f,e)*
     $                        area(i,1,f,e)   
          intgvar(2) = intgvar(2) + lvar(j1,j2,1,e)*uny(i,1,f,e)*
     $                        area(i,1,f,e)
          if (ndim.eq.3) then  
            intgvar(3) = intgvar(3) + lvar(j1,j2,1,e)*unz(i,1,f,e)*
     $                        area(i,1,f,e)
          endif
  
        enddo
      enddo
   
      return
      end subroutine CALC_SURF_INTEGRAL

!----------------------------------------------------------------------

      subroutine rbcastn(buf,len,sid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer sid
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,sid,nekcomm,ierr)

      return
      end
!-----------------------------------------------------------------------
