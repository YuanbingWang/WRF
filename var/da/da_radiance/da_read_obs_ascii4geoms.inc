subroutine da_read_obs_ascii4geoms (iv, infile_tb)
   !--------------------------------------------------------
   !  Purpose: read in geoms(rttov id: 17,0,19, ATMS) Level-1 data in ASCII format
   !           and form innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read the file twice
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2019/08/15 - Creation         Yuanbing Wang, NUIST/CAS
   !------------------------------------------------------------------------------
   
   implicit none

   character(len=*), intent(in)    :: infile_tb
   type(iv_type),    intent(inout) :: iv

! fixed parameter values
   integer,parameter::nfile_max = 7     ! each netcdf file contains 
   
! interface variable
   integer iret, rcode, ncid                      ! return status

! array data
   real(4), allocatable :: vlatitude(:,:)  !  value for latitude 
   real(4), allocatable :: vlongitude(:,:) !  value for longitude

   real(4), allocatable    :: tbb(:,:,:)  ! tb for band 1-22
   real(4), allocatable    :: sat_zenith(:,:),sat_azimuth(:,:),sun_zenith(:,:),sun_azimuth(:,:)
   
   real(r_kind),parameter  :: tbmin  = 50._r_kind
   real(r_kind),parameter  :: tbmax  = 550._r_kind

   real(kind=8)                   :: obs_time
   type (datalink_type),pointer   :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   integer(i_kind), allocatable   :: idate5(:,:)

   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id
   real(r_kind)      :: tb, crit
   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse

   integer           :: i,j,k,l,m,n,ifile,landsea_mask,sfcflag
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth
   integer(i_kind)  :: num_geoms_local, num_geoms_global, num_geoms_used, num_geoms_thinned
   integer(i_kind)  :: num_geoms_used_tmp, num_geoms_file
   integer(i_kind)  :: num_geoms_local_local, num_geoms_global_local, num_geoms_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1
   integer          :: nchan,nlongitude,nlatitude,ifov,iscan,ichannels
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max)
   logical          :: fexist
   
! Allocatable arrays
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)

   if (trace_use) call da_trace_entry("da_read_obs_ascii4geoms")

!  0.0  Initialize variables
!-----------------------------------
   head_allocated = .false.
   platform_id  = 17  ! Table-2 Col 1 corresponding to 'fy3'
   satellite_id = 0   ! Table-2 Col 3
   sensor_id    = 19  ! Table-3 Col 2 corresponding to 'mwhs-2'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_geoms_file    = 0
   num_geoms_local   = 0
   num_geoms_global  = 0
   num_geoms_used    = 0
   num_geoms_thinned = 0

   do i = 1, rtminit_nsensor
      if (platform_id  == rtminit_platform(i) &
          .and. satellite_id == rtminit_satid(i)    &
          .and. sensor_id    == rtminit_sensor(i)) then
         inst = i
         exit
      end if
   end do
   if (inst == 0) then
      call da_warning(__FILE__,__LINE__, &
          (/"The combination of Satellite_Id and Sensor_Id for geoms is not found"/))
      if (trace_use) call da_trace_exit("da_read_obs_ascii4geoms")
      return
   end if

   nchan = iv%instid(inst)%nchan
   write(unit=stdout,fmt=*)'GEOMS nchan: ',nchan
   allocate(data_all(1:nchan))

! 1.0 Assign file names and prepare to read geoms files
!-------------------------------------------------------------------------
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   
   ! first check if geoms file is available
   filename1 = trim(infile_tb)
   inquire (file=filename1, exist=fexist)
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
   else
      ! check if ascii files are available for multiple input files
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max
         write(filename1,fmt='(A,2I1)') trim(infile_tb),0,i
         inquire (file=trim(filename1), exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_tb(nfile)  = filename1
         else
            exit
         end if
      end do
   end if

   if ( nfile == 0 ) then
      call da_warning(__FILE__,__LINE__, &
         (/"No valid geoms file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_ascii4geoms")
      return
   end if 

   infile_loop:  do ifile = 1, nfile
      num_geoms_file_local    = 0
      num_geoms_local_local   = 0
      num_geoms_global_local  = 0

   ! open ASCII L1 file for read
      open(unit=2019,file=fname_tb(ifile),status='old',iostat=iret)
      if(iret /= 0)then
         call da_warning(__FILE__,__LINE__, &
              (/"Cannot open ASCII file "//trim(fname_tb(ifile))/))
         cycle infile_loop
      endif

   ! read dimensions: 	nlatitude and nlongitude
	  read(2019,*) nlatitude,nlongitude 
      write(unit=stdout,fmt=*)'nlongitude,nlatitude: ',nlongitude,nlatitude   	  

   ! read array: tb for band 1-22
   ! read	 
      allocate(tbb(nlongitude,nlatitude,nchan)) 
	  allocate(vlatitude(nlongitude,nlatitude))
	  allocate(vlongitude(nlongitude,nlatitude))
	  allocate(sat_zenith(nlongitude,nlatitude))
	  allocate(sun_zenith(nlongitude,nlatitude))
	  allocate(sat_azimuth(nlongitude,nlatitude))
	  allocate(sun_azimuth(nlongitude,nlatitude))
	  allocate(idate5(6,nlatitude))
	  do j=1,nlatitude
		do i=1,nlongitude
           read(2019,*) sun_azimuth(i,j),sun_zenith(i,j),sat_azimuth(i,j),sat_zenith(i,j),&
						idate5(1,j),idate5(2,j),idate5(3,j),idate5(4,j),idate5(5,j),idate5(6,j),&
						sfcflag,vlatitude(i,j),vlongitude(i,j),tbb(i,j,1),tbb(i,j,2),tbb(i,j,3),&
						tbb(i,j,4),tbb(i,j,5),tbb(i,j,6),tbb(i,j,7),tbb(i,j,8),tbb(i,j,9),tbb(i,j,10),&
						tbb(i,j,11),tbb(i,j,12),tbb(i,j,13),tbb(i,j,14),tbb(i,j,15),,tbb(i,j,16),&
						tbb(i,j,17),tbb(i,j,18),tbb(i,j,19),tbb(i,j,20),tbb(i,j,21),,tbb(i,j,22),
		end do
	  end do	
 
   ! close infile_tb file 		 
      close(2019) 

! 2.0 Loop to read netcdf and assign information to a sequential structure
!-------------------------------------------------------------------------   
   ! Allocate arrays to hold data
      if ( .not. head_allocated ) then
         allocate (head)
         nullify  ( head % next )
         p => head
         head_allocated = .true.
      end if

   ! start scan_loop
      scan_loop:     do iscan=1, nlatitude 
         call da_get_julian_time(idate5(1,iscan),idate5(2,iscan),idate5(3,iscan),&
		                         idate5(4,iscan),idate5(5,iscan),obs_time)
         if (obs_time<time_slots(0) .or. obs_time>=time_slots(num_fgat_time)) cycle scan_loop
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and. obs_time < time_slots(ifgat)) exit
         end do

      ! start fov_loop
         fov_loop:      do ifov=1, nlongitude
			
			if ( sat_zenith(ifov,iscan) > 65.0 ) cycle fov_loop
			
            num_geoms_file       = num_geoms_file + 1
            num_geoms_file_local = num_geoms_file_local + 1
            info%lat  =  vlatitude(ifov,iscan)
            info%lon  =  vlongitude(ifov,iscan)

            call da_llxy (info, loc, outside, outside_all)
            if (outside_all) cycle fov_loop

            num_geoms_global       = num_geoms_global + 1
            num_geoms_global_local = num_geoms_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE

            num_geoms_local       = num_geoms_local + 1
            num_geoms_local_local = num_geoms_local_local + 1
            write(unit=info%date_char, &
            fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
               idate5(1,iscan), '-', idate5(2,iscan), '-', idate5(3,iscan), '_', idate5(4,iscan), &
               ':', idate5(5,iscan), ':', idate5(6,iscan)
            info%elv = 0.0

! 3.0  Make Thinning
! Map obs to thinning grid
!-------------------------------------------------------------------
            if (thinning) then
               dlat_earth = info%lat !degree
               dlon_earth = info%lon
               if (dlon_earth<zero)  dlon_earth = dlon_earth+r360
               if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
               dlat_earth = dlat_earth*deg2rad !radian
               dlon_earth = dlon_earth*deg2rad
               crit = 1.
               call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
               if (.not. iuse) then
                  num_geoms_thinned = num_geoms_thinned+1
                  cycle fov_loop
               end if
            end if

            num_geoms_used = num_geoms_used + 1
            data_all = missing_r		

            do k=1,nchan
               tb = tbb(ifov,iscan,k)
               if( tb < tbmin .or. tb > tbmax ) tb = missing_r
               data_all(k)= tb
            end do
			
! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            p%info             = info
            p%loc              = loc
            p%landsea_mask     = 1
            p%scanpos          = ifov
            p%satzen           = sat_zenith(ifov,iscan)
            p%satazi           = sat_azimuth(ifov,iscan)
            p%solzen           = sun_zenith(ifov,iscan)
            p%solazi           = sun_azimuth(ifov,iscan)
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
         end do fov_loop
      end do scan_loop

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_geoms_file    : ',num_geoms_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_geoms_global  : ',num_geoms_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_geoms_local   : ',num_geoms_local_local

      deallocate(vlatitude)
      deallocate(vlongitude)
      deallocate(tbb)
      deallocate(idate5)  
      deallocate(sat_zenith)
	  deallocate(sun_zenith)
	  deallocate(sat_azimuth)
	  deallocate(sun_azimuth)
   
   end do infile_loop  

   deallocate(data_all) ! Deallocate data arrays  

   if (thinning .and. num_geoms_global > 0 ) then
#ifdef DM_PARALLEL
   ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         j = j + thinning_grid(inst,ifgat)%itxmax
      end do 

      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            in(j) = thinning_grid(inst,ifgat)%score_crit(i)
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            if ( ABS(out(j)-thinning_grid(inst,ifgat)%score_crit(i)) > 1.0E-10 ) &
            thinning_grid(inst,ifgat)%ibest_obs(i) = 0
         end do
      end do

      deallocate( in  )
      deallocate( out )

#endif

   ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_geoms_used_tmp = num_geoms_used
      do j = 1, num_geoms_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

      ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_geoms_thinned = num_geoms_thinned + 1
            num_geoms_used = num_geoms_used - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_geoms_used
   iv%total_rad_channel = iv%total_rad_channel + num_geoms_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_geoms_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_geoms_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'geoms data counts: '
   write(stdout,fmt='(a,i7)') ' In file: ',num_geoms_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_geoms_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_geoms_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_geoms_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_geoms_thinned

!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------

   if (num_geoms_used > 0) then
      iv%instid(inst)%num_rad  = num_geoms_used
      iv%instid(inst)%info%nlocal = num_geoms_used
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
         'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
      call da_allocate_rad_iv (inst, nchan, iv)
   end if

!  6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------
   p => head

   do n = 1, num_geoms_used
      i = p%sensor_index 
      call da_initialize_rad_iv (i, n, iv, p)
      current => p
      p => p%next
   ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
   deallocate ( p )
   deallocate (ptotal)

   if (trace_use) call da_trace_exit("da_read_obs_ascii4geoms")
   
end subroutine da_read_obs_ascii4geoms
