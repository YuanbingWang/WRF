subroutine da_read_obs_netcdf4giirs (iv, infile_tb)
   !--------------------------------------------------------
   !  Purpose: read in FY4 GIIRS Level-1 data in NETCDF4 format
   !           and form innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read the file twice
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2018/01/22 - Creation         Yuanbing Wang, NUIST/CAS
   !  To be devoloped: 1.time information; 2.dimension sequence
   !------------------------------------------------------------------------------
   
   implicit none

   character(len=*), intent(in)    :: infile_tb
   type(iv_type),    intent(inout) :: iv

! fixed parameter values
   integer,parameter::nfile_max = 7     ! each netcdf file contains
   real,parameter::scale_factor_tb=0.01          ! Maximum allowed NumberOfScans
   real,parameter::add_offset_tb=273.15          ! low  resolution pixel width   
   
! interface variable
   integer iret, rcode, ncid                      ! return status

! array data
   real(4), allocatable :: vlatitude(:)  !  value for latitude 
   real(4), allocatable :: vlongitude(:) !  value for longitude

   real(4), allocatable    :: tbb(:,:,:)  ! tb for band 7-16
   real(4), allocatable    :: sat_zenith(:,:)
   
   real(r_kind),parameter  :: tbmin  = 50._r_kind
   real(r_kind),parameter  :: tbmax  = 550._r_kind

   real(kind=8)                   :: obs_time
   type (datalink_type),pointer   :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   integer(i_kind)    :: idate5(6)
   character(len=80) :: filename,str1,str2

   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id
   real(r_kind)      :: tb, crit
   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse

   integer           :: i,j,k,l,m,n, ifile, landsea_mask
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth
   integer(i_kind)  :: num_giirs_local, num_giirs_global, num_giirs_used, num_giirs_thinned
   integer(i_kind)  :: num_giirs_used_tmp, num_giirs_file
   integer(i_kind)  :: num_giirs_local_local, num_giirs_global_local, num_giirs_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1
   integer          :: nchan,nsatobs,isatobs,ichannels
   integer          :: satobsid,latid,lonid,tbb_id,sazid
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max)
   logical          :: fexist
   
! Allocatable arrays
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)

   if (trace_use) call da_trace_entry("da_read_obs_netcdf4giirs")

!  0.0  Initialize variables
!-----------------------------------
   head_allocated = .false.
   platform_id  = 52  ! Table-2 Col 1 corresponding to 'fy4'
   satellite_id = 1   ! Table-2 Col 3
   sensor_id    = 98  ! Table-3 Col 2 corresponding to 'giirs'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_giirs_file    = 0
   num_giirs_local   = 0
   num_giirs_global  = 0
   num_giirs_used    = 0
   num_giirs_thinned = 0

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
          (/"The combination of Satellite_Id and Sensor_Id for AHI is not found"/))
      if (trace_use) call da_trace_exit("da_read_obs_netcdf4giirs")
      return
   end if

   nchan = iv%instid(inst)%nchan
   nchan=961
   write(unit=stdout,fmt=*)'AHI nchan: ',nchan
   allocate(data_all(1:nchan))

! 1.0 Assign file names and prepare to read ahi files
!-------------------------------------------------------------------------
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   
   ! first check if giirs nc file is available
   filename1 = trim(infile_tb)
   inquire (file=filename1, exist=fexist)
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
   else
      ! check if netcdf4 files are available for multiple input files
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
         (/"No valid AHI file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_netcdf4giirs")
      return
   end if 

  !open the data info file 
   open(unit=1990,file='giirs_info',status='old',iostat=iret)
   if(iret /= 0)then
      call da_warning(__FILE__,__LINE__,(/"data_info file read error"/))
   endif
   read(1990,*) 
   read(1990,*) 

   infile_loop:  do ifile = 1, nfile
      num_giirs_file_local    = 0
      num_giirs_local_local   = 0
      num_giirs_global_local  = 0

   ! open NETCDF4 L1 file for read
      iret = nf90_open(fname_tb(ifile), NF90_NOWRITE, ncid)
      if(iret /= 0)then
         call da_warning(__FILE__,__LINE__, &
              (/"Cannot open NETCDF4 file "//trim(fname_tb(ifile))/))
         cycle infile_loop
      endif

   ! read dimensions: 	nsatobs
      iret = nf90_inq_dimid(ncid, "ncl1", satobsid)
      iret = nf90_inquire_dimension(ncid, satobsid, len=nsatobs)
	  
      write(unit=stdout,fmt=*)'nsatobs: ',nsatobs   
 
 ! read array: time
      iret = nf90_get_att(ncid, nf90_global, "id", filename)
      if(iret /= 0)then
         call da_warning(__FILE__,__LINE__, &
             (/"NETCDF4 read error for: observation date"/))
      endif	  
	  read(filename,"(A7,I4,I2,I2,A1,I2,I2)") str1,idate5(1),idate5(2),idate5(3), &
						                      str2,idate5(4),idate5(5)
	  idate5(6)=00	  
	  write(unit=stdout,fmt=*)'observation date: ', idate5
   
   ! read array: latlon
   ! read lat
	  iret = nf90_inq_varid(ncid, 'latMW', latid)	  
	  allocate( vlatitude(nsatobs))
      iret = nf90_get_var(ncid, latid, vlatitude)
      if(iret /= 0)then
         call da_warning(__FILE__,__LINE__, &
             (/"NETCDF4 read error for: Latitude of Observation Point"/))
      endif
   ! read lon
	  iret = nf90_inq_varid(ncid, 'lonMW', lonid)
	  allocate( vlongitude(nsatobs))
      iret = nf90_get_var(ncid, lonid, vlongitude)
      if(iret /= 0)then
          call da_warning(__FILE__,__LINE__, &
              (/"NETCDF4 read error for: Longitude of Observation Point"/))
          call da_trace_exit("da_read_obs_netcdf4giirs")
      endif
   ! sample display
      !write(unit=stdout,fmt=*)'latitude,longitude(pixel=1,scan=1): ',vlatitude(1),vlongitude(1) 	  

   ! read array: tb for band 7-16
   ! read	 
      allocate( tbb(nsatobs,nchan)) 
      do k=1,nchan	  
		iret = nf90_inq_varid(ncid, 'ES_RealLW_TBB', tbb_id)	  
		iret = nf90_get_var(ncid, tbb_id, tbb(:,k))
		if(iret /= 0)then
		  call da_warning(__FILE__,__LINE__, &
			   (/"NETCDF4 read error for: Brightness Temperature"/))
		endif	  		
      end do
   ! sample display
	  !write(unit=stdout,fmt=*)'tbb(pixel=1,scan=1,chan=',k,'): ', tbb(1,1,1)	  
	 
   ! read array: satellite zenith angle
   ! read
	  iret = nf90_inq_varid(ncid, 'SaZMW', sazid)	  
	  allocate( sat_zenith(nsatobs))
	  iret = nf90_get_var(ncid, sazid, sat_zenith)
      if(iret /= 0)then
         call da_warning(__FILE__,__LINE__, &
            (/"NETCDF4 read error for: satellite zenith angle"/))
      endif	  	  
   ! sample display
      !write(unit=stdout,fmt=*)'satellite zenith angle(pixel=1,scan=1): ',sat_zenith(1,1)
 
   ! close infile_tb file 		 
      iret = nf90_close(ncid) 

	 !read date information	  
	  read(1990,*) idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),idate5(6)	  

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
      scan_loop:     do isatobs=1, nsatobs
	  
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)
         if (obs_time<time_slots(0) .or. obs_time>=time_slots(num_fgat_time)) cycle scan_loop
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and. obs_time < time_slots(ifgat)) exit
         end do
			
			if ( sat_zenith(isatobs) > 65.0 ) cycle fov_loop
			
            num_giirs_file       = num_giirs_file + 1
            num_giirs_file_local = num_giirs_file_local + 1
            info%lat  =  vlatitude(isatobs)
            info%lon  =  vlongitude(isatobs)

            call da_llxy (info, loc, outside, outside_all)
            if (outside_all) cycle fov_loop

            num_giirs_global       = num_giirs_global + 1
            num_giirs_global_local = num_giirs_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE

            num_giirs_local       = num_giirs_local + 1
            num_giirs_local_local = num_giirs_local_local + 1
            write(unit=info%date_char, &
            fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
               idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
               ':', idate5(5), ':', idate5(6)
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
                  num_giirs_thinned = num_giirs_thinned+1
                  cycle fov_loop
               end if
            end if

            num_giirs_used = num_giirs_used + 1
            data_all = missing_r		

            do k=1,nchan
               tb = tbb(isatobs,k)
               if( tb < tbmin .or. tb > tbmax ) tb = missing_r
               data_all(k)= tb
            end do
			
! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            p%info             = info
            p%loc              = loc
            p%landsea_mask     = 1
            p%scanpos          = sat_zenith(isatobs)
            p%satzen           = sat_zenith(isatobs)
            p%satazi           = 0
            p%solzen           = 0
            p%solazi           = 0
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
      end do scan_loop

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_file    : ',num_giirs_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_global  : ',num_giirs_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_giirs_local   : ',num_giirs_local_local

      deallocate(vlatitude)
      deallocate(vlongitude)
      deallocate(tbb) 
      deallocate(sat_zenith)
   
   end do infile_loop  

   close(1990) !close date information file   
   deallocate(data_all) ! Deallocate data arrays  

   if (thinning .and. num_giirs_global > 0 ) then
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
      num_giirs_used_tmp = num_giirs_used
      do j = 1, num_giirs_used_tmp
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
            num_giirs_thinned = num_giirs_thinned + 1
            num_giirs_used = num_giirs_used - 1
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

   iv%total_rad_pixel   = iv%total_rad_pixel + num_giirs_used
   iv%total_rad_channel = iv%total_rad_channel + num_giirs_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_giirs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_giirs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'AHI data counts: '
   write(stdout,fmt='(a,i7)') ' In file: ',num_giirs_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_giirs_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_giirs_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_giirs_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_giirs_thinned

!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------

   if (num_giirs_used > 0) then
      iv%instid(inst)%num_rad  = num_giirs_used
      iv%instid(inst)%info%nlocal = num_giirs_used
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
         'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
      call da_allocate_rad_iv (inst, nchan, iv)
   end if

!  6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------
   p => head

   do n = 1, num_giirs_used
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

   if (trace_use) call da_trace_exit("da_read_obs_netcdf4giirs")
   
end subroutine da_read_obs_netcdf4giirs
