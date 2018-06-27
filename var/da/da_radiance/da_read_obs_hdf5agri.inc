subroutine da_read_obs_hdf5agri (iv,infile_tb,infile_clp)
   !--------------------------------------------------------
   !  Purpose: read in CMA AGRI Level-1B HDF5 and Level-2 data in NETCDF format
   !           and form innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read the file twice
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2018/06/10 - Creation         Yuanbing Wang, NUIST/CAS
   !------------------------------------------------------------------------------

   implicit none

   character(len=*), intent(in)    :: infile_tb, infile_clp
   type(iv_type),    intent(inout) :: iv

#if defined(HDF5)
! fixed parameter values
   integer,parameter::nlatitude=2748          ! Maximum allowed NumberOfScans
   integer,parameter::nlongitude=2748          ! low  resolution pixel width
   integer,parameter::time_dims=6       ! Time dimension
   integer,parameter::nfile_max = 8     ! each hdf file contains ~50min of data
                                        ! at most 8 files for a 6-h time window
! interface variable
   integer iret                         ! return status
   integer(HID_T) fhnd1       ! file handle
   integer(HID_T) dhnd1         ! dataset handle
   integer(HSIZE_T) sz1(2)              ! array size 1
   
! array data
   real(8) :: vlatitude(nlongitude,nlatitude)  !  value for latitude 
   real(8) :: vlongitude(nlongitude,nlatitude) !  value for longitude

   real(4) :: tb07(nlongitude,nlatitude)  ! tb for band 7
   real(4) :: tb08(nlongitude,nlatitude)  ! tb for band 8
   real(4) :: tb09(nlongitude,nlatitude)  ! tb for band 9
   real(4) :: tb10(nlongitude,nlatitude)  ! tb for band 10
   real(4) :: tb11(nlongitude,nlatitude)  ! tb for band 11
   real(4) :: tb12(nlongitude,nlatitude)  ! tb for band 12
   real(4) :: tb13(nlongitude,nlatitude)  ! tb for band 13
   real(4) :: tb14(nlongitude,nlatitude)  ! tb for band 14

   real(4) :: sat_zenith(nlongitude,nlatitude)  ! satellite_zenith_angle
   integer(4)::ColumnNumber(nlongitude,nlatitude),LineNumber(nlongitude,nlatitude)
   integer(4) :: cloud_mask(nlongitude,nlatitude)     !obs cloud mask
   
   real(4) :: latlonx,latlony,latlonsd1,latlonsd2,latlonsd3,latlonsd
   real(4) :: latlonsn,latlons1,latlons2,latlons3,latlonsxy
			  
   real(r_kind),parameter  :: tbmin  = 50._r_kind
   real(r_kind),parameter  :: tbmax  = 550._r_kind
   real(r_kind),parameter  :: tb_scale  = 100._r_kind

   real(kind=8)                   :: obs_time
   type (datalink_type),pointer   :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   integer(i_kind)   :: idate5(6)
   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id,clmid,ncid
   real(r_kind)      :: tb, crit
   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse

   integer           :: i,j,k,l,m,n, ifile,landsea_mask
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth
   integer(i_kind)  :: num_agri_local, num_agri_global, num_agri_used, num_agri_thinned
   integer(i_kind)  :: num_agri_used_tmp, num_agri_file
   integer(i_kind)  :: num_agri_local_local, num_agri_global_local, num_agri_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1, filename2
   integer          :: nchan,ilongitude,ilatitude,ichannels
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max)
   character(80)    :: fname_clp(nfile_max)
   logical          :: fexist, got_clp_file

! Allocatable arrays
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)

   if (trace_use) call da_trace_entry("da_read_obs_hdf5agri")

!  0.0  Initialize variables
!-----------------------------------
   head_allocated = .false.
   platform_id  = 52  ! Table-2 Col 1 corresponding to 'himawari'
   satellite_id = 1   ! Table-2 Col 3
   sensor_id    = 99  ! Table-3 Col 2 corresponding to 'agri'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_agri_file    = 0
   num_agri_local   = 0
   num_agri_global  = 0
   num_agri_used    = 0
   num_agri_thinned = 0
   
   sz1(1)=nlongitude
   sz1(2)=nlatitude  

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
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if

! Initialize HDF5 library and Fortran90 interface
   call H5open_f(iret)
   if(iret.lt.0)then
      call da_warning(__FILE__,__LINE__,(/"Problems initializing HDF5 Lib, can't read AHI data."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   endif

   nchan = iv%instid(inst)%nchan
   write(unit=stdout,fmt=*)'AHI nchan: ',nchan
   allocate(data_all(1:nchan))

! 1.0 Assign file names and prepare to read agri files
!-------------------------------------------------------------------------
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   ! first check if hdf file is available
   filename1 = trim(infile_tb)
   filename2 = trim(infile_clp)
   inquire (file=filename1, exist=fexist)
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
      fname_clp(nfile) = filename2
   else
      ! check if L1SGRTBR-0x.h5 is available for multiple input files
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max
         write(filename1,fmt='(A,A,I2.2,A)') trim(infile_tb),'-',i
         write(filename2,fmt='(A,A,I2.2,A)') trim(infile_clp),'-',i
         inquire (file=filename1, exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_tb(nfile)  = filename1
            fname_clp(nfile) = filename2
         else
            exit
         end if
      end do
   end if

   if ( nfile == 0 ) then
      call da_warning(__FILE__,__LINE__,(/"No valid AHI file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if

  !open the data info file 
   open(unit=1990,file='agri_info',status='old',iostat=iret)
   if(iret /= 0)then
      call da_warning(__FILE__,__LINE__,(/"data_info file read error"/))
   endif
   read(1990,*) 
   read(1990,*)
   read(1990,*)
   read(1990,*)
   read(1990,*)
   read(1990,*)
   read(1990,*)
   
   infile_loop:  do ifile = 1, nfile
      num_agri_file_local    = 0
      num_agri_local_local   = 0
      num_agri_global_local  = 0

   ! open infile_tb HDF5 file for read
      call H5Fopen_f(fname_tb(ifile),H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"Cannot open HDF5 file"//trim(fname_tb(ifile))/))
         cycle infile_loop
      endif  

   ! read tb for band 7
      call H5Dopen_f(fhnd1,'NOMChannel07',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb07,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 7"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb07(pixel=1,scan=1): ',tb07(1,1)	  
	  
   ! read tb for band 8
      call H5Dopen_f(fhnd1,'NOMChannel08',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb08,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 8"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb08(pixel=1,scan=1): ',tb08(1,1)

   ! read tb for band 9
      call H5Dopen_f(fhnd1,'NOMChannel09',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb09,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 9"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb09(pixel=1,scan=1): ',tb09(1,1)

   ! read tb for band 10
      call H5Dopen_f(fhnd1,'NOMChannel10',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb10,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 10"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb10(pixel=1,scan=1): ',tb10(1,1)

   ! read tb for band 11
      call H5Dopen_f(fhnd1,'NOMChannel11',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb11,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 11"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb11(pixel=1,scan=1): ',tb11(1,1)		 

   ! read tb for band 12
      call H5Dopen_f(fhnd1,'NOMChannel12',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb12,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 12"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb12(pixel=1,scan=1): ',tb12(1,1)		

   ! read tb for band 13
      call H5Dopen_f(fhnd1,'NOMChannel13',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb13,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: 13"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb13(pixel=1,scan=1): ',tb13(1,1)	 

   ! read tb for band 14
      call H5Dopen_f(fhnd1, 'NOMChannel14',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_REAL,tb14,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Band 14"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*) 'tb14(pixel=1,scan=1): ',tb14(1,1)	

   ! close infile_tb and HDF5
      call H5Fclose_f(fhnd1,iret)	  

	  
   ! open infile for reading geo data
      call H5Fopen_f('L1AGRIGEO',H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)
	  
   ! read array: satellite_zenith_angle
      call H5Dopen_f(fhnd1,'NOMSatelliteZenith',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_IEEE_F32BE,sat_zenith,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: satellite_zenith_angle"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*)'sat_zenith(pixel=1,scan=1): ',sat_zenith(1,1)
	  
   ! read array: LineNumber
      call H5Dopen_f(fhnd1,'LineNumber',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_SHORT,LineNumber,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: LineNumber"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*)'LineNumber(pixel=1,scan=1): ',LineNumber(1,1)	 

   ! read array: ColumnNumber
      call H5Dopen_f(fhnd1,'ColumnNumber',dhnd1,iret)
      call H5Dread_f(dhnd1,H5T_NATIVE_SHORT,ColumnNumber,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: ColumnNumber"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*)'ColumnNumber(pixel=1,scan=1): ',ColumnNumber(1,1)	  

   ! close infile_geo and HDF5	  
	  call H5Fclose_f(fhnd1,iret)

	  
     !open infile_clp file 
      got_clp_file = .false.
      iret = nf90_open(fname_clp(ifile), NF90_NOWRITE, ncid)
      if ( iret == 0 ) then
         got_clp_file = .true.
      end if
	  
      if ( got_clp_file ) then
		! read array: cloud mask 
          rcode = nf90_inq_varid(ncid, 'CLM', clmid)
		  iret = nf90_get_var(ncid, clmid, cloud_mask)
          if(iret /= 0)then
            call da_warning(__FILE__,__LINE__,(/"NETCDF4 read error for: cloud_mask data"/))
          endif
        ! sample display
          !write(unit=stdout,fmt=*)'cloud_mask(pixel=1,scan=1): ',cloud_mask(1,1)		 
		! close infile_clp file 		 
		  iret = nf90_close(ncid)		 	 
      end if 
	  
   ! read binary latlon
      open(520,file='L1AGRILATLON',access='direct',form='unformatted', &
	           recl=8,convert='little_endian',iostat=iret) 
      if(iret.lt.0)then
         call da_warning(__FILE__,__LINE__,(/"HDF5 read error for: Latitude"/))
      end if	  
	  k = 0
	  do j = 1,2748
		do i = 1,2748
			read(520,rec=k) vlongitude(i,j)
			k = k + 1
			read(520,rec=k) vlatitude(i,j)
		end do
	  end do
      close(520)  
   ! sample display
      write(unit=stdout,fmt=*)'latitude,longitude(pixel=1,scan=1): ',vlatitude(1,1),vlongitude(1,1)	  

	 !read date information	  
	  read(1990,*) idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),idate5(6)
	  
! 2.0 Loop to read hdf file and assign information to a sequential structure
!-------------------------------------------------------------------------

   ! Allocate arrays to hold data
      if ( .not. head_allocated ) then
         allocate (head)
         nullify  ( head % next )
         p => head
         head_allocated = .true.
      end if

   ! start scan_loop
      scan_loop:     do ilatitude=1, nlatitude	 

         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)
		 
         if ( obs_time < time_slots(0) .or. obs_time >= time_slots(num_fgat_time) ) cycle scan_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and. obs_time  < time_slots(ifgat) ) exit
         end do

      ! start fov_loop: longitude
         fov_loop:      do ilongitude=1, nlongitude
			
			if ( sat_zenith(ilongitude,ilatitude) > 65.0 ) cycle fov_loop
			
            num_agri_file       = num_agri_file + 1
            num_agri_file_local = num_agri_file_local + 1
			
			!caculation the latitude and longitude
			latlonx=(ColumnNumber(ilongitude,ilatitude)-1373.5)*(65536*3.1415926)/(180*10233137)
			latlony=(LineNumber(ilongitude,ilatitude)-1373.5)*(65536*3.1415926)/(180*10233137)
			latlonsd1=cos(latlony)**2+sin(latlony)**2*1.0067395
			latlonsd2=42164*cos(latlonx)*cos(latlony)
			latlonsd3=latlonsd2**2-(42164**2-6378.137**2)*(latlonsd1)
			latlonsd=sqrt(latlonsd3)
			latlonsn=(42164*cos(latlonx)*cos(latlony)-latlonsd)/(latlonsd00)
			latlons1=42164-latlonsn*cos(latlonx)*cos(latlony)
			latlons2=latlonsn*sin(latlonx)*cos(latlony)
			latlons3=-latlonsn*sin(latlony)
			latlonsxy=sqrt(latlons1**2+latlons2**2)
			vlatitude(ilongitude,ilatitude)=57.2957805*atand(latlons3/latlonsxy*1.0067395)
			vlongitude(ilongitude,ilatitude)=57.2957805*atand(latlons2/latlons1)+104.7
			
            info%lat  =  vlatitude(ilongitude,ilatitude)
            info%lon  =  vlongitude(ilongitude,ilatitude)

            call da_llxy (info, loc, outside, outside_all)
            if (outside_all) cycle fov_loop

            num_agri_global       = num_agri_global + 1
            num_agri_global_local = num_agri_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE

            num_agri_local       = num_agri_local + 1
            num_agri_local_local = num_agri_local_local + 1
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
                  num_agri_thinned = num_agri_thinned+1
                  cycle fov_loop
               end if
            end if

            num_agri_used = num_agri_used + 1
            data_all = missing_r
			
			tb_scale = 1.0

            tb = tb07(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(1)= tb 

            tb = tb08(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(2)= tb 
			
			tb = tb09(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(3)= tb 

            tb = tb10(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(4)= tb 

            tb = tb11(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(5)= tb
			
			tb = tb12(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(6)= tb

            tb = tb13(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(7)= tb

            tb = tb14(ilongitude,ilatitude) / tb_scale
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(8)= tb		

! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            p%info             = info
            p%loc              = loc
            p%landsea_mask     = 1
            p%scanpos          = ilongitude !nint(sat_zenith(ilongitude,ilatitude))+1.001_r_kind !
            p%satzen           = sat_zenith(ilongitude,ilatitude)
            p%satazi           = 0
            p%solzen           = 0
            p%solazi           = 0
			p%cloudflag        = cloud_mask(ilongitude,ilatitude)
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
         end do fov_loop
      end do scan_loop

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_file    : ',num_agri_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_global  : ',num_agri_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_local   : ',num_agri_local_local
   end do infile_loop
 
   close(1990) !close date information file
   call H5close_f(iret)
   deallocate(data_all) ! Deallocate data arrays

   if (thinning .and. num_agri_global > 0 ) then
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
      num_agri_used_tmp = num_agri_used
      do j = 1, num_agri_used_tmp
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
            num_agri_thinned = num_agri_thinned + 1
            num_agri_used = num_agri_used - 1
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

   iv%total_rad_pixel   = iv%total_rad_pixel + num_agri_used
   iv%total_rad_channel = iv%total_rad_channel + num_agri_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_agri_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_agri_global

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
   write(stdout,fmt='(a,i7)') ' In file: ',num_agri_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_agri_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_agri_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_agri_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_agri_thinned


!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------

   if (num_agri_used > 0) then
      iv%instid(inst)%num_rad  = num_agri_used
      iv%instid(inst)%info%nlocal = num_agri_used
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
         'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
      call da_allocate_rad_iv (inst, nchan, iv)
   end if

!  6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------
   p => head

   do n = 1, num_agri_used
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

   if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
#else
   call da_error(__FILE__,__LINE__,(/"Needs to be compiled with HDF5 library"/))
#endif
end subroutine da_read_obs_hdf5agri
