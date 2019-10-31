subroutine da_read_obs_hdf5agri (iv, infile_tb, infile_geo, infile_clm)        !da_setup_radiance_structures(iv, 'agri', 'geo','clm')
   !--------------------------------------------------------
   !  Purpose: read in FY4 AGRI data in HDF5 format
   !           and form innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read the file twice
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !------------------------------------------------------------------------------

   implicit none

   character(len=*), intent(in)    :: infile_tb,infile_geo,infile_clm                                                
   type(iv_type),    intent(inout) :: iv

#if defined(HDF5)                                                                             
! fixed parameter values
   integer,parameter::nscan=2748        ! Maximum allowed NumberOfScans                          
   integer,parameter::fov=2748          ! high resolution pixel width                            
   integer,parameter::time_dims=6       ! Time dimension                                         
   integer,parameter::nfile_max = 7    ! each hdf file contains ~15min of data                  
                                        ! at most 12 files for a 6-h time window
   integer,parameter::h=42164
   integer,parameter::CFAC=10233137
   integer,parameter::LFAC=10233137
   real,parameter::pai=3.1415926536
   real,parameter::ea=6378.137
   real,parameter::eb=6356.7523
   real,parameter::cetlon=104.7
   real,parameter::COFF=1375.5
   real,parameter::LOFF=1375.5

! interface variable                    
   integer iret
   integer(HID_T) fhnd1                 ! file handle                                            
   integer(HID_T) ahnd1                 ! attribute handle                                       
   integer(HID_T) dhnd1                 ! dataset handle                                         
   integer(HID_T) shnd1                 ! dataspace handle                                       
   integer(HSIZE_T) sz1(3)              ! array size 1                                           
   integer(HSIZE_T) sz2(3)              ! array size 2
    
! array data            
   integer(4) l(fov,nscan),c(fov,nscan)
   real(8) x1,y1,sd,sn,s1,s2,s3,sxy
   character(4) year_c,month_c,day_c,hour_c,minute_c,second_c
   real(8) timehdf(2,nscan)
   character(19) trans_time
   integer(8) int_time(2,nscan)
   integer::year=0,month=0,day=0,hour=0,minute=0,second=0 
   real(4) :: lat(fov,nscan)            ! lat                             
   real(4) :: lon(fov,nscan)            ! lon 
!  integer(4) :: lof(fov,nscan)         ! land ocean flag                             
   integer(1) :: clm(fov,nscan)
   real(4)    :: sun_az(fov,nscan)      ! sun_azimuth                                          
   real(4)    :: sat_az(fov,nscan)                                                    
   real(4)    :: sun_zen(fov,nscan)     ! sun_zenith                                            
   real(4)    :: sat_zen(fov,nscan)                                                     
   real(r_kind),parameter  :: tbmin  = 50._r_kind                                                   
   real(r_kind),parameter  :: tbmax  = 550._r_kind

   real(kind=8)                   :: obs_time                                                         
   type (datalink_type),pointer   :: head, p, current, prev
   type(info_type)                :: info                                                                                          
   type(model_loc_type)           :: loc                                   
   integer(i_kind)    :: idate5(6)                                                                   
   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id                                    
   real(r_kind)      :: tb1, crit
   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse
   integer           :: i,j,k,m,n,a,ifile, landsea_mask ,ncid, dimid, varid,ierr
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth                                                           
   integer(i_kind)  :: num_agri_local, num_agri_global, num_agri_used, num_agri_thinned                   
   integer(i_kind)  :: num_agri_used_tmp, num_agri_file
   integer(i_kind)  :: num_agri_local_local, num_agri_global_local, num_agri_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1,filename2,filename3  
   integer          :: nchan,ifov,iscan,ichannels
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max),fname_geo(nfile_max),fname_clm(nfile_max)
   logical          :: fexist

! Allocatable arrays                                           
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)
   real,allocatable             :: obstime(:,:)
   real(4),allocatable          :: sunglint(:,:)
   integer,parameter  :: num_chan=7
   integer(4)         :: nom(fov,nscan,num_chan)
   character(len=80) nom_name(num_chan)
   data nom_name/'NOMChannel08','NOMChannel09','NOMChannel10',&
                 'NOMChannel11','NOMChannel12',&
                 'NOMChannel13','NOMChannel14'/
   real(4)            :: cal(4096,num_chan)
   character(len=80) cal_name(num_chan)
   data cal_name/'CALChannel08','CALChannel09','CALChannel10',&
                 'CALChannel11','CALChannel12',&
                 'CALChannel13','CALChannel14'/
   real(4)            :: tb(fov,nscan,num_chan)             ! tb 
   include "netcdf.inc"
   if (trace_use) call da_trace_entry("da_read_obs_hdf5agri")

!  0.0  Initialize variables                                 
!-----------------------------------
   head_allocated = .false.
   platform_id  = 52  ! rttov user guide Table-2 Col 1 corresponding to 'fy4'         
   satellite_id = 1   ! Table-2 Col 3   use fy-4a id 
   sensor_id    = 99 ! Table-3 Col 2 corresponding to 'fy4-1-agri'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_agri_file    = 0
   num_agri_local   = 0
   num_agri_global  = 0
   num_agri_used    = 0
   num_agri_thinned = 0

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
          (/"The combination of Satellite_Id and Sensor_Id for fy-4 is not found"/))            
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if

! Initialize HDF5 library and Fortran90 interface                                            
   call H5open_f(iret)                                                                       
   if(iret.lt.0)then
      call da_warning(__FILE__,__LINE__, &
           (/"Problems in Initializing HDF5 library. Can not read fy-4 AGRI HDF5 data. "/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if

   nchan = iv%instid(inst)%nchan                                                              
   write(unit=stdout,fmt=*)'fy4 agri nchan: ',nchan                                             
   allocate(data_all(1:nchan))                                                                 

! 1.0 Assign file names and prepare to read fy4 agri files                                        
!-------------------------------------------------------------------------
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   ! first check if agri.hdf is available                                                  
   filename1 = trim(infile_tb)
   inquire (file=filename1, exist=fexist)                                                   
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
   else
      ! check if agri-0x.hdf is available for multiple input files                         
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max                                                                    
         write(filename1,fmt='(A,2I1)') trim(infile_tb),0,i	 
         inquire (file=filename1, exist=fexist)
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
         (/"No valid fy-4 agri.hdf or agri-01.hdf file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if
   
  ! fy4 agri geo.hdf 
   nfile       = 0
   fname_geo = '' !initialize
   ! check if geo.hdf is available 
   filename2 = trim(infile_geo)                                                        
   inquire (file=filename2, exist=fexist)                                                    
   if ( fexist ) then
       nfile=1
       fname_geo  = filename2
   else
      ! check if geo-0x.hdf is available for multiple input files                         
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max 
		 write(filename2,fmt='(A,2I1)') trim(infile_geo),0,i                    
         inquire (file=filename2, exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_geo(nfile)  = filename2
         else
            exit
         end if
      end do
   end if  
	  
   if ( nfile == 0 ) then                                                                
      call da_warning(__FILE__,__LINE__, &
         (/"No valid fy-4 agri geo.hdf or geo-01.hdf file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if
   
  ! clm.nc 
   nfile       = 0
   fname_clm = '' !initialize
   ! check if clm.nc is available 
   filename3 = trim(infile_clm)                                                      
   inquire (file=filename3, exist=fexist)                                                    
   if ( fexist ) then
       nfile=1
       fname_clm  = filename3
   else
      ! check if clm-01.nc is available for multiple input files                         
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max    
         write(filename3,fmt='(A,2I1)') trim(infile_clm),0,i                        
         inquire (file=filename3, exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_clm(nfile)  = filename3
         else
            exit
         end if
      end do
   end if  
	  
   if ( nfile == 0 ) then                                                                
      call da_warning(__FILE__,__LINE__, &
         (/"No valid fy-4 clm.nc or clm-01.nc file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5agri")
      return
   end if

   ! read file: agri and geo and clm for fy4 agri data
  
   infile_loop:  do ifile = 1, nfile                                                    
      num_agri_file_local    = 0
      num_agri_local_local   = 0
      num_agri_global_local  = 0

   ! open HDF5 file for read
      call H5Fopen_f(fname_tb(ifile),H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)                                    
      if(iret.lt.0)then                                                                                            
         call da_warning(__FILE__,__LINE__, &
              (/"Cannot open HDF5 file "//trim(fname_tb(ifile))/))
         cycle infile_loop
      end if
      
   ! read array: scantime
   ! read
      call H5Dopen_f(fhnd1,'NOMObsTime',dhnd1,iret)                                                
        sz1(1)=nscan
        sz1(2)=2
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_DOUBLE,timehdf(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
           if(iret.lt.0)then
              call da_warning(__FILE__,__LINE__, &
               (/"HDF5 read error for: Obs Time"/))
           endif
      call H5Dclose_f(dhnd1,iret)
      allocate  (obstime(1:time_dims,1:nscan))  ! year, month, day, hour, min, sec             
      do iscan=1,nscan
          int_time(1,iscan)=int(timehdf(1,iscan),8)
          if(int_time(1,iscan) .ne. 9999) then
             write(trans_time,"(I17)") int_time(1,iscan)
             year_c=trans_time(1:4)
             month_c=trans_time(5:6)
             day_c=trans_time(7:8)
             hour_c=trans_time(9:10)
             minute_c=trans_time(11:12)
             second_c=trans_time(13:14)
             read(year_c,*) year 
             read(month_c,*) month
             read(day_c,*) day
             read(hour_c,*) hour
             read(minute_c,*) minute
             read(second_c,*) second
             obstime(1,iscan) = year
             obstime(2,iscan) = month
             obstime(3,iscan) = day
             obstime(4,iscan) = hour
             obstime(5,iscan) = minute
             obstime(6,iscan) = second
		   end if
      end do
      write(unit=stdout,fmt=*)'time(scan=1500) year: ',obstime(1,1500),' month:',obstime(2,1500),' day: ',obstime(3,1500),&
         ' hour: ',obstime(4,1500),' minute: ',obstime(5,1500),' second: ',obstime(6,1500)                 
      
   !read brightness temperature               !
   ! read array: nom for gray scale                                    
      do k=1,num_chan
         call H5Dopen_f(fhnd1,nom_name(k),dhnd1,iret)                           
            sz1(1)=nscan
            sz1(2)=fov
         call H5Dread_f(dhnd1, &                                                         
            H5T_NATIVE_INTEGER,nom(:,:,k),sz1,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning(__FILE__,__LINE__, &
                  (/"HDF5 read error for: Gray Scale"/))
         endif
         call H5Dclose_f(dhnd1,iret)                                                   
         if (print_detail_rad) then
             write(unit=message(1),fmt='(A,I6,A,I6)')&
               'nom(pixel=1499,scan=1500,chan=',k,'): ',nom(1500,1499,k)                                     
             call da_message(message(1:1))
         endif
      end do
      
   ! read array: cal for calibration                                  
      do k=1,num_chan
         call H5Dopen_f(fhnd1,cal_name(k),dhnd1,iret)                           
            sz1(1)=4096
         call H5Dread_f(dhnd1, &                                                         
            H5T_NATIVE_REAL,cal(:,k),sz1,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning(__FILE__,__LINE__, &
                  (/"HDF5 read error for: calibration"/))
         endif
         call H5Dclose_f(dhnd1,iret)                                                   
         if (print_detail_rad) then
            write(unit=message(1),fmt='(A,I6,A,F10.4)')&
               'cal(number=500,chan=',k,'): ',cal(500,k)                    
            call da_message(message(1:1))
         endif
      end do

   ! calculate IR brightness temperature 
     do k=1,7
       do j=1,2748
         do i=1,2748
           a=nom(i,j,k)
           if (a.eq.65535) then
             tb(i,j,k)=-999.
           else
             tb(i,j,k)=cal(a+1,k)
           end if
         end do
       end do
     end do
   ! sample play
      write(unit=message(1),fmt='(A,I6,A,F10.4)')&
        'tb(fov=1000,scan=1000,chan=',2,'): ',tb(1000,1000,2)
      call da_message(message(1:1))
   ! close agri file 
      call H5Fclose_f(fhnd1,iret)   
   
   ! open geo HDF5 file for read
      call H5Fopen_f(fname_geo(ifile),H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)                                    
        if(iret.lt.0)then                                                                                            
           call da_warning(__FILE__,__LINE__, &
              (/"Cannot open HDF5 file "//trim(fname_geo(ifile))/))
           cycle infile_loop
        end if  
   ! read line number 
      call h5dopen_f(fhnd1,"LineNumber",dhnd1,iret)
        sz1(1)=nscan
        sz1(2)=fov
      call h5dread_f(dhnd1,H5T_NATIVE_Integer,l(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
           call da_warning(__FILE__,__LINE__, &
                  (/"HDF5 read error for:line number "/))   
        end if
      call h5dclose_f(dhnd1,iret)
   ! read column number  
      call h5dopen_f(fhnd1,"ColumnNumber",dhnd1,iret)
        sz1(1)=nscan
        sz1(2)=fov
      call h5dread_f(dhnd1,H5T_NATIVE_Integer,c(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
           call da_warning(__FILE__,__LINE__, &
                  (/"HDF5 read error for: column number"/)) 
        end if
      call h5dclose_f(dhnd1,iret)
      write(unit=stdout,fmt=*) c(1000,1009),c(999,1009),l(1000,1008),l(999,1008)
   ! calculate lat lon
    do iscan=1,nscan
      do ifov=1,fov
        x1=pai*(c(ifov,iscan)-COFF)/(180*CFAC)*(2**16)
        y1=pai*(l(ifov,iscan)-LOFF)/(180*LFAC)*(2**16)
        sd=((h*cos(x1)*cos(y1))**2-(cos(y1)**2+ea**2/(eb**2)*(sin(y1)**2))*(h**2-ea**2))**(1./2.)
        sn=(h*cos(x1)*cos(y1)-sd)/(cos(y1)**2+ea**2/(eb**2)*sin(y1)**2)
        s1=h-sn*cos(x1)*cos(y1)
        s2=sn*sin(x1)*cos(y1)
        s3=-sn*sin(y1)
        sxy=(s1**2+s2**2)**(1./2.)
        lon(ifov,iscan)=180./pai*atan(s2/s1)+cetlon
        lat(ifov,iscan)=180./pai*atan(ea**2/(eb**2)*(s3/sxy))
      end do
    end do
  ! sample play 
      write(unit=stdout,fmt=*) lat(1373,1374),lon(1373,1374)
      write(unit=stdout,fmt=*) lat(1250,1250),lon(1250,1250)  
  
  ! read array: satellite zenith                                                 
  ! read
      call H5Dopen_f(fhnd1, &
         'NOMSatelliteZenith',dhnd1,iret)
       sz1(1)=nscan
       sz1(2)=fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sat_zen(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
           call da_warning(__FILE__,__LINE__, &
             (/"HDF5 read error for: satellite zenith"/))
        endif
      call H5Dclose_f(dhnd1,iret)
  ! sample play
      write(unit=stdout,fmt=*)'sat_zen(pixel=1000,scan=1000): ',sat_zen(1000,1000)

  ! read array: sun zenith                                          
      call H5Dopen_f(fhnd1, &
         'NOMSunZenith',dhnd1,iret)
        sz1(1)=nscan
        sz1(2)=fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sun_zen(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
          call da_warning(__FILE__,__LINE__, &
            (/"HDF5 read error for: Sun zenith"/))
        endif
      call H5Dclose_f(dhnd1,iret)
   ! sample display
      write(unit=stdout,fmt=*)'sun_zen(pixel=1000,scan=1000): ',sun_zen(1000,1000)
     
   ! read array: Satellite azimuth                                                     
   ! read
      call H5Dopen_f(fhnd1, &
         'NOMSatelliteAzimuth',dhnd1,iret)
        sz1(1)=nscan
        sz1(2)=fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sat_az(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
          call da_warning(__FILE__,__LINE__, &
            (/"HDF5 read error for: NOMSatelliteAzimuth"/))
        endif
      call H5Dclose_f(dhnd1,iret)
   ! sample play     
       write(unit=stdout,fmt=*)'sat_az(pixel=1000,scan=1000): ',sat_az(1000,1000)
     
   ! read array: Sun azimuth                                                     
   ! read
      call H5Dopen_f(fhnd1, &
         'NOMSunAzimuth',dhnd1,iret)
        sz1(1)=nscan
        sz1(2)=fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sun_az(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
        if(iret.lt.0)then
          call da_warning(__FILE__,__LINE__, &
            (/"HDF5 read error for: NOMSunAzimuth"/))
        endif
      call H5Dclose_f(dhnd1,iret)
   ! sample play      
       write(unit=stdout,fmt=*)'sun_az(pixel=1000,scan=1000): ',sun_az(1000,1000)
	
   ! read array: sun glint angle
     allocate(sunglint(fov,nscan))
	call h5dopen_f(fhnd1,"NOMSunGlintAngle",dhnd1,iret)
      sz1(1)=nscan
      sz1(2)=fov
    call h5dread_f(dhnd1,H5T_NATIVE_REAL,sunglint(:,:),sz1,iret,H5S_ALL_F,H5S_ALL_F)
       if(iret.lt.0)then
	   call da_warning(__FILE__,__LINE__, &
            (/"HDF5 read error for:sun glint angle"/))     
       endif
    call h5dclose_f(dhnd1,iret)
   ! sample play      
       write(unit=stdout,fmt=*)'sunglint(pixel=1000,scan=1000): ',sunglint(1000,1000)

   ! close file and HDF5
      call H5Fclose_f(fhnd1,iret)  
	  
   ! open clm.nc file for read
   ! read CLM
      ierr=nf_open(fname_clm(ifile), nf_nowrite, ncid)
      ierr=nf_inq_varid (ncid, 'CLM', varid)
      ierr=nf_get_var (ncid,varid,clm)
         if(ierr.ne.0)then
           call da_warning(__FILE__,__LINE__, &
             (/"nc file read error for: CLM"/))
         end if
   !sample play
      write(unit=stdout,fmt=*) 'clm(scan=1544,pixel=605)',clm(605,1544)
     ierr=nf_close(varid)
     ierr=nf_close(ncid)

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
      scan_loop:     do iscan=1, nscan
	  
         if(int_time(1,iscan) .eq. 9999) cycle scan_loop
         do i = 1, 6
           idate5(i)=obstime(i, iscan)
         end do
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)                       
         if ( obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time) ) cycle scan_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
               obs_time  < time_slots(ifgat) ) exit
         end do

      ! start fov_loop                                   
         fov_loop:      do ifov=1,fov
            num_agri_file       = num_agri_file + 1
            num_agri_file_local = num_agri_file_local + 1
            info%lat  =  lat(ifov,iscan)                              
            info%lon  =  lon(ifov,iscan)
            call da_llxy (info, loc, outside, outside_all)                        
            if (outside_all) cycle fov_loop
            
            num_agri_global       = num_agri_global + 1
            num_agri_global_local = num_agri_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE
            
         ! Discard data over Land (before:landmask =1 -->Land or =0 -->Sea;after:landmask =0 -->Land or =1 -->Sea)                                              !   landsea_mask = 0                                         !!!!!land
         !   if(lof(ifov,iscan) .eq. 0 ) landsea_mask = 1             !!!!!sea
         !   if( landsea_mask == 0 ) cycle fov_loop

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
               dlat_earth = info%lat               !degree
               dlon_earth = info%lon
               if (dlon_earth<zero)  dlon_earth = dlon_earth+r360
               if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
               dlat_earth = dlat_earth*deg2rad             !radian           
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
      
            do k=1,7
               tb1 = tb(ifov,iscan,k)
               if( tb1 < tbmin .or. tb1 > tbmax ) tb1 = missing_r
               data_all(k)= tb1
            end do
! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            p%info             = info            
            p%loc              = loc         
            p%landsea_mask     = 1                
            p%scanpos          = ifov          
            p%satzen           = sat_zen(ifov,iscan)
            p%solzen           = sun_zen(ifov,iscan)
            p%satazi           = sat_az(ifov,iscan)
            p%solazi           = sun_az(ifov,iscan)
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat
            p%cloudflag        = clm(ifov,iscan)
            p%sun_glint        = sunglint(ifov,iscan)
			
            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
         end do fov_loop
      end do scan_loop
   ! Dellocate arrays
      deallocate  (obstime)
      deallocate  (sunglint)

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_file    : ',num_agri_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_global  : ',num_agri_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_agri_local   : ',num_agri_local_local
   end do infile_loop

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

   write(unit=stdout,fmt='(a)') 'fy4 agri data counts: '
   write(stdout,fmt='(a,i10)') ' In file: ',num_agri_file
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
