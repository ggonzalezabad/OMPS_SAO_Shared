MODULE OMSAO_metadata_module

  ! ==============================
  ! Module for Metadata parameters
  ! ==============================
  IMPLICIT NONE

  ! NULL element value
  CHARACTER(LEN=4), PARAMETER :: null_element='NULL'

  ! Define metadata type
  INTEGER(KIND=4), PARAMETER :: maxlen = 1300
  TYPE met_type
     CHARACTER(LEN=maxlen) :: element, element_value
     CHARACTER(LEN=3) :: origin
  END type met_type

  ! Variable holding number of elements
  INTEGER(KIND=4) :: number_of_metadata_elements
  ! Number of metadata elements
  TYPE(met_type), ALLOCATABLE :: metadata_struct(:)

  ! String for date and time data
  INTEGER(KIND=4), PARAMETER :: dt_len = 25
  CHARACTER(LEN=dt_len) :: dt_str


CONTAINS
  SUBROUTINE read_metadata_table ()

    USE OMSAO_variables_module, ONLY: pcfvar
    IMPLICIT NONE

    ! Local variables
    INTEGER(KIND=4), PARAMETER :: funit = 37
    CHARACTER(LEN=1), PARAMETER :: hdr_chr = '#', sep_chr = ','
    CHARACTER(LEN=maxlen) :: tmp_line
    INTEGER(KIND=4) :: num_lines, num_hdr, iline, ichr
    INTEGER(KIND=4) :: fchr, lchr

    ! Create date and time string
    CALL generate_date_and_time()

    ! Get number of elements
    num_lines=0; num_hdr=0
    OPEN(funit, file=TRIM(pcfvar%mtd_fname), status='old', action='read')
    DO
       READ(funit,'(A)',END=10) tmp_line
       IF (tmp_line(1:1) .EQ. hdr_chr) num_hdr = num_hdr + 1
       num_lines = num_lines + 1 
    END DO
    10 REWIND(funit)
    ! Fill up initial metadata_struct values    
    number_of_metadata_elements = num_lines-num_hdr
    ALLOCATE(metadata_struct(1:number_of_metadata_elements))
    DO iline = 1, num_lines
       READ(funit,'(A)') tmp_line
       IF (iline .LE. num_hdr) CYCLE
       ! Split tmp_line in different fields (comma separated) and save them
       ! to metadata_struct
       fchr=1;lchr=1
       DO ichr = 1, maxlen
          IF (tmp_line(ichr:ichr) .EQ. sep_chr) THEN
             lchr = ichr-1
             CONTINUE
          END IF
       END DO
       metadata_struct(iline-num_hdr)%element=tmp_line(fchr:lchr)
       metadata_struct(iline-num_hdr)%element_value=null_element
       fchr=lchr+2; lchr=maxlen
       metadata_struct(iline-num_hdr)%origin=TRIM(ADJUSTL(tmp_line(fchr:lchr)))
    END DO
    CLOSE(funit)
    
  END SUBROUTINE read_metadata_table

  SUBROUTINE check_metadata()

    IMPLICIT NONE
    CHARACTER(LEN=14), PARAMETER :: location = 'check_metadata'
    INTEGER(KIND=4) :: imet

    DO imet = 1, number_of_metadata_elements
       IF (TRIM(metadata_struct(imet)%element_value) .EQ. 'NULL') &
            write(*,10) dt_str, location, TRIM(metadata_struct(imet)%element)
    END DO
10  FORMAT(A,'---',A,': ',A,' not allocated at runtime!!!')

  END SUBROUTINE check_metadata

  SUBROUTINE write_metadata()

    USE HDF5
    USE OMSAO_he5_module, ONLY: pge_swath_file_id, HE5_SWclose
    USE OMSAO_variables_module, ONLY: pcfvar
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=14), PARAMETER :: location = 'write_metadata'
    CHARACTER(len=maxlen) :: sublocation
    INTEGER(KIND=4) :: hdferr, fid, grpid, imet, dspace_id, dset_id, attr_id, atype_id

    ! Close L2 outfile
    hdferr = HE5_SWclose (pge_swath_file_id, hdferr)

    ! Open L2 file using HDF5 interface
    CALL H5OPEN_f(hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5OPEN_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF
    CALL H5FOPEN_f(TRIM(ADJUSTL(pcfvar%l2_fname)), H5F_ACC_RDWR_F, fid, hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5FOPEN_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    ! Open HDFEOS INFORMATION group
    CALL H5GOPEN_f(fid,'HDFEOS INFORMATION',grpid,hdferr)    
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5GOPEN_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    ! Create GranuleMetadata dataset
    CALL h5screate_simple_f(1,(/INT(1,KIND=8)/),dspace_id,hdferr)
    CALL H5DCREATE_f(grpid,'GranuleMetadata',H5T_NATIVE_CHARACTER,dspace_id,dset_id,hdferr)    
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5DCREATE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    
    ! Write metadata to GranuleMetadata group
    DO imet = 1, number_of_metadata_elements
       CALL H5TCOPY_F(H5T_NATIVE_CHARACTER,atype_id,hdferr)
       IF (hdferr .NE. 0) THEN
          sublocation = 'H5TCOPY_f'
          write(*,10) dt_str,location,TRIM(sublocation)
       END IF
       CALL H5TSET_SIZE_F(atype_id,LEN(TRIM(metadata_struct(imet)%element_value),KIND=8), &
            hdferr)
       IF (hdferr .NE. 0) THEN
          sublocation = 'H5TSET_SIZE_f'
          write(*,10) dt_str,location,TRIM(sublocation)
       END IF
       CALL H5ACREATE_f(dset_id, TRIM(metadata_struct(imet)%element), &
            atype_id,dspace_id,attr_id,hdferr)
       IF (hdferr .NE. 0) THEN
          sublocation = 'H5ACREATE_f'
          write(*,10) dt_str,location,TRIM(sublocation)
       END IF
       CALL H5AWRITE_F(attr_id,atype_id, &
            TRIM(metadata_struct(imet)%element_value), (/INT(1,KIND=8)/),hdferr)
       IF (hdferr .NE. 0) THEN
          sublocation = 'H5AWRITE_f'
          write(*,10) dt_str,location,TRIM(sublocation)
       END IF
    END DO

    ! Close GranuleMetadata dateset
    CALL H5DCLOSE_F(dset_id,hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5DCLOSE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    ! Close GranuleMetadata dateset
    CALL H5SCLOSE_F(dspace_id,hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5SCLOSE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    ! Close HDFEOS INFORMATION group
    CALL H5GCLOSE_F(grpid,hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5GCLOSE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

    ! Close L2 file
    CALL H5FCLOSE_f(fid,hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5FCLOSE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF
    
    ! Close hdf interface
    CALL H5CLOSE_f(hdferr)
    IF (hdferr .NE. 0) THEN
       sublocation = 'H5CLOSE_f'
       write(*,10) dt_str,location,TRIM(sublocation)
    END IF

10  FORMAT(A,'---',A,': ',A,' error !!!')
    
  END SUBROUTINE write_metadata

  SUBROUTINE deallocate_metadata()

    IMPLICIT NONE
    IF (ALLOCATED(metadata_struct)) DEALLOCATE(metadata_struct)

  END SUBROUTINE deallocate_metadata

  SUBROUTINE generate_date_and_time ()

    IMPLICIT NONE
    
    !Local variables
    INTEGER(KIND=4), DIMENSION(8) :: values
    INTEGER(KIND=4) :: hour, min

    ! Get local time values
    CALL date_and_time(values=values)

    ! Get hours form UTC
    hour = INT(values(4)/60,KIND=4)
    min  = MOD(values(4),60)

    ! Define output format
    IF (hour .LT. 0 .OR. min .LT. 0) THEN
       WRITE(dt_str,10) values(1),values(2),values(3),values(5),values(6),values(7),ABS(hour),ABS(min)
    ELSE
       WRITE(dt_str,20) values(1),values(2),values(3),values(5),values(6),values(7),hour,min
    END IF

10  FORMAT(I0.4,'-',I0.2,'-',I0.2,'T',I0.2,':',I0.2,':',I0.2,'+',I0.2,':',I0.2)
20  FORMAT(I0.4,'-',I0.2,'-',I0.2,'T',I0.2,':',I0.2,':',I0.2,'-',I0.2,':',I0.2)


  END SUBROUTINE generate_date_and_time

  SUBROUTINE fill_l1b_metadata_values( omps_data )

    USE OMSAO_OMPS_READER, ONLY: omps_nmev_type
    IMPLICIT NONE

    TYPE(omps_nmev_type), INTENT(IN) :: omps_data
    INTEGER(KIND=4) :: imet, il1b

    DO imet = 1, number_of_metadata_elements
       IF (TRIM(ADJUSTL(metadata_struct(imet)%origin)) .NE. 'L1B') CYCLE
       DO il1b = 1, omps_data%nattr
          IF ( TRIM(ADJUSTL(metadata_struct(imet)%element)) .EQ. &
               TRIM(ADJUSTL(omps_data%global_attributes(il1b))) ) THEN
             metadata_struct(imet)%element_value = &
                  TRIM(ADJUSTL(omps_data%global_attributes_values(il1b)))
          END IF
       END DO
    END DO

  END SUBROUTINE fill_l1b_metadata_values

  SUBROUTINE fill_run_metadata_values(omps_data)

    USE OMSAO_OMPS_READER, ONLY: omps_nmev_type
    IMPLICIT NONE

    TYPE(omps_nmev_type), INTENT(IN) :: omps_data
    INTEGER(KIND=4) :: imet

    DO imet = 1, number_of_metadata_elements
       IF (TRIM(ADJUSTL(metadata_struct(imet)%origin)) .NE. 'RUN') CYCLE
       SELECT CASE (TRIM(ADJUSTL(metadata_struct(imet)%element)))
       CASE('ProductionDateTime')
          metadata_struct(imet)%element_value=dt_str
       CASE('date_created')
          metadata_struct(imet)%element_value=dt_str
       CASE('WestBoundingCoordinate')
          WRITE(metadata_struct(imet)%element_value,10) &
               MINVAL(omps_data%LongitudeCorner)
       CASE('NorthBoundingCoordinate')
          WRITE(metadata_struct(imet)%element_value,10) &
               MAXVAL(omps_data%LatitudeCorner)
       CASE('EastBoundingCoordinate')
          WRITE(metadata_struct(imet)%element_value,10) &
               MAXVAL(omps_data%LongitudeCorner)
       CASE('SouthBoundingCoordinate')
          WRITE(metadata_struct(imet)%element_value,10) &
               MINVAL(omps_data%LatitudeCorner)
       END SELECT
    END DO

10 FORMAT(F7.2)

  END SUBROUTINE fill_run_metadata_values

  SUBROUTINE fill_pcf_metadata_values ()

    USE OMSAO_indices_module, ONLY: &
         config_lun_values, config_lun_strings, n_config_luns

    IMPLICIT NONE
    INTEGER(KIND=4) :: imet, icon

    DO imet = 1, number_of_metadata_elements
       IF (TRIM(ADJUSTL(metadata_struct(imet)%origin)) .NE. 'PCF') CYCLE
       DO icon = 1, n_config_luns
          IF ( TRIM(ADJUSTL(metadata_struct(imet)%element)) .EQ. &
               TRIM(ADJUSTL(config_lun_strings(icon))) ) THEN
             metadata_struct(imet)%element_value = TRIM(ADJUSTL(config_lun_values(icon)))
             CONTINUE
          END IF
       END DO
       SELECT CASE ( TRIM(ADJUSTL(metadata_struct(imet)%element)) )
       CASE ('input_granules')
          CALL get_input_granules (imet)
       CASE ('input_ancilliary_data')
       END SELECT
    END DO
  END SUBROUTINE fill_pcf_metadata_values

  SUBROUTINE get_input_granules (imet)

    USE OMSAO_variables_module, ONLY: pcfvar, ctrvar
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN) :: imet
    INTEGER(KIND=4), PARAMETER :: n_lun_inp=6
    CHARACTER(LEN=5), DIMENSION(n_lun_inp), PARAMETER :: &
         input_str = (/'L1bRa','L1bRe','L1bIr','L2Cld','L2O3T','L2Pre'/)
    CHARACTER(LEN=maxlen), DIMENSION(n_lun_inp) :: input_val
    CHARACTER(LEN=maxlen) :: out_string
    INTEGER(KIND=4) :: i, ppos

    ! L1bRa (processed granule)
    input_val(1) = TRIM(ADJUSTL(pcfvar%l1b_rad_fname))
    ! L1bRe (radiance reference granule)
    input_val(2) = TRIM(ADJUSTL(pcfvar%l1b_radref_fname))
    ! L1bIr (solar irradiance)
    input_val(3) = TRIM(ADJUSTL(pcfvar%l1b_irrad_fname))
    ! L2Cld (L2 cloud product)
    input_val(4) = TRIM(ADJUSTL(pcfvar%cld_fname))
    ! L2OT3 (L2 ozone product)
    input_val(5) = TRIM(ADJUSTL(pcfvar%l2_to3_fname))
    ! L2Pre (Prefit granule)
    input_val(6) = TRIM(ADJUSTL(pcfvar%prefit_fname))

    ! Generate final string
    ppos = SCAN(input_val(1),'/', BACK=.true.)
    if (ppos > 0) THEN
       out_string = input_str(1)//': '//TRIM(input_val(1)(ppos+1:maxlen))
    else
       out_string = input_str(1)//': '//TRIM(input_val(1))
    end if
    DO i = 2, n_lun_inp
       IF ( (i .EQ. 2) .AND. (.NOT. ctrvar%yn_radiance_reference) ) CYCLE
       IF ( (i .EQ. 3) .AND. (ctrvar%yn_solar_comp .OR. ctrvar%yn_solmonthave) ) CYCLE
       IF ( (i .EQ. 6) .AND. (.NOT. ctrvar%yn_prefit(1)) ) CYCLE
       ppos = SCAN(input_val(i),'/', BACK=.true.)
       if (ppos > 0) THEN
          out_string = TRIM(out_string)//'; '//input_str(i)//': '//TRIM(input_val(i)(ppos+1:maxlen))
       else
          out_string = TRIM(out_string)//'; '//input_str(i)//': '//TRIM(input_val(i))
       end if
    END DO
    metadata_struct(imet)%element_value = TRIM(out_string)

  END SUBROUTINE get_input_granules

END MODULE OMSAO_metadata_module
