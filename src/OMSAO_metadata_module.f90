MODULE OMSAO_metadata_module

  ! ==============================
  ! Module for Metadata parameters
  ! ==============================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: &
       maxchlen, str_missval, i4_missval, blank23, blank25, blank27, blank30
  USE OMSAO_indices_module,    ONLY: n_sao_pge, sao_pge_min_idx, sao_pge_max_idx

  IMPLICIT NONE

  INCLUDE "PGS_MET.f"
  INCLUDE "PGS_PC.f"

  ! NULL element value
  CHARACTER(LEN=4), PARAMETER :: null_element='NULL'

  ! Define metadata type
  TYPE met_type
     CHARACTER(LEN=maxchlen) :: element, element_description, element_value
     LOGICAL :: yn_dynamic
  END type met_type

  ! Variable holding number of elements
  INTEGER(KIND=4) :: number_of_metadata_elements
  ! Number of metadata elements
  TYPE(met_type), ALLOCATABLE :: metadata_struct(:)

  ! String for date and time data
  INTEGER(KIND=4), PARAMETER :: dt_len = 25
  CHARACTER(LEN=dt_len) :: dt_str

  ! Define metadata elements. If yn_dynamic element_value will be filled up somewhere else

  ! -----------------------
  ! MetaData STRING Objects
  ! ----------------------------------------------------------------------
  ! NOTE: "l1r" refers to the processed earthshine granule, while
  !       "l1R" to the possibly used earthshine radiance reference granule
  ! ----------------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_str = 24
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_str), PARAMETER :: &
       mdata_string_fields = RESHAPE ( (/                 &
       "ShortName                        ", "inv"//blank30, "mcf"//blank30, &
       "AssociatedSensorShortName.1      ", "inv"//blank30, "mcf"//blank30, &
       "AssociatedPlatformShortName.1    ", "inv"//blank30, "mcf"//blank30, &
       "AssociatedInstrumentShortName.1  ", "inv"//blank30, "mcf"//blank30, &
       "LongName                         ", "arc"//blank30, "mcf"//blank30, &
       "ESDTDescriptorRevision           ", "arc"//blank30, "mcf"//blank30, &
       "EquatorCrossingDate.1            ", "inv"//blank30, "l1r"//blank30, &
       "EquatorCrossingTime.1            ", "inv"//blank30, "l1r"//blank30, &
       "PGEVERSION                       ", "inv"//blank30, "l1i"//blank30, &
       "PGEVERSION                       ", "inv"//blank30, "l1r"//blank30, &
       "PGEVERSION                       ", "inv"//blank30, "l1R"//blank30, &
       "RangeBeginningDate               ", "inv"//blank30, "l1r"//blank30, &
       "RangeEndingDate                  ", "inv"//blank30, "l1r"//blank30, &
       "RangeBeginningTime               ", "inv"//blank30, "l1r"//blank30, &
       "RangeEndingTime                  ", "inv"//blank30, "l1r"//blank30, &
       "ShortName                        ", "inv"//blank30, "l1i"//blank30, &
       "ShortName                        ", "inv"//blank30, "l1r"//blank30, &
       "ShortName                        ", "inv"//blank30, "l1R"//blank30, &
       "AssociatedSensorShortName.1      ", "inv"//blank30, "l1r"//blank30, &
       "AssociatedPlatformShortName.1    ", "inv"//blank30, "l1r"//blank30, &
       "AssociatedInstrumentShortName.1  ", "inv"//blank30, "l1r"//blank30, &
       "OrbitData                        ", "arc"//blank30, "l1r"//blank30, &
       "PGEVERSION                       ", "inv"//blank30, "pge"//blank30, &
       "ShortName                        ", "inv"//blank30, "pge"//blank30 /), &
       (/ 3, n_mdata_str /) )
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_str) :: mdata_string_values

  ! --------------------------------------------------------------------------
  ! MetaData STRING Objects specific to OMHCHO: Pre-fitted slant columns of
  ! Bro ("BRO") and O3 ("OOO").
  ! --------------------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_omhcho = 4
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_omhcho), PARAMETER :: &
       mdata_omhcho_fields = RESHAPE ( (/                 &
       "PGEVERSION                       ", "inv"//blank30, "BRO"//blank30,    &
       "PGEVERSION                       ", "inv"//blank30, "OOO"//blank30,    &
       "ShortName                        ", "inv"//blank30, "BRO"//blank30,    &
       "ShortName                        ", "inv"//blank30, "OOO"//blank30 /), &
       (/ 3, n_mdata_omhcho /) )
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_omhcho) :: mdata_omhcho_values

  INTEGER   (KIND=i4), PARAMETER :: n_mdata_omchocho = 2
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_omchocho), PARAMETER :: &
       mdata_omchocho_fields = RESHAPE ( (/                 &
       "PGEVERSION                       ", "inv"//blank30, "LQH2O "//blank27,    &
       "ShortName                        ", "inv"//blank30, "LQH2O "//blank27 /), &
       (/ 3, n_mdata_omchocho /) )
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_omchocho) :: mdata_omchocho_values


  ! --------------------------------------------------------------------------
  ! MetaData STRING Objects specific (and common) to OMHCHO and OMCHOCHO: 
  ! cloud parameters from a cloud PGE ("CLD").
  ! --------------------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_voc = 2
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_voc), PARAMETER :: &
       mdata_voc_fields = RESHAPE ( (/                 &
       "PGEVERSION                       ", "inv"//blank30, "CLD"//blank30, &
       "ShortName                        ", "inv"//blank30, "CLD"//blank30 /), &
       (/ 3, n_mdata_voc /) )
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_mdata_voc) :: mdata_voc_values

  ! ---------------------------------------------------------------------------
  ! MetaData STRING Product Specific Attributes (STRING)
  ! ---------------------------------------------------------------------------
  ! We keep these separate from the other MetaData fields because the PSAs will
  ! be copied between .met files only, but this involves jumping through some
  ! more hoops than just the use of the Toolkit routines. They also don't have
  ! a VALUES arrays associtated with them, since we read and write on the fly.
  ! ---------------------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_psa = 14
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_psa), PARAMETER :: &
       mdata_psa_fields = RESHAPE ( (/                                   &
       "EndBlockNr                    ", "psa"//blank27, "l1r"//blank27, &
       "ExpeditedData                 ", "psa"//blank27, "l1r"//blank27, &
       "ExposureTimes                 ", "psa"//blank27, "l1r"//blank27, &
       "InstrumentConfigurationIDs    ", "psa"//blank27, "l1r"//blank27, &
       "MasterClockPeriods            ", "psa"//blank27, "l1r"//blank27, &
       "NrMeasurements                ", "psa"//blank27, "l1r"//blank27, &
       "NrSpatialZoom                 ", "psa"//blank27, "l1r"//blank27, &
       "NrSpectralZoom                ", "psa"//blank27, "l1r"//blank27, &
       "NrZoom                        ", "psa"//blank27, "l1r"//blank27, &
       "PathNr                        ", "psa"//blank27, "l1r"//blank27, &
       "SolarEclipse                  ", "psa"//blank27, "l1r"//blank27, &
       "SouthAtlanticAnomalyCrossing  ", "psa"//blank27, "l1r"//blank27, &
       "SpacecraftManeuverFlag        ", "psa"//blank27, "l1r"//blank27, &
       "StartBlockNr                  ", "psa"//blank27, "l1r"//blank27  /), &
       (/ 3, n_mdata_psa /) )

  ! ----------------------------------------------
  !!! Removed because OMOCLO is retrieved from VIS 
  ! ----------------------------------------------
  !     "QAStatPctMeasErrorUV2         ", "psa"//blank27, "l1r"//blank27, &
  !     "QAStatPctPixBadUV2            ", "psa"//blank27, "l1r"//blank27, &
  !     "QAStatPctPixProcessingErrorUV2", "psa"//blank27, "l1r"//blank27, &

  ! ------------------------
  ! MetaData INTEGER Objects
  ! ------------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_int = 2
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_int), PARAMETER :: &
       mdata_integer_fields = RESHAPE ( (/                 &
       "VersionID                 ", "inv"//blank23, "pcf"//blank23, &
       "OrbitNumber.1             ", "inv"//blank23, "l1r"//blank23  /), &
       (/ 3, n_mdata_int /) )

  INTEGER (KIND=i4), DIMENSION(n_mdata_int) :: mdata_integer_values


  ! -----------------------
  ! MetaData DOUBLE Objects
  ! -----------------------
  INTEGER   (KIND=i4), PARAMETER :: n_mdata_dbl = 3
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (3, n_mdata_dbl), PARAMETER :: &
       mdata_double_fields = RESHAPE ( (/                 &
       "EquatorCrossingLongitude.1", "inv"//blank23, "l1r"//blank23,    &
       !"SpaceCraftMeanAltitude    ", "arc"//blank23, "pge"//blank23,    &
       "SpaceCraftMinAltitude     ", "arc"//blank23, "l1r"//blank23,    &
       "SpaceCraftMaxAltitude     ", "arc"//blank23, "l1r"//blank23 /), &
       (/ 3, n_mdata_dbl /) )

  REAL (KIND=r8), DIMENSION(n_mdata_dbl) :: mdata_double_values


  ! ---------------------------------------
  ! Base names for SAO L2 output parameters
  ! ------------------------------------------------------------------
  ! These Base Name are appended with the PGE target molecules during
  ! the initialization of the L2 MetaData parameter values. They will
  ! be accessed as pge_parameter_names(n_para, pge_idx), and so have
  ! to correspond to the order of PGE indices defined in
  ! OMSAO_indices_module:
  !
  !        (1) OClO  (2) BrO  (3) HCHO  (4) O3
  !
  ! ------------------------------------------------------------------
  INTEGER   (KIND=I4),  PARAMETER :: n_l2_output_paras = 1
  CHARACTER (LEN=37), DIMENSION (n_l2_output_paras, sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       pge_parameter_names = RESHAPE ( (/ &
       "Slant Column Chlorine Dioxide        ",   &
       "Total Column Bromine Oxide           ",   &
       "Total Column Formaldehyde            ",   &
       "Slant Column Ozone                   ",   &
       "Slant Column Nitrogen Dioxide        ",   &
       "Slant Column Sulfur Dioxide          ",   &
       "Slant Column Glyoxal                 ",   &
       "Slant Column Iodine Monoxide         ",   &
       "Slant Column Water Vapor             ",   &
       "Slant Column Nitrous Acid            ",   &
       "Slant Column Oxygen Collision Complex",   &
       "Slant Column Liquid Water Absorption ",   &
       "Slant Column Nitrogen Dioxide        " /),&
       (/ n_l2_output_paras, sao_pge_max_idx-sao_pge_min_idx+1 /) )

  ! ---------------------------------------------------
  ! PGE L2 Inventory MetaData fields/values
  ! --------------------------------------------------
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L) :: &
       AutomaticQualityFlag            = str_missval, &
       AutomaticQualityFlagExplanation = str_missval
  INTEGER   (KIND=i4)             :: &
       QAPercentMissingData            = i4_missval
  INTEGER   (KIND=i4)             :: &
       QAPercentOutofBoundsData        = i4_missval

  ! -------------------------------------------------------------
  ! Arrays and variables to hold Metatdata values read from files
  ! -------------------------------------------------------------
  INTEGER (KIND=i4) :: mcf_versionid

  ! -----------------------------------------------------------------------
  ! Some fixed strings. They will be used as short hand for special task,
  ! like extracting the date in the required format from the Metadata field
  ! 'RageBeginningDate', or are recurring strings, like 'CoreMetadata.xxx'
  ! that we don't want to spell out all the time they appear.
  ! -----------------------------------------------------------------------
  CHARACTER (LEN=15), PARAMETER :: amd_str  = "ArchiveMetadata"
  CHARACTER (LEN=12), PARAMETER :: cmd_str  = "CoreMetadata"
  CHARACTER (LEN=18), PARAMETER :: rbd_str  = "RangeBeginningDate"
  CHARACTER (LEN=18), PARAMETER :: rbt_str  = "RangeBeginningTime"
  CHARACTER (LEN=27), PARAMETER :: apsn_str = "AssociatedPlatformShortName"
  CHARACTER (LEN=29), PARAMETER :: aisn_str = "AssociatedInstrumentShortName"
  CHARACTER (LEN= 9), PARAMETER :: vid_str  = "VersionID"
  CHARACTER (LEN= 9), PARAMETER :: sn_str   = "ShortName"
  CHARACTER (LEN=14), PARAMETER :: lgid_str = "LocalGranuleID"
  CHARACTER (LEN=13), PARAMETER :: par_str  = "ParameterName"
  CHARACTER (LEN=11), PARAMETER :: orbn_str = "OrbitNumber"

  ! -------------------------------------------------------
  ! Parameters and variables for L2 Metadata initialization
  ! -------------------------------------------------------
  CHARACTER (LEN=PGSd_MET_GROUP_NAME_L), DIMENSION (PGSd_MET_NUM_OF_GROUPS) :: mcf_groups
  CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX)                                     :: l2_local_granule_id

  ! -------------------------------------------------
  ! STRING variables for various purposes:
  !  * ShortName from MCF
  !  * Granule Start/End Time from PCF
  !    (will be compared against RangeBeginningTime
  !     and RangeEndTime)
  ! -------------------------------------------------
  CHARACTER (LEN=maxchlen) :: mcf_shortname
  CHARACTER (LEN=maxchlen) :: pcf_granule_s_time, pcf_granule_e_time

CONTAINS
  SUBROUTINE read_metadata_table ()

    USE OMSAO_variables_module, ONLY: pcfvar
    IMPLICIT NONE

    ! Local variables
    INTEGER(KIND=4), PARAMETER :: funit = 37, nfld=4
    CHARACTER(LEN=1), PARAMETER :: hdr_chr = '#', sep_chr = ','
    CHARACTER(LEN=maxchlen) :: tmp_line
    INTEGER(KIND=4) :: num_lines, num_hdr, iline, ichr, ifld
    INTEGER(KIND=4), DIMENSION(4) :: fchr, lchr

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
       fchr=1;lchr=1;ifld=1
       DO ichr = 1, maxchlen
          IF (tmp_line(ichr:ichr) .EQ. sep_chr) THEN
             lchr(ifld) = ichr-1
             ifld=ifld+1
             fchr(ifld) = ichr+1
             IF (ifld .EQ. nfld) THEN
                lchr(ifld) = maxchlen
                CONTINUE
             END IF
          END IF
       END DO
       metadata_struct(iline-num_hdr)%element=tmp_line(fchr(1):lchr(1))
       metadata_struct(iline-num_hdr)%element_value=tmp_line(fchr(2):lchr(2))
       IF (tmp_line(fchr(3):lchr(3)) .EQ. 'static') metadata_struct(iline-num_hdr)%yn_dynamic=.FALSE.
       IF (tmp_line(fchr(3):lchr(3)) .EQ. 'dynamic') metadata_struct(iline-num_hdr)%yn_dynamic=.TRUE.
       metadata_struct(iline-num_hdr)%element_description=tmp_line(fchr(4):lchr(4))
    END DO
    CLOSE(funit)

  END SUBROUTINE read_metadata_table

  SUBROUTINE check_metadata()

    IMPLICIT NONE
    CHARACTER(LEN=14), PARAMETER :: location = 'check_metadata'
    INTEGER(KIND=4) :: imet

    CALL generate_date_and_time()
    DO imet = 1, number_of_metadata_elements

       IF (metadata_struct(imet)%yn_dynamic) THEN
          IF (TRIM(metadata_struct(imet)%element_value) .EQ. 'NULL') &
          write(*,10) dt_str, location, TRIM(metadata_struct(imet)%element)
       ENDIF

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
    CHARACTER(len=maxchlen) :: sublocation
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

END MODULE OMSAO_metadata_module
