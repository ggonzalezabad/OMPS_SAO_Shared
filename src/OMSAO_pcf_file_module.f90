MODULE OMSAO_pcf_file_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_pcf_file

  CONTAINS

SUBROUTINE read_pcf_file ( pge_error_status )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_indices_module, ONLY: config_lun_array, config_lun_values, &
       config_lun_strings, pge_static_input_luns, pge_l2_output_lun, &
       solmonthave_lun, solcomp_lun, refsec_lun, refsec_cld_lun, &
       albedo_lun, slitfunc_lun, n_config_luns, max_rs_idx, &
       l1b_radiance_lun, l1b_radianceref_lun, l1b_irradiance_lun, &
       icf_idx, prefit_lun, pge_molid_lun, &
       versionid_lun, swathname_lun, instrument_name_lun, &
       pge_version_lun, proclevel_lun, granule_e_lun, granule_s_lun, &
       orbitnumber_lun, verbosity_lun, amf_table_lun, cld_lun, cld_climatology_lun
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, f_sep, &
       omsao_f_get_molindex, omsao_f_getlun, omsao_s_get_molindex, omsao_w_getlun, &
       omsao_w_subroutine, pge_errstat_error, pge_errstat_fatal, pge_errstat_warning, &
       pgs_smf_mask_lev_s, vb_lev_default, vb_lev_stmdebug, PGSd_PC_VALUE_LENGTH_MAX, &
       error_check
  USE OMSAO_metadata_module, ONLY: pcf_granule_s_time,  pcf_granule_e_time, &
       n_mdata_int, mdata_integer_fields, mdata_integer_values
  USE OMSAO_parameters_module, ONLY: zerospec_string, str_missval, maxchlen
  USE OMSAO_he5_module, ONLY: pge_swath_name, process_level, &
       instrument_name, pge_version
  USE OMSAO_variables_module, ONLY: refspecs_original, pcfvar
  USE OMSAO_wfamf_module, ONLY: climatology_lun
  USE OMSAO_control_file_module, ONLY: read_fitting_control_file

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'read_pcf_file'

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                      :: i, j, strlen
  CHARACTER (LEN=maxchlen)                 :: lunstr
  CHARACTER (LEN=PGSd_PC_VALUE_LENGTH_MAX) :: tmpchar

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER (KIND=i4) :: errstat, version  

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=maxchlen), EXTERNAL :: int2string
  INTEGER   (KIND=i4),      EXTERNAL :: &
       pgs_pc_getnumberoffiles, pgs_pc_getreference, pgs_pc_getconfigdata, &
       pgs_smf_teststatuslevel

  errstat = pge_errstat_ok

  ! --------------------------------------------------------------------
  ! Read all ConfigData from PCF and set variables associated with them. 
  ! Some, like MoleculeID, might lead to PGE termination if in error.
  ! --------------------------------------------------------------------
  DO i = 1, n_config_luns
     errstat = &
          PGS_PC_GetConfigData ( config_lun_array(i), config_lun_values(i) )
     errstat = PGS_SMF_TestStatusLevel(errstat)
     CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), vb_lev_default, &
          pge_error_status )
     ! -------------------------------------------------------------
     ! Remove any quotes (") that might have made in into the string
     ! -------------------------------------------------------------
     CALL remove_quotes ( config_lun_values(i) )

     SELECT CASE ( config_lun_array(i) )
        ! --------------------------------------------------------------------
        ! Get the verbosity threshold. The setting of this variable determines
        ! how "chatty" in terms of screen output the PGE will be. If the READ
        ! was in error, proceed with DEBUG level
        ! --------------------------------------------------------------------
     CASE ( verbosity_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "1"
        READ (config_lun_values(i), '(I1)') pcfvar%verb_thresh_lev
        ! ----------------------------------------------------------------------
        ! Get the orbit number. This will be checked against L1B MetaData later.
        ! Note that the source for OrbitNumber is the L1B file, NOT the PCF. But
        ! this only matters in case they are different, which will be checked
        ! later. If the ARE different, then the L1B value will take precidence.
        ! ----------------------------------------------------------------------
     CASE ( orbitnumber_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "00000"
        READ (config_lun_values(i), '(I5)') pcfvar%orbit_number

        ! -----------------------------------------------------------------
        ! Get Granule Start and End time. This is provided by the Scheduler
        ! and has to be checked against the L1B MetaData fields for
        ! RangeBeginningTime and RangeEndingTime.
        ! -----------------------------------------------------------------
     CASE ( granule_s_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_s_time = TRIM(ADJUSTL(config_lun_values(i)))
     CASE ( granule_e_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_e_time = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Process Level
        ! ------------------------------------------------------------------
     CASE ( proclevel_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        process_level = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the PGE Version
        ! ------------------------------------------------------------------
     CASE ( pge_version_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_version = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Instrument Name
        ! ------------------------------------------------------------------
     CASE ( instrument_name_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        instrument_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the HE5 Swath Name; set to MISSING VALUE if it can't be found.
        ! ------------------------------------------------------------------
     CASE ( swathname_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_swath_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ----------------------------------------------------------------------
        ! Get ECS Collection version number and assign it to a Metadata variable
        ! ----------------------------------------------------------------------
     CASE ( versionid_lun )
        READ (config_lun_values(i), '(I1)') pcfvar%ecs_version_id
        getidx: DO j = 1, n_mdata_int
           IF ( mdata_integer_fields(3,j) == "pcf"              .AND. ( &
                INDEX(mdata_integer_fields(1,j), 'VERSIONID') /= 0 .OR. &
                INDEX(mdata_integer_fields(1,j), 'VersionID') /= 0     )  ) THEN
              mdata_integer_values(j) = pcfvar%ecs_version_id
              EXIT getidx
           END IF
        END DO getidx

        ! -------------------------------------------------------------------------
        ! Get the SAO PGE Name string. This string is converted to the PGE
        ! index number (10, 11, 12), which in turn is used to access LUNs and other
        ! elements that are PGE specific. All of this has to be done here, because
        ! we require the PGE index to identify LUNs and other PGE specific items.
        ! >>> We MUST know this, or else we can't execute the PGE <<<
        ! -------------------------------------------------------------------------
     CASE ( pge_molid_lun )
        CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, &
             OMSAO_F_GETLUN, modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), &
             vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_error ) RETURN

        ! ---------------------------------------------------------------
        ! Translate molecule string into index (required for LUN look-up)
        ! ---------------------------------------------------------------
        errstat = pge_errstat_ok
        CALL get_pge_ident ( &
             TRIM(ADJUSTL(config_lun_values(i))), pcfvar%pge_name, pcfvar%pge_idx, errstat )
        CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, &
             OMSAO_F_GET_MOLINDEX, modulename, vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_fatal ) RETURN
        CALL error_check ( 0, 1, pge_errstat_ok, OMSAO_S_GET_MOLINDEX, &
             TRIM(ADJUSTL(pcfvar%pge_name)), vb_lev_stmdebug, errstat )  
     END SELECT
  END DO
  errstat = pge_errstat_ok

  ! ---------------------------------------------------------
  ! Static input file with tabulated OMI slit function values
  ! ---------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (slitfunc_lun, version, pcfvar%slitfunc_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"SLITFUNC_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ------------------------------------------------------------------
  ! Read names of static input files from PCF. The first entry is the 
  ! algorithm control file (initial fitting values, etc.), the other
  ! entries are reference spectra.
  ! ------------------------------------------------------------------
  pcfvar%static_input_fnames = zerospec_string
  DO i = icf_idx, max_rs_idx
     version = 1
     errstat = PGS_PC_GetReference (pge_static_input_luns(i), version, tmpchar)
     tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))

     errstat = PGS_SMF_TestStatusLevel(errstat)
     IF ( errstat /= pgs_smf_mask_lev_s .OR. strlen == 0 ) THEN
        lunstr = int2string ( pge_static_input_luns(i), 1 )
        CALL error_check ( 0, 1, pge_errstat_fatal, OMSAO_F_GETLUN, &
             modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
             vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_error ) RETURN
     ELSE IF ( INDEX ( TRIM(ADJUSTL(tmpchar)), zerospec_string ) == 0 ) THEN
        pcfvar%static_input_fnames(i) = TRIM(ADJUSTL(tmpchar))
     END IF
  END DO
  ! ------------------------------------------------------
  ! Save file names with reference spectra to global array
  ! ------------------------------------------------------
  refspecs_original(1:max_rs_idx)%FileName = pcfvar%static_input_fnames(1:max_rs_idx)

  ! ----------------------------------------------------
  ! Read fitting conrol parameters from input file
  ! ----------------------------------------------------
  errstat = pge_errstat_ok
  CALL read_fitting_control_file ( errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"READ_FITTING_CONTROL_FILE.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  ! -----------------------------
  ! Read Irradiance L1B file name
  ! -----------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_irradiance_lun, version, pcfvar%l1b_irrad_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_IRRADIANCE_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ---------------------------
  ! Read Radiance L1B file name
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_radiance_lun, version, pcfvar%l1b_rad_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_RADIANCE_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -------------------------------------
  ! Read Radiance Reference L1B file name
  ! -------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_radianceref_lun, version, pcfvar%l1b_radref_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_RADIANCEREF_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ------------------------------------------
  ! Read name of file with the radiance clouds
  ! ------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( cld_lun, version, pcfvar%cld_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"CLD_FILENAME_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------------------
  ! Read name of AMF table file. Remember that a missing AMF table is
  ! not a fatal problem, since in that case the slant columns will be written
  ! to the output file.
  ! -------------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( amf_table_lun, version, pcfvar%amf_table_fname)
  errstat = PGS_SMF_TestStatusLevel ( errstat )
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"AMF_TABLE_FILENAME_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -------------------------------
  ! Read cloud climatology filename
  ! -------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( cld_climatology_lun, version, pcfvar%cld_climatology_fname )
  errstat = PGS_SMF_TestStatusLevel ( errstat )
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"CLIMATOLOGY_FILENAME_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  ! -------------------------
  ! Read climatology filename
  ! -------------------------
  version = 1
  errstat = PGS_PC_GetReference ( climatology_lun, version, pcfvar%climatology_fname )
  errstat = PGS_SMF_TestStatusLevel ( errstat )
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"CLIMATOLOGY_FILENAME_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------------------
  ! Read name of file with Solar Spectrum Composite 
  ! (whether we use it or not, since we find that out only after we read the
  !  fitting control file)
  ! -------------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( solcomp_lun, version, pcfvar%solcomp_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"SOLCOMP_FILENAME_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------------------
  ! Read name of file with Monthly Average Irradiace
  ! (whether we use it or not, since we find that out only after we read the
  !  fitting control file) !gga
  ! -------------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( solmonthave_lun, version, pcfvar%solmonthave_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"SOLMONTHAVE_FILENAME_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ------------------------------------------------------------
  ! Read name of file with GEOS-Chem background Reference Sector
  ! concentrations !gga 
  ! ------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( refsec_lun, version, pcfvar%refsec_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
          modulename//f_sep//"REFSEC_FILENAME_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ----------------------------------------------------
  ! Read name of file with the radiance reference clouds
  ! ----------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( refsec_cld_lun, version, pcfvar%refsec_cld_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"REFSEC_CLD_FILENAME_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ----------------------------------
  ! Read name of OMLER albedo file gga
  ! ----------------------------------
  version = 1
  errstat = PGS_PC_GetReference (albedo_lun, version, pcfvar%albedo_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
          modulename//f_sep//"ALBEDO_filename", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ---------------------------
  ! Read name of L2 output file
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_l2_output_lun, version, pcfvar%l2_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"PGE_L2_OUTPUT_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ---------------------------------------
  ! Read name of cloud climatology filename
  ! ---------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_l2_output_lun, version, pcfvar%l2_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"PGE_L2_OUTPUT_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! ---------------------------
  ! Read prefit column filename
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (prefit_lun, version, pcfvar%prefit_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"PREFIT_LUN ", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  RETURN
END SUBROUTINE read_pcf_file

END MODULE OMSAO_pcf_file_module
