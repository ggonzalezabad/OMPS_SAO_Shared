SUBROUTINE omi_pge_fitting_process ( pge_idx, n_max_rspec,             &
                                     pge_error_status )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_data_module, ONLY: allocate_variables, deallocate_variables
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, &
       pge_errstat_fatal, error_check, f_sep, omsao_f_subroutine, &
       vb_lev_default, pge_error_status_exit
  USE OMSAO_he5_module, ONLY: NrofScanLines, NrofCrossTrackPixels
  USE OMSAO_variables_module, ONLY: pcfvar
  USE OMSAO_solcomp_module, ONLY: soco_pars_deallocate
  USE OMSAO_omps_reader, ONLY: omps_nmev_type, omps_nmto3_type, &
       omps_nmev_reader, omps_nmto3_reader

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE(omps_nmev_type), ALLOCATABLE :: omps_data
  TYPE(omps_nmto3_type):: omps_to3_data
  INTEGER (KIND=i4) :: nTimesRad, nXtrackRad, nWvlCCD
  INTEGER (KIND=i4) :: errstat
  CHARACTER(LEN=23) :: modulename = 'omi_pge_fitting_process'

  ! ---------------------
  ! OMPS reader variables
  ! ---------------------
  INTEGER (KIND=i4) :: omps_reader_status
  
  pge_error_status = pge_errstat_ok

  ! ----------------------------------------------------------------------------------
  ! Since the OMPS NM files are not that big I'm going to read here the whole file and
  ! assign the values needed to the significant variables.
  ! ----------------------------------------------------------------------------------
  ALLOCATE (omps_data)
  omps_reader_status = OMPS_NMEV_READER(omps_data,TRIM(ADJUSTL(pcfvar%l1b_rad_fname)))
  CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, &
       OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Read OMPS radiance data.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  ! Allocate data variables
  CALL allocate_variables (omps_data%nxtrack,omps_data%nlines,omps_data%nwavel,errstat)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Allocate radiance variables.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  ! Copy and assing values to data variables from omps_data
  CALL omps_to_data_variables ( omps_data, &
       omps_data%nlines, omps_data%nxtrack, omps_data%nwavel)
  
  omps_reader_status = OMPS_NMTO3_READER(omps_to3_data,TRIM(ADJUSTL(pcfvar%l2_to3_fname)))
  CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, &
       OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Read OMPS total ozone data.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  ! ---------------------------------------------
  ! Fudge to fill up the OMI variables to be used
  ! later on (irradiance, radiance, ...)
  ! After this don't need to use other
  ! readers.
  ! Substitutes: omi_read_irradiance
  !              omi_read_radiance
  ! ---------------------------------------------
  NrofScanLines        = omps_data%nLines
  NrofCrossTrackPixels = omps_data%nXtrack
  nTimesRad            = omps_data%nLines
  nXtrackRad           = omps_data%nXtrack
  nWvlCCD              = omps_data%nWavel  

  ! Now I can deallocate omps_data)
  DEALLOCATE(omps_data)

  CALL omi_fitting ( &
       pge_idx, nTimesRad,   nXtrackRad, n_max_rspec,  pge_error_status )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"omi_fitting", vb_lev_default, pge_error_status )

  ! -------------------------------------------------
  ! Deallocation of some potentially allocated memory
  ! -------------------------------------------------
  CALL soco_pars_deallocate (errstat)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"soco_pars_deallocate.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  CALL deallocate_variables (errstat)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"deallocate_variables.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666
  
  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
666 CONTINUE
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  RETURN
END SUBROUTINE omi_pge_fitting_process


SUBROUTINE omi_fitting ( &
       pge_idx, nTimesRad, nXtrackRad, n_max_rspec, pge_error_status )

  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_fatal, &
       pge_errstat_warning, pge_errstat_error, f_sep, omsao_w_nopixel, &
       omsao_w_subroutine, vb_lev_default, error_check, omsao_f_subroutine
  USE OMSAO_indices_module, ONLY: max_rs_idx
  USE OMSAO_variables_module, ONLY: ctrvar, pcfvar, n_rad_wvl
  USE OMSAO_data_module, ONLY: latitude, column_amount, &
       cross_track_skippix, longitude, column_uncert, &
       n_comm_wvl, szenith, fit_rms, vzenith, fitconv_flag, height, &
       ins_database, ins_database_wvl, n_ins_database_wvl
  USE OMSAO_he5_datafields_module
  USE OMSAO_solar_wavcal_module, ONLY: xtrack_solar_calibration_loop
  USE OMSAO_radiance_ref_module, ONLY: &
       xtrack_radiance_reference_loop, create_radiance_reference
  USE OMSAO_wfamf_module, ONLY: read_climatology, amf_calculation_bis
  USE OMSAO_Reference_sector_module, ONLY: Reference_Sector_Correction
  USE OMSAO_omps_reader, ONLY: omps_nmev_type, omps_nmev_reader
  USE OMSAO_database_module, ONLY: xtrack_prepare_database, deallocate_hr_database

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nTimesRad, nXtrackRad, n_max_rspec
  
  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  TYPE(omps_nmev_type), ALLOCATABLE :: omps_data_radiance_reference
  INTEGER (KIND=i4) :: nTimesRadRR, nWvlCCDrr, nXtrackRadRR
  INTEGER (KIND=i4) :: first_line, last_line, errstat, first_wc_pix, last_wc_pix, &
       omps_reader_status
  REAl (KIND=r8), DIMENSION (1) :: zerovec = 0.0_r8
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRad,0:nTimesRad-1) :: saomqf
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: refmqf
  REAL    (KIND=r8), DIMENSION (1:nXtrackRad,0:nTimesRad-1) :: saoamf
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: refamf
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRad,0:nTimesRad-1) :: amfflg
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: refamfflg

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! and Common Mode spectrum
  ! ----------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (0:nTimesRad-1,2) :: omi_xtrpix_range
  LOGICAL,           DIMENSION (0:nTimesRad-1)   :: yn_common_range, yn_radfit_range

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'omi_fitting'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4), EXTERNAL :: &
       he5_init_swath, he5_define_fields, he5_close_output_file, &
       he5_set_field_attributes, he5_write_global_attributes,    &
       he5_write_swath_attributes, he5_open_readwrite

  errstat = pge_errstat_ok

  ! -------------------------------------------------------------------
  ! Range of cross-track pixels to fit. This is based on the selection
  ! in the fitting control file.
  !
  ! NOTE that we set OMI_XTRPIX_RANGE for all swath lines because the
  ! choice of swath lines to process may not contain the radiance 
  ! reference/calibration line.
  ! -------------------------------------------------------------------
  CALL omi_set_xtrpix_range ( &
     nTimesRad, nXtrackRad, ctrvar%pixnum_lim(3:4), &
     omi_xtrpix_range(0:nTimesRad-1,1:2), &
     first_wc_pix, last_wc_pix, errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"omi_set_xtrpix_range", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  ! -----------------------------------------------------
  ! Obtain number of levels in climatology and if present
  ! read values for molecule of intetest (<--FIXME)  
  ! -----------------------------------------------------
  CALL read_climatology ( errstat )
  
  ! ----------------------------------------
  ! Initialization of HE5 output data fields
  ! ----------------------------------------
  errstat = HE5_Init_Swath ( TRIM(ADJUSTL(pcfvar%l2_fname)) )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"HE5_Init_Swath", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  CALL he5_initialize_datafields ( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"HE5_Initialize_datafields", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  errstat = HE5_Define_Fields ( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"HE5_Define_Fields", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  ! ------------------------
  ! Write geolocation fields
  ! ------------------------
  CALL he5_write_geolocation ( nTimesRad, nXtrackRad, errstat)

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------
  CALL xtrack_solar_calibration_loop ( first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! ---------------------------------------
  ! If we are not using radiance reference
  ! this call to prepare database is enough
  ! otherwise we have call it again after
  ! calculating the radiance reference
  ! ---------------------------------------
  CALL xtrack_prepare_database( &
       first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )
  
  ! ------------------------------------------------
  ! If use radiance reference we read the radiance
  ! reference file and work out an efective radiance
  ! reference that is saved to the solar array
  ! ------------------------------------------------
  IF (ctrvar%yn_radiance_reference) THEN

     ALLOCATE (omps_data_radiance_reference)
     omps_reader_status = OMPS_NMEV_READER(omps_data_radiance_reference, &
          TRIM(ADJUSTL(pcfvar%l1b_radref_fname)))
     CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, &
          OMSAO_F_SUBROUTINE, &
          modulename//f_sep//"Read OMPS radiance reference data.", vb_lev_default, pge_error_status )
     IF (pge_error_status >= pge_errstat_error ) GO TO 666

     ! Copy dimensions from radiance reference in to variables
     nTimesRadRR = omps_data_radiance_reference%nLines
     nXtrackRadRR = omps_data_radiance_reference%nXtrack
     nWvlCCDrr = omps_data_radiance_reference%nWavel  

     ! --------------------------------------------------------------------
     ! If the radiance reference is obtained from the same L1b file, we can
     ! simply copy the variables we have just read to the corresponding 
     ! "rr" ones (in this case, the dimensions are the same). Otherwise we
     ! have to read them from the radiance reference granule.
     ! --------------------------------------------------------------------
     CALL create_radiance_reference (omps_data_radiance_reference, nTimesRadRR, nXtrackRadRR, &
          nWvlCCDrr, errstat)
     CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
          modulename//f_sep//"create_radiance_reference", vb_lev_default, &
          pge_error_status )
     IF (pge_error_status >= pge_errstat_error ) GO TO 666

  END IF

  ! --------------------------------------------------------------
  ! Terminate on not having any cross-track pixels left to process
  ! --------------------------------------------------------------
  IF ( ALL ( cross_track_skippix ) ) THEN
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_NOPIXEL, &
          modulename, vb_lev_default, errstat )
     GO TO 666
  END IF

  ! ---------------------------------------------------------------------
  ! If we are using a radiance reference AND want to remove the target
  ! gas from it (important for BrO, for example), we have to run through
  ! all spectra that go into the reference, compute the average column,
  ! and then remove that from the averaged radiance reference spectrum
  ! (owing to the fact that the average of hundreds of OMI spectra still
  !  doesn't produce a decently fitted column).
  ! ---------------------------------------------------------------------
  IF ( ctrvar%yn_radiance_reference) THEN 

     ! --------------------------------------------------
     ! Retrieve target column in radiance reference using
     ! the solar irradiance. If ctrvar%yn_remove_target
     ! then the target column is substracted from
     ! radiance reference.
     ! -------------------------------------------------
     CALL xtrack_radiance_reference_loop ( .TRUE., ctrvar%yn_remove_target, &
          nXtrackRadRR, nWvlCCDrr, first_wc_pix, last_wc_pix, pge_error_status )

     ! -------------------------------------------
     ! Once we have the radiance reference ready
     ! we have to update the database of reference 
     ! cross sections.
     ! -------------------------------------------
     CALL xtrack_prepare_database( &
          first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )

     IF (ctrvar%yn_remove_target) &
          ! -----------------------------------------------------------------
          ! Because we have updated the radiance reference after substracting
          ! the target column we can now retrieve the target column against
          ! the radiance reference stored in the ins_database as a security
          ! check.
          ! -----------------------------------------------------------------
          CALL xtrack_radiance_reference_loop ( .FALSE.,.FALSE., &
          nXtrackRadRR, nWvlCCDrr, first_wc_pix, last_wc_pix, pge_error_status )
     
  END IF

  ! ------------------------------------
  ! Output database of reference spectra
  ! to L2 file
  ! ------------------------------------
  IF( ctrvar%yn_diagnostic_run ) THEN
     n_rad_wvl = MAXVAL(n_ins_database_wvl)
     CALL he5_write_ins_database(ins_database(1:max_rs_idx,1:n_rad_wvl,1:nxtrackrad), &
          ins_database_wvl(1:n_rad_wvl, 1:nxtrackrad), &
          max_rs_idx, n_rad_wvl, nxtrackrad, errstat) 
  ENDIF
  
  ! -----------------------------------------------------------------
  ! Now we enter the on-line computation of the common mode spectrum.
  ! -----------------------------------------------------------------
  IF ( ctrvar%yn_common_iter ) THEN
     
     ! ----------------------------------------------------------
     ! Set the logical YN array that determines which swath lines
     ! will be used in the common mode
     ! ----------------------------------------------------------
     yn_common_range(0:nTimesRad-1) = .FALSE.
     CALL find_swathline_range ( &
          nTimesRad, nXtrackRad, latitude(1:nXtrackRad,0:nTimesRad-1),       &
          ctrvar%common_latrange(1:2), yn_common_range(0:nTimesRad-1), errstat             )

     ! ------------------------------------------------------------------
     ! Interface to the loop over all swath lines
     ! ------------------------------------------
     CALL omi_pge_swathline_loop (                     &
          nTimesRad, nxtrackRad, n_max_rspec, &
          yn_common_range(0:nTimesRad-1),              &
          omi_xtrpix_range(0:nTimesRad-1,1:2),         &
          .FALSE., .TRUE., pge_error_status ) !Logicals for writing out and common mode fit

     ! ---------------------------------------------------
     ! Set the index value of the Common Mode spectrum and
     ! assign values to the fitting parameter arrays
     ! ---------------------------------------------------
     CALL compute_common_mode ( .FALSE., nXtrackRad, 1, zerovec, zerovec, .TRUE. )
     
     ! -------------------------------------------
     ! Write the just computed common mode to file
     ! -------------------------------------------
     IF ( ctrvar%yn_diagnostic_run ) &
          CALL he5_write_common_mode ( nXtrackRad, n_comm_wvl, pge_error_status )

     ! -----------------------------------------------
     ! Write the database with the common mode to file
     ! -----------------------------------------------
     IF( ctrvar%yn_diagnostic_run ) THEN
        n_rad_wvl = MAXVAL(n_ins_database_wvl)
        CALL he5_write_ins_database(ins_database(1:max_rs_idx,1:n_rad_wvl,1:nxtrackrad), &
             ins_database_wvl(1:n_rad_wvl, 1:nxtrackrad), &
             max_rs_idx, n_rad_wvl, nxtrackrad, errstat) 
     END IF
  END IF
  
  ! ----------------------------------------------------------
  ! Now into the proper fitting, with or without common mode.
  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! Set the logical YN array that determines which swath lines
  ! will be processed. Unless we have constrained either swath
  ! line numbers or the latitude range, this can be set to
  ! .TRUE. universically.
  ! ----------------------------------------------------------
  ! First, set the range of swath lines to process
  ! ----------------------------------------------
  first_line = 0  ;  last_line = nTimesRad-1
  IF ( ctrvar%pixnum_lim(1) > 0 ) first_line = MIN(ctrvar%pixnum_lim(1), last_line)
  IF ( ctrvar%pixnum_lim(2) > 0 ) last_line  = MAX( MIN(ctrvar%pixnum_lim(2), last_line), first_line )

  yn_radfit_range = .FALSE.
  IF ( first_line > 0 .OR. &
       last_line < nTimesRad-1 .OR. &
       ctrvar%radfit_latrange(1) > -90.0_r4 .OR. &
       ctrvar%radfit_latrange(2) < +90.0_r4 ) THEN
     IF ( ctrvar%radfit_latrange(1) > -90.0_r4 .OR. &
          ctrvar%radfit_latrange(2) < +90.0_r4 ) THEN
        CALL find_swathline_range ( &
             nTimesRad, nXtrackRad, latitude(1:nXtrackRad,0:nTimesRad-1), &
             ctrvar%radfit_latrange(1:2), yn_radfit_range(0:nTimesRad-1), errstat )
     ELSE
        yn_radfit_range = .TRUE.
        IF ( first_line > 0           ) yn_radfit_range(0:first_line-1)          = .FALSE.
        IF ( last_line  < nTimesRad-1 ) yn_radfit_range(last_line+1:nTimesRad-1) = .FALSE.
     END IF
  ELSE
     yn_radfit_range = .TRUE.
  END IF

  ! ------------------------------------------
  ! Interface to the loop over all swath lines
  ! ------------------------------------------
  CALL omi_pge_swathline_loop (                     &
       nTimesRad, nXtrackRad, n_max_rspec, &
       yn_radfit_range(0:nTimesRad-1),              &
       omi_xtrpix_range(0:nTimesRad-1,1:2),         &
       ctrvar%yn_diagnostic_run, .FALSE., pge_error_status ) ! Logical for commiting and common mode

  IF ( ctrvar%yn_radiance_reference) &
       DEALLOCATE (omps_data_radiance_reference)

  CALL deallocate_hr_database(errstat)
  IF (errstat /= pge_errstat_ok) THEN
     WRITE(*,*) 'Error deallocating high resolution database'
     STOP
  END IF

  ! ---------------------------
  ! SCD to VCD (AMF calculation
  ! ---------------------------
!!$  CALL amf_calculation_bis ( pge_idx, pcfvar%cld_fname, &
!!$       nTimesRad, nXtrackRad,                                &
!!$       latitude(1:nXtrackRad,0:nTimesRad-1),             &
!!$       longitude(1:nXtrackRad,0:nTimesRad-1),            &
!!$       szenith(1:nXtrackRad,0:nTimesRad-1),              &
!!$       vzenith(1:nXtrackRad,0:nTimesRad-1),              &
!!$       omi_xtrpix_range(0:nTimesRad-1,1:2),                  &
!!$       column_amount(1:nXtrackRad,0:nTimesRad-1),        &
!!$       column_uncert(1:nXtrackRad,0:nTimesRad-1),        &
!!$       saoamf(1:nXtrackRad,0:nTimesRad-1),                   &
!!$       amfflg(1:nXtrackRad,0:nTimesRad-1),                   &
!!$       height(1:nXtrackRad,0:nTimesRad-1), .TRUE.,       &
!!$       errstat)       

  ! ---------------------------
  ! Work out fitting statistics
  ! ---------------------------
  CALL compute_fitting_statistics (nTimesRad, nXtrackRad, &
       omi_xtrpix_range(0:nTimesRad-1,1:2), &
       column_amount(1:nXtrackRad,0:nTimesRad-1), &
       column_uncert(1:nXtrackRad,0:nTimesRad-1), &
       fit_rms(1:nXtrackRad,0:nTimesRad-1), &
       fitconv_flag(1:nXtrackRad,0:nTimesRad-1), &
       saomqf(1:nXtrackRad,0:nTimesRad-1), .TRUE., errstat)  

  ! -------------------------------------------------------------
  ! Write to L2 file results:
  !   column_amount(r8) -------> Column Amount(r8)
  !   column_uncertainty(r8) --> Column Uncertainty(r8)
  !   fitconv_flag(i2) --------> Fit Convergence Flag(i2)
  !   fit_rms(r8) -------------> Fitting RMS(r8)
  ! -------------------------------------------------------------
  CALL he5_write_results ( nTimesRad, nXtrackRad, errstat)
  stop
  ! -----------------------------------------------------
  ! Compute reference sector correction:
  ! 0. Save radiance reference data in to omi variables
  ! 1. Perform retrieval over the Ocean for whole granule
  ! 2. AMF calculation over the Ocean
  ! 3. Work out reference sector correction
  ! -----------------------------------------------------
  IF ( ctrvar%yn_refseccor) THEN
     CALL omps_data_to_omi_variables ( OMPS_data_radiance_reference, &
          nTimesRadRR, nXtrackRadRR, nWvlCCDrr)
     
     yn_radfit_range = .TRUE.
     omi_xtrpix_range(0:nTimesRadRR-1,1) = 1
     omi_xtrpix_range(0:nTimesRadRR-1,2) = nXtrackRadRR
     
     ! -----------------
     ! Perform retrieval
     ! -----------------
     CALL omi_pge_swathline_loop (                         &
          nTimesRadRR, nXtrackRadRR, n_max_rspec, &
          yn_radfit_range(0:nTimesRadRR-1),                &
          omi_xtrpix_range(0:nTimesRadRR-1,1:2),           &
          ctrvar%yn_radiance_reference, .FALSE.,           &
          .FALSE., pge_error_status )
     
     ! -----------------------
     ! Perform AMF calculation
     ! ------------------------------------
     ! Change name of the file used for the 
     ! clouds to reference sector one
     ! ------------------------------------
     
     CALL amf_calculation_bis ( pge_idx, pcfvar%refsec_cld_fname, &
          nTimesRadRR, nXtrackRadRR,                                   &
          latitude(1:nXtrackRadRR,0:nTimesRadRR-1),                &
          longitude(1:nXtrackRadRR,0:nTimesRadRR-1),               &
          szenith(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
          vzenith(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
          omi_xtrpix_range(0:nTimesRadRR-1,1:2),                       &
          column_amount(1:nXtrackRadRR,0:nTimesRadRR-1),           &
          column_uncert(1:nXtrackRadRR,0:nTimesRadRR-1),           &
          refamf(1:nXtrackRadRR,0:nTimesRadRR-1),                      &
          refamfflg(1:nXtrackRadRR,0:nTimesRadRR-1),                   &
          height(1:nXtrackRadRR,0:nTimesRadRR-1), .FALSE.,         &
          errstat)       

     ! ---------------------------
     ! Work out fitting statistics
     ! ---------------------------
     CALL compute_fitting_statistics (pge_idx, nTimesRadRR, nXtrackRadRR, &
          omi_xtrpix_range(0:nTimesRadRR-1,1:2), &
          column_amount(1:nXtrackRadRR,0:nTimesRadRR-1), &
          column_uncert(1:nXtrackRadRR,0:nTimesRadRR-1), &
          fit_rms(1:nXtrackRadRR,0:nTimesRadRR-1), &
          fitconv_flag(1:nXtrackRadRR,0:nTimesRadRR-1), &
          refmqf(1:nXtrackRadRR,0:nTimesRadRR-1), .FALSE., errstat)

     CALL Reference_sector_Correction ( nTimesRad, nXtrackRad, nTimesRadRR, nXtrackRadRR, &
          refmqf(1:nXtrackRadRR,0:nTimesRadRR-1), refamf(1:nXtrackRadRR,0:nTimesRadRR-1), &
          refamfflg(1:nXtrackRadRR,0:nTimesRadRR-1), pge_error_status)

  ENDIF

  ! ---------------------
  ! Write some attributes
  ! ---------------------
  errstat = pge_errstat_ok
  errstat = he5_write_global_attributes( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_WRITE_GLOBAL_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_write_swath_attributes ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_WRITE_SWATH_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_set_field_attributes   ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_SET_FIELD_ATTRIBUTES", vb_lev_default, pge_error_status )

  ! -----------------
  ! Close output file
  ! -----------------
  errstat = pge_errstat_ok
  errstat = he5_close_output_file ( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_CLOSE_OUTPUT_FILE", vb_lev_default, pge_error_status )

  RETURN
  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
666 CONTINUE
  IF ( pge_error_status >= pge_errstat_error ) RETURN

END SUBROUTINE omi_fitting
