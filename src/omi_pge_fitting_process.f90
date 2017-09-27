SUBROUTINE omi_pge_fitting_process ( pge_idx, n_max_rspec,             &
                                     pge_error_status )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_omidata_module, ONLY: allocate_radiance_variables
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
  TYPE(omps_nmev_type) :: omps_data, omps_data_radiance_reference
  TYPE(omps_nmto3_type):: omps_to3_data
  INTEGER (KIND=i4) :: nTimesRad,   nXtrackRad,   nWvlCCD
  INTEGER (KIND=i4) :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr
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
  omps_reader_status = OMPS_NMEV_READER(omps_data,TRIM(ADJUSTL(pcfvar%l1b_rad_fname)))
  CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Read OMPS radiance data.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  ! Allocate program variables
  CALL allocate_radiance_variables (omps_data%nxtrack,omps_data%nlines,omps_data%nwavel,errstat)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Allocate radiance variables.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  stop
  omps_reader_status = OMPS_NMEV_READER(omps_data_radiance_reference,TRIM(ADJUSTL(pcfvar%l1b_radref_fname)))
  CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"Read OMPS radiance reference data.", vb_lev_default, pge_error_status )
  IF (pge_error_status >= pge_errstat_error ) GO TO 666

  omps_reader_status = OMPS_NMTO3_READER(omps_to3_data,TRIM(ADJUSTL(pcfvar%l2_to3_fname)))
  CALL error_check ( INT(omps_reader_status,KIND=i4), pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
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
  nTimesRadRR          = omps_data_radiance_reference%nLines
  nXtrackRad           = omps_data%nXtrack
  nXtrackRadRR         = omps_data_radiance_reference%nXtrack
  nWvlCCD              = omps_data%nWavel  
  nWvlCCDrr            = omps_data_radiance_reference%nWavel  

  CALL omps_data_to_omi_variables ( omps_data, &
       nTimesRad, nXtrackRad, nWvlCCD)

  stop

  CALL omi_fitting (                                  &
       pge_idx, omps_data_radiance_reference,         &
       nTimesRad,   nXtrackRad, n_max_rspec, &
       nTimesRadRR, nXtrackRadRR, nWvlCCDrr, pge_error_status       )

  IF ( pge_error_status >= pge_errstat_fatal ) GO TO 666

  ! -------------------------------------------------
  ! Deallocation of some potentially allocated memory
  ! -------------------------------------------------
  CALL soco_pars_deallocate (errstat)  
  
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


SUBROUTINE omi_fitting (                                  &
       pge_idx, omps_data_radiance_reference, &
       nTimesRad,   nXtrackRad, n_max_rspec, &
       nTimesRadRR, nXtrackRadRR, nWvlCCDrr, &
       pge_error_status )

  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_fatal, &
       pge_errstat_warning, pge_errstat_error, f_sep, omsao_w_nopixel, &
       omsao_w_subroutine, vb_lev_default, error_check
  USE OMSAO_indices_module, ONLY: sao_molecule_names
  USE OMSAO_parameters_module, ONLY: i2_missval
  USE OMSAO_variables_module, ONLY: ctrvar, pcfvar
  USE OMSAO_omidata_module, ONLY: latitute, column_amount, &
       cross_track_skippix, radcal_itnum, radcal_xflag, &
       solcal_itnum, solcal_xflag, &
       longitude, column_uncert, n_comm_wvl, szenith, &
       fit_rms, vzenith, fitconv_flag, height
  USE OMSAO_he5_module, ONLY:  pge_swath_name
  USE OMSAO_he5_datafields_module
  USE OMSAO_solar_wavcal_module, ONLY: xtrack_solar_calibration_loop
  USE OMSAO_radiance_ref_module, ONLY: &
       radiance_reference_lnums, xtrack_radiance_reference_loop, create_radiance_reference
  USE OMSAO_wfamf_module, ONLY: omi_read_climatology, CmETA, amf_calculation_bis
  USE OMSAO_pixelcorner_module, ONLY: compute_pixel_corners
  USE OMSAO_Reference_sector_module, ONLY: Reference_Sector_Correction
  USE OMSAO_omps_reader, ONLY: omps_nmev_type

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  TYPE(omps_nmev_type), INTENT(IN) :: omps_data_radiance_reference
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nTimesRad, nXtrackRad, n_max_rspec
  INTEGER (KIND=i4), INTENT (IN) :: nTimesRadRR, nWvlCCDrr, nXtrackRadRR
  
  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  INTEGER   (KIND=i4) ::                                                 &
       iline, first_line, last_line, errstat, first_wc_pix, last_wc_pix, &
       first_pix, last_pix
  REAl (KIND=r8), DIMENSION (1) :: zerovec = 0.0_r8
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRad,0:nTimesRad-1)     :: saomqf
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRadRR,0:nTimesRadRR-1) :: refmqf
  REAL    (KIND=r8), DIMENSION (1:nXtrackRad,0:nTimesRad-1)     :: saoamf
  REAL    (KIND=r8), DIMENSION (1:nXtrackRadRR,0:nTimesRadRR-1) :: refamf
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRadRR,0:nTimesRadRR-1) :: amfflg
  INTEGER (KIND=i2), DIMENSION (1:nXtrackRadRR,0:nTimesRadRR-1) :: refamfflg

  ! ----------------------------------------------------------------------
  ! Swath dimensions and variables that aren't passed from calling routine
  ! ----------------------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: molname

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! and Common Mode spectrum
  ! ----------------------------------------------------------
  INTEGER (KIND=i1), DIMENSION (0:nTimesRad-1)   :: omi_binfac
  INTEGER (KIND=i4), DIMENSION (0:nTimesRad-1,2) :: omi_xtrpix_range
  LOGICAL,           DIMENSION (0:nTimesRad-1)   :: &
       omi_yn_szoom, yn_common_range, yn_radfit_range

  INTEGER (KIND=i1), DIMENSION (0:nTimesRadRR-1)   :: omi_binfac_rr
  INTEGER (KIND=i4), DIMENSION (0:nTimesRadRR-1,2) :: omi_xtrpix_range_rr
  LOGICAL,           DIMENSION (0:nTimesRadRR-1)   :: omi_yn_szoom_rr

  ! ----------------------------------------------------------
  ! OMI L1b latitudes
  ! ----------------------------------------------------------
  REAL (KIND=r4), DIMENSION (1:nXtrackRad, 0:nTimesRad-1) :: l1b_latitudes

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

  ! ---------------------------------------------------------------
  ! Some initializations that will save us headaches in cases where 
  ! a proper set-up of those variables failes or is bypassed.
  ! ---------------------------------------------------------------
  first_pix = 1 ; last_pix = 1

  errstat = pge_errstat_ok

  ! --------------------------------
  ! Name of the main output molecule
  ! --------------------------------
  molname = sao_molecule_names(pge_idx)

  ! ------------------------------------
  ! For OMPS we are never in a zoom mode
  ! ------------------------------------
  omi_yn_szoom = .FALSE.

  ! -------------------------------------------------------------------
  ! Range of cross-track pixels to fit. This is based on the selection
  ! in the fitting control file and whether the granule being processed
  ! is in global or spatial zoom mode, or even a mixture thereof.
  !
  ! NOTE that we set OMI_XTRPIX_RANGE for all swath lines because the
  ! choice of swath lines to process may not contain the radiance 
  ! reference/calibration line.
  ! -------------------------------------------------------------------
  CALL omi_set_xtrpix_range ( &
     nTimesRad, nXtrackRad, ctrvar%pixnum_lim(3:4), &
     omi_xtrpix_range(0:nTimesRad-1,1:2), &
     first_wc_pix, last_wc_pix, errstat )

  ! --------------------------------------------------------------------
  ! If the radiance reference is obtained from the same L1b file, we can
  ! simply copy the variables we have just read to the corresponding 
  ! "rr" ones (in this case, the dimensions are the same). Otherwise we
  ! have to read them from the radiance reference granule.
  ! --------------------------------------------------------------------
  IF ( TRIM(ADJUSTL(pcfvar%l1b_radref_fname)) /= TRIM(ADJUSTL(pcfvar%l1b_rad_fname)) ) THEN
    CALL omi_set_xtrpix_range ( &
          nTimesRadRR, nXtrackRadRR, ctrvar%pixnum_lim(3:4), &
           omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2), &
          first_wc_pix, last_wc_pix, errstat )
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ELSE
     omi_binfac_rr      (0:nTimesRad-1)       = omi_binfac      (0:nTimesRad-1)
     omi_yn_szoom_rr    (0:nTimesRadRR-1)     = omi_yn_szoom    (0:nTimesRad-1)
     omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2) = omi_xtrpix_range(0:nTimesRad-1,1:2)
  END IF

  ! ---------------------------------------------------------------
  ! The Climatology is going to be read here and kept in memory. If
  ! this has a bad impact in the efficiency of the application then
  ! I will find a different way. We are doing this to be able to in
  ! itialize the output he5 with the correct number of levels for
  ! the Scattering weights and Gas_profiles output. We are going to
  ! use the number of levels in the climatology as the number of le
  ! vels of the reported scattering weights.
  ! ---------------------------------------------------------------
  CALL omi_read_climatology ( errstat )
  
  ! ----------------------------------------
  ! Initialization of HE5 output data fields
  ! ----------------------------------------
  errstat = HE5_Init_Swath ( TRIM(ADJUSTL(pcfvar%l2_fname)), pge_swath_name, nTimesRad, nXtrackRad, CmETA )
  CALL he5_initialize_datafields ( )

  errstat = HE5_Define_Fields ( pge_idx, pge_swath_name, nTimesRad, nXtrackRad, CmETA )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  
  ! ------------------------
  ! Write geolocation fields
  ! ------------------------
  CALL he5_write_geolocation ( nTimesRad, nXtrackRad, &
       first_wc_pix, last_wc_pix, errstat)

  ! ---------------------------------------------------------------
  ! Work out pixel corners so we can plot and write them to L2 file
  ! ---------------------------------------------------------------
  CALL compute_pixel_corners (nTimesRad, nXtrackRad, &
       latitute(1:nXtrackRad, 0:nTimesRad-1),    &
       longitude(1:nXtrackRad, 0:nTimesRad-1),   &
       errstat)

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------
  solcal_itnum = i2_missval ; solcal_xflag = i2_missval
  CALL xtrack_solar_calibration_loop ( first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! ------------------------------------------------
  ! If use radiance reference we read the radiance
  ! reference file and work out an efective radiance
  ! reference that is saved to the solar array
  ! ------------------------------------------------
  IF (ctrvar%yn_radiance_reference) THEN
     CALL create_radiance_reference (nTimesRadRR, nXtrackRadRR, nWvlCCDrr, errstat)
  END IF

  ! -----------------------------------------------------
  ! Across-track loop for radiance wavelength calibration
  ! -----------------------------------------------------
  radcal_itnum = i2_missval ; radcal_xflag = i2_missval
  CALL xtrack_radiance_wvl_calibration (                          &
       ctrvar%yn_radiance_reference, ctrvar%yn_solar_comp,                      &
       first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

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

     iline = SUM ( radiance_reference_lnums(1:2) ) / 2
     IF ( iline < 0 .OR. iline > nTimesRadRR ) iline = nTimesRadRR / 2
     first_pix = omi_xtrpix_range_rr(iline,1)
     last_pix  = omi_xtrpix_range_rr(iline,2)
     CALL xtrack_radiance_reference_loop (                                   &
          ctrvar%yn_radiance_reference, ctrvar%yn_remove_target,                           &
          nXtrackRadRR, nWvlCCDrr, first_pix, last_pix, pge_idx, pge_error_status )

     ! -------------------------------------------------------------
     ! Write the output from solar/earthshine wavelength calibration
     ! and radiance reference to file. The latter results will be
     ! overwritten in the call to XTRACK_RADIANCE_REFERENCE_LOOP
     ! below, hence we need to write them out here.
     ! -------------------------------------------------------------
     CALL he5_write_wavcal_output ( nXtrackRad, first_pix, last_pix, errstat )

  END IF

  ! -----------------------------------------------------------------
  ! Before we go any further we need to read the L1b latitude values,
  ! since we base our screening of which swath lines to process on
  ! those values. Both common mode, if used, and the radiance fit
  ! uses the same arrays, so we read this only ones.
  !
  ! We could shave off some fractional minute from the run time by
  ! not reading the latitudes in cases where no radiance reference
  ! is used, i.e., where both radiance granule and radiance reference
  ! granule are the same. The down-side is an increase in virtual
  ! memory program uses, plus some more logic to find out whether to
  ! read the latitudes or not. For now we are going with a second
  ! read, particularly since the current algorithm settings would
  ! require it anyway.
  ! -----------------------------------------------------------------
  ! For OMPS it is enough to copy the values of latitute in
  ! l1b_latitudes
  ! -----------------------------------------------------------
  l1b_latitudes(1:nXtrackRad,0:nTimesRad-1) = latitute(1:nXtrackRad,0:nTimesRad-1)

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
          nTimesRad, nXtrackRad, l1b_latitudes(1:nXtrackRad,0:nTimesRad-1),       &
          ctrvar%common_latrange(1:2), yn_common_range(0:nTimesRad-1), errstat             )

     ! -------------------------------------------------------------
     ! First and last swath line number will be overwritten. Hence
     ! we save them for further reference.
     ! -------------------------------------------------------------
     !first_line_save = first_line
     !last_line_save  = last_line

     ! ------------------------------------------------------------------
     ! There is a time saver catch built into the assignment below, 
     ! which we probably want to rethink. If we are only doing a few
     ! lines but the common mode extends over a wide range of a
     ! latitudes, then we would be processing a lot of lines to derive
     ! the common mode. This is currently being excluded, but the better
     ! way may be to remember that we can control everything through
     ! the fitting control file.
     ! ------------------------------------------------------------------
     ! Interface to the loop over all swath lines
     ! ------------------------------------------
     CALL omi_pge_swathline_loop (                     &
          pge_idx, nTimesRad, nxtrackRad, n_max_rspec, &
          yn_common_range(0:nTimesRad-1),              &
          omi_xtrpix_range(0:nTimesRad-1,1:2),         &
          ctrvar%yn_radiance_reference, .FALSE.,       &
          .FALSE., pge_error_status )
 
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
  IF ( first_line         > 0           .OR. &
       last_line          < nTimesRad-1 .OR. &
       ctrvar%radfit_latrange(1) > -90.0_r4    .OR. &
       ctrvar%radfit_latrange(2) < +90.0_r4           ) THEN

     IF ( ctrvar%radfit_latrange(1) > -90.0_r4    .OR. &
          ctrvar%radfit_latrange(2) < +90.0_r4           ) THEN
        CALL find_swathline_range ( &
             nTimesRad, nXtrackRad, l1b_latitudes(1:nXtrackRad,0:nTimesRad-1),       &
             ctrvar%radfit_latrange(1:2), yn_radfit_range(0:nTimesRad-1), errstat             )
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
       pge_idx, nTimesRad, nXtrackRad, n_max_rspec, &
       yn_radfit_range(0:nTimesRad-1),              &
       omi_xtrpix_range(0:nTimesRad-1,1:2),         &
       ctrvar%yn_radiance_reference, .FALSE.,       &
       .TRUE., pge_error_status )

  ! ---------------------------
  ! SCD to VCD (AMF calculation
  ! ---------------------------
  CALL amf_calculation_bis ( pge_idx, pcfvar%cld_fname, &
       nTimesRad, nXtrackRad,                                &
       latitute(1:nXtrackRad,0:nTimesRad-1),             &
       longitude(1:nXtrackRad,0:nTimesRad-1),            &
       szenith(1:nXtrackRad,0:nTimesRad-1),              &
       vzenith(1:nXtrackRad,0:nTimesRad-1),              &
       omi_xtrpix_range(0:nTimesRad-1,1:2),                  &
       column_amount(1:nXtrackRad,0:nTimesRad-1),        &
       column_uncert(1:nXtrackRad,0:nTimesRad-1),        &
       saoamf(1:nXtrackRad,0:nTimesRad-1),                   &
       amfflg(1:nXtrackRad,0:nTimesRad-1),                   &
       height(1:nXtrackRad,0:nTimesRad-1), .TRUE.,       &
       errstat)       

  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
666 CONTINUE
  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  ! ---------------------------
  ! Work out fitting statistics
  ! ---------------------------
  CALL compute_fitting_statistics (pge_idx, nTimesRad, nXtrackRad, &
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
          pge_idx, nTimesRadRR, nXtrackRadRR, n_max_rspec, &
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
          latitute(1:nXtrackRadRR,0:nTimesRadRR-1),                &
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
END SUBROUTINE omi_fitting
