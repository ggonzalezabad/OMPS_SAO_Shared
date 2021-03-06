MODULE OMSAO_data_module

  USE OMSAO_precision_module, ONLY: i1, i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: maxchlen, max_spec_pts
  USE OMSAO_indices_module, ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx
  USE OMSAO_variables_module, ONLY: ctrvar, database

  IMPLICIT NONE

  ! -----------------------------
  ! Maximum data/swath dimensions
  ! -----------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       nxtrack_max    =  60, nwavel_max = 300, &
       nwavelcoef_max =   5, nlines_max = 500, &
       nUTCdim        =   27

  ! ---------------------------------------------------
  ! Maximum dimension of data fields in SAO L2 products
  ! ---------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_field_maxdim = 3

  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: NISE_snowfree =   0, NISE_allsnow = 100, NISE_permice = 101, &
       NISE_drysnow = 103, NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! Values to go into the diagnostic array that shows how the AMF was
  ! computed, with values that indicate missing cloud products, glint,
  ! and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: cfr_addmiss = 1000, ctp_addmiss = 2000, glint_add = 10000, &
       geo_amf = -1, oobview_amf = -2, wfmod_amf = -9, bigsza_amf = -3

  ! ------------------------------------------------------------------
  ! Arrays for L1b data (radiance, radiance reference, and irradiance)
  ! ------------------------------------------------------------------
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:) :: spacecraft_alt, spacecraft_alt_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: time, time_reference       
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: radiance_errstat, instrument_flag, &
       radiance_errstat_reference, instrument_flag_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: geoflg, xtrflg, geoflg_reference, xtrflg_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: height, land_water_flg, snowicefraction, &
       height_reference, land_water_flg_reference, snowicefraction_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: latitude, longitude, latitude_reference, longitude_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: szenith, sazimuth, szenith_reference, sazimuth_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: vzenith, vazimuth, vzenith_reference, vazimuth_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: razimuth, azimuth_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:,:) :: latitudecorner, longitudecorner
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_spec, radiance_spec_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_prec, radiance_prec_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_wavl, radiance_wavl_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_qflg, radiance_qflg_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: time_utc, time_utc_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: irradiance_prec, irradiance_wavl, irradiance_spec, &
       irradiance_wght, radref_spec, radref_wavl, radref_wght
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yn_process_pixel
  CHARACTER(LEN=nUTCdim), ALLOCATABLE, DIMENSION(:) :: utc_time

  ! -------------------------------------------------------------------
  ! Arrays for irradiance or radiance reference and ccd pixel selection
  ! -------------------------------------------------------------------
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: irradiance_qflg, radref_qflg
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:) :: radref_sza, radref_vza
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: ccdpix_selection
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: ccdpix_exclusion
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:,:) :: ccdpix_selection_rad
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:,:) :: ccdpix_exclusion_rad

  ! ----------------------------------
  ! Arrays for common mode calculation
  ! ----------------------------------
  INTEGER (KIND=i4) :: n_comm_wvl
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: common_cnt
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: common_spc, common_wvl

  ! -----------------------------------
  ! Arrays for spectral fitting results
  ! -----------------------------------
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: column_amount, column_uncert, fit_rms, radfit_chisq
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: fitconv_flag, itnum_flag

  ! -------------------------------------------------------
  ! Albedo array (to be put together with other AMF arrays)
  ! -------------------------------------------------------
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: albedo

  ! ----------------------------------------------------------------------------
  ! Correlations with main output product. Due to a bug in the HDF-EOS5 routines
  ! (non-TLCF implementation), STRING fields cannot be written to file directly.
  ! A work-around solution is to convert the CHARACTERs to INTEGERs. Thus the
  ! need for the additional array CORRELATION_NAMES_INT.
  ! ----------------------------------------------------------------------------
  CHARACTER (LEN=maxchlen), DIMENSION (n_max_fitpars) :: correlation_names
  CHARACTER (LEN=n_max_fitpars*maxchlen)              :: correlation_names_concat

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  INTEGER (KIND=i4) :: nxtrack_rad, ntime_rad, ntimes_loop
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: nwav_irrad, nwav_radref
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: nwav_rad, nwav_rad_reference

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: n_ins_database_wvl
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:) :: solcal_itnum, radcal_itnum, radref_itnum, &
       solcal_xflag, radcal_xflag, radref_xflag
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: solcal_pars, radcal_pars, radref_pars
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: ins_database
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: ins_database_wvl
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: ins_sol_wav_avg
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: radcal_chisq, radref_chisq, &
       radref_col, radref_dcol, radref_rms, radref_xtrcol

  ! ------------------------------
  ! Distance between Earth and Sun
  ! ------------------------------
  REAL (KIND=r4) :: EarthSunDistance

  ! ---------------------------------------------------------
  ! Scan line, block line, and across-track pixel numbers
  ! ---------------------------------------------------------
  INTEGER (KIND=i4) :: scanline_no, blockline_no

  ! -----------------------
  ! L2 output data QA flags
  ! ----------------------------------------------------------------
  ! Per Centages of output column data that are used to classify the
  ! scientific data quality:
  !         Good data >= QA_PERCENT_PASS   : "Passed"
  !         Good data >= QA_PERCENT_SUSPECT: "Suspect"
  !         Good data <  QA_PERCENT_SUSPECT: "Failed"
  ! ----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: qa_percent_passed = 75, qa_percent_suspect = 50


  ! ------------------------------------------------------
  ! Finally some variables that will be initialized in the
  ! course of the processing.
  ! ------------------------------------------------------
  INTEGER (KIND=i4) :: n_radwvl, n_radrefwvl, n_irradwvl, nclenfit

  ! --------------------------------
  ! Current cross-track pixel number
  ! --------------------------------
  INTEGER (KIND=i4) :: curr_xtrack_pixnum

  LOGICAL, ALLOCATABLE, DIMENSION(:) ::  cross_track_skippix

  ! ------------------------------------
  ! Number of significant digits to keep
  ! ------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_roff_dig = 5

CONTAINS

  SUBROUTINE allocate_variables (nx,nt,nw,errstat)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER(KIND=i4), INTENT(IN) :: nx,nt,nw
    ! Modified varaible
    INTEGER(KIND=i4), INTENT(INOUT) :: errstat

    errstat = 0
    IF (.NOT. ALLOCATED(spacecraft_alt)) THEN
       ALLOCATE(spacecraft_alt(0:nt-1),time(0:nt-1),xtrflg(1:nx,0:nt-1),instrument_flag(0:nt-1), &
            latitude(1:nx,0:nt-1), latitudecorner(1:4,1:nx,0:nt-1), longitude(1:nx,0:nt-1), &
            longitudecorner(1:4,1:nx,0:nt-1), szenith(1:nx,0:nt-1), &
            sazimuth(1:nx,0:nt-1),vzenith(1:nx,0:nt-1),vazimuth(1:nx,0:nt-1), &
            razimuth(1:nx,0:nt-1), nwav_irrad(1:nx), nwav_rad(1:nx,0:nt-1), &
            irradiance_wavl(1:nw,1:nx),irradiance_spec(1:nw,1:nx),irradiance_wght(1:nw,1:nx), &
            ins_sol_wav_avg(1:nx), ccdpix_selection(1:nx,1:4), ccdpix_exclusion(1:nx,1:2), &
            ccdpix_selection_rad(1:nx,0:nt-1,1:4), ccdpix_exclusion_rad(1:nx,0:nt-1,1:2), &
            radiance_wavl(1:nw,1:nx,0:nt-1), radiance_spec(1:nw,1:nx,0:nt-1), &
            radiance_prec(1:nw,1:nx,0:nt-1), radiance_qflg(1:nw,1:nx,0:nt-1), &
            yn_process_pixel(1:nx,0:nt-1), utc_time(0:nt-1), &
            solcal_itnum(1:nx), radcal_itnum(1:nx), solcal_xflag(1:nx), &
            radcal_xflag(1:nx), radcal_chisq(1:nx), solcal_pars(1:max_calfit_idx,1:nx), &
            radcal_pars(1:max_calfit_idx,1:nx), n_ins_database_wvl(1:nx), &
            ins_database(1:max_rs_idx,1:nw,1:nx), ins_database_wvl(1:nw,1:nx), &
            database(max_rs_idx, 1:nw), cross_track_skippix(1:nx), itnum_flag(1:nx,0:nt-1), &
            fitconv_flag(1:nx,0:nt-1), column_amount(1:nx,0:nt-1), column_uncert(1:nx,0:nt-1), &
            fit_rms(1:nx,0:nt-1), radfit_chisq(1:nx,0:nt-1), stat=errstat)
       ! Initialize some variables
       cross_track_skippix = .FALSE. !Only bad pixels will be changed to TRUE
    ENDIF
    IF (ctrvar%yn_radiance_reference) THEN
       IF (.NOT. ALLOCATED(radref_sza)) THEN
          ALLOCATE (radref_sza(1:nx), radref_vza(1:nx), nwav_radref(1:nx), &
               radref_spec(1:nw,1:nx), radref_wavl(1:nw,1:nx), radref_wght(1:nw,1:nx), &
               radref_qflg(1:nw,1:nx), radref_pars(1:max_calfit_idx,1:nx), &
               radref_xflag(1:nx), radref_itnum(1:nx), radref_chisq(1:nx), &
               radref_col(1:nx), radref_dcol(1:nx), radref_rms(1:nx), radref_xtrcol(1:nx), stat=errstat)
       END IF
    END IF
    IF (ctrvar%yn_common_iter) THEN
       IF (.NOT. ALLOCATED(common_cnt)) THEN
          ALLOCATE (common_cnt(1:nx), common_wvl(1:nx,1:nw), common_spc(1:nx,1:nw), stat=errstat)
       END IF
    END IF

  END SUBROUTINE allocate_variables  

  SUBROUTINE deallocate_variables (errstat)

    IMPLICIT NONE

    ! Modified varaible
    INTEGER(KIND=i4), INTENT(INOUT) :: errstat

    errstat = 0
    IF (ALLOCATED(spacecraft_alt)) THEN
       DEALLOCATE(spacecraft_alt,time,xtrflg,instrument_flag, latitude, latitudecorner, &
            longitude, longitudecorner, szenith, sazimuth, vzenith, vazimuth, &
            razimuth, nwav_irrad, nwav_rad, irradiance_wavl, irradiance_spec, irradiance_wght, &
            ins_sol_wav_avg, ccdpix_selection, ccdpix_exclusion, ccdpix_selection_rad, &
            ccdpix_exclusion_rad, radiance_wavl, &
            radiance_spec, radiance_prec, radiance_qflg, yn_process_pixel, utc_time, radcal_chisq, &
            solcal_itnum, radcal_itnum, solcal_xflag, radcal_xflag, &
            solcal_pars, radcal_pars, n_ins_database_wvl, ins_database, &
            ins_database_wvl, database, cross_track_skippix, itnum_flag, fitconv_flag, &
            column_amount, column_uncert, fit_rms, radfit_chisq, stat=errstat)
    ENDIF
    IF (ctrvar%yn_radiance_reference) THEN
       IF (ALLOCATED(radref_sza)) THEN
          DEALLOCATE(radref_sza, radref_vza, nwav_radref, radref_spec, radref_wavl, radref_wght, &
               radref_qflg, radref_pars, radref_xflag, radref_itnum, radref_chisq, radref_col, &
               radref_dcol, radref_rms, radref_xtrcol, stat=errstat)
       END IF
    END IF
    IF (ctrvar%yn_common_iter) THEN
       IF (ALLOCATED(common_cnt)) THEN
          DEALLOCATE(common_cnt, common_wvl, common_spc, stat=errstat)
       END IF
    END IF

  END SUBROUTINE deallocate_variables  


END MODULE OMSAO_data_module
