MODULE OMSAO_omidata_module

  USE OMSAO_precision_module, ONLY: i1, i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: maxchlen, max_spec_pts
  USE OMSAO_indices_module, ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx, o3_t1_idx, o3_t3_idx
  IMPLICIT NONE

  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       nxtrack_max    =  60, nwavel_max = 300, &
       nwavelcoef_max =   5, nlines_max = 500, &
       nUTCdim        =   27

  ! -------------------------------------------------------
  ! Maximum dimension of OMI data fields in SAO L2 products
  ! -------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_field_maxdim = 3

  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       NISE_snowfree =   0, NISE_allsnow = 100, NISE_permice = 101, NISE_drysnow = 103, &
       NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! Values to go into the diagnostic array that shows how the AMF was
  ! computed, with values that indicate missing cloud products, glint,
  ! and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       omi_cfr_addmiss = 1000, omi_ctp_addmiss = 2000, omi_glint_add = 10000, &
       omi_geo_amf = -1, omi_oobview_amf = -2, omi_wfmod_amf = -9, omi_bigsza_amf = -3

  ! ------------------------------------------------------------------
  ! Arrays for L1b data (radiance, radiance reference, and irradiance)
  ! ------------------------------------------------------------------
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:) :: spacecraft_alt, spacecraft_alt_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: time, time_reference
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: radiance_errstat, instrument_flag, &
       radiance_errstat_reference, instrument_flag_reference
  INTEGER (KIND=i1), ALLOCATABLE, DIMENSION (:,:) :: xtrflg_l1b, xtrflg_l1b_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: geoflg, xtrflg, geoflg_reference, xtrflg_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: height, land_water_flg, snowicefraction, &
       height_reference, land_water_flg_reference, snowicefraction_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: latitute, longitude, latitute_reference, longitude_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: szenith, sazimuth, szenith_reference, sazimuth_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: vzenith, vazimuth, vzenith_reference, vazimuth_reference
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: razimuth, azimuth_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_spec, radiance_spec_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_prec, radiance_prec_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_wavl, radiance_wavl_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_qflg, radiance_qflg_reference
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:,:) :: radiance_ccdpix, radiance_ccdpix_reference
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: time_utc, time_utc_reference
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: irradiance_prec, irradiance_wavl, irradiance_spec, &
       irradiance_wght, radref_spec, radref_wavl, radref_wght

  ! -------------------------------------------------------------------
  ! Arrays for irradiance or radiance reference and ccd pixel selection
  ! -------------------------------------------------------------------
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: irradiance_ccdpix
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: irradiance_qflg, radref_qflg
  REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:) :: radref_sza, radref_vza
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION(:,:) :: ccdpix_selection
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION(:,:) :: ccdpix_exclusion

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

  ! --------------------------------------------------------
  ! Ozone is a special case: We can have up to 3 temperatues
  ! --------------------------------------------------------
  REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx, nxtrack_max,0:nlines_max-1) :: o3_amount, o3_uncert

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  INTEGER (KIND=i4) :: nwavel, ntimes_loop
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: nwav_irrad, nwav_radref
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: nwav_rad, nwav_rad_reference

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: n_ins_database_wvl
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:) :: solcal_itnum, radcal_itnum, radref_itnum, &
       solcal_xflag, radcal_xflag, radref_xflag
  REAL    (KIND=r8), DIMENSION (max_calfit_idx, nxtrack_max) :: solcal_pars, radcal_pars, radref_pars
  REAL    (KIND=r8), DIMENSION (max_rs_idx, nwavel_max, nxtrack_max) :: ins_database
  REAL    (KIND=r8), DIMENSION (            nwavel_max, nxtrack_max) :: ins_database_wvl
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: omi_sol_wav_avg
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:) :: solcal_chisq, radcal_chisq, radref_chisq, &
       radref_col, radref_dcol, radref_rms, radref_xtrcol
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max,0:nlines_max-1)        :: omi_wavwin_rad, omi_fitwin_rad
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max)                       :: omi_wavwin_sol, omi_fitwin_sol

  ! ------------------------------
  ! Distance between Earth and Sun
  ! ------------------------------
  REAL (KIND=r4) :: EarthSunDistance

  ! ---------------------------------------------------------
  ! OMI scan line, block line, and across-track pixel numbers
  ! ---------------------------------------------------------
  INTEGER (KIND=i4) :: omi_scanline_no, omi_blockline_no, omi_xtrackpix_no

  ! ---------------------------
  ! OMI L2 output data QA flags
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
  INTEGER (KIND=i4) :: n_omi_radwvl, n_radrefwvl, n_omi_irradwvl, nclenfit

  ! --------------------------------
  ! Current cross-track pixel number
  ! --------------------------------
  INTEGER (KIND=i4) :: curr_xtrack_pixnum

  LOGICAL, DIMENSION (nxtrack_max) ::  omi_cross_track_skippix = .FALSE.

  ! ------------------------------------
  ! Number of significant digits to keep
  ! ------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_roff_dig = 5

CONTAINS

  SUBROUTINE allocate_radiance_variables (nx,nt,nw,errstat)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER(KIND=i4), INTENT(IN) :: nx,nt,nw
    ! Modified varaible
    INTEGER(KIND=i4), INTENT(INOUT) :: errstat

    errstat = 0
    IF (.NOT. ALLOCATED(spacecraft_alt)) THEN
       print*, 'Allocating...'
       ALLOCATE(spacecraft_alt(0:nt-1),xtrflg(1:nx,0:nt-1),instrument_flag(0:nt-1), &
            latitute(1:nx,0:nt-1),longitude(1:nx,0:nt-1), szenith(1:nx,0:nt-1), &
            sazimuth(1:nx,0:nt-1),vzenith(1:nx,0:nt-1),vazimuth(1:nx,0:nt-1), &
            razimuth(1:nx,0:nt-1), nwav_irrad(1:nx), nwav_rad(1:nx,0:nt-1), &
            irradiance_wavl(1:nw,nx),irradiance_spec(1:nw,1:nx),omi_sol_wav_avg(1:nx), &
            ccdpix_selection(nx,4), ccdpix_exclusion(nx,2), radiance_wavl(1:nw,1:nx,0:nt-1), &
            radiance_spec(1:nw,1:nx,0:nt-1), radiance_prec(1:nw,1:nx,0:nt-1), &
            radiance_qflg(1:nw,1:nx,0:nt-1), stat=errstat)
    ENDIF

  END SUBROUTINE allocate_radiance_variables  

END MODULE OMSAO_omidata_module
