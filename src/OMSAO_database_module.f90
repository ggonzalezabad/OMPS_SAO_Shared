MODULE OMSAO_database_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: normweight, r8_missval
  USE OMSAO_indices_module, ONLY: max_rs_idx, hwe_idx, asy_idx, &
       sha_idx, solar_idx, ring_idx
  USE OMSAO_variables_module, ONLY: hw1e, e_asym, g_shap, &
       sol_wav_avg, database, ctrvar, fit_winwav_idx, refspecs_original
  USE OMSAO_errstat_module, ONLY: omsao_e_interpol, omsao_w_interpol_range, &
       pge_errstat_error, pge_errstat_ok, pge_errstat_warning, vb_lev_default, &
       vb_lev_develop, error_check
  

  IMPLICIT NONE

  ! Variables for high resolution convolved reference spectra
  INTEGER (KIND=i4) :: nhr
  REAL (KIND=r8) :: hrd
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: hr_database
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:) :: hr_grid

CONTAINS
  SUBROUTINE xtrack_prepare_database ( &
       first_pix, last_pix, n_max_rspec, n_comm_wvl_out, errstat )
    
    USE OMSAO_data_module, ONLY: cross_track_skippix, nwav_radref, &
         n_ins_database_wvl, solcal_pars, ins_sol_wav_avg,  &
         nwav_irrad, irradiance_wght, irradiance_wavl, &
         irradiance_spec, ins_database, ins_database_wvl, &
         n_irradwvl
    
    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix, n_max_rspec
    
    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (OUT) :: n_comm_wvl_out
    
    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, n_ref_wvl, igrid
    REAL (KIND=r8), DIMENSION (n_max_rspec) :: ref_wvl, ref_spc, ref_wgt
    
    locerrstat = pge_errstat_ok
    
    ! -------------------------------------------------
    ! If not allocated highres_database variables do it
    ! To be used later on interpolation for each pixel
    ! -------------------------------------------------
    IF (.NOT. ALLOCATED(hr_grid)) THEN
       ! Work out number of points of the grid
       hrd = refspecs_original(solar_idx)%RefSpecWavs(2) - &
            refspecs_original(solar_idx)%RefSpecWavs(1)
       nhr = INT((ctrvar%fit_winwav_lim(4) - ctrvar%fit_winwav_lim(1)) / &
            hrd, KIND=i4) + 1
       ! Allocate variables
       ALLOCATE (hr_grid(1:nhr), &
            hr_database(1:max_rs_idx,1:nhr), stat = locerrstat)
       ! Fill up values of highres_grid
       DO igrid = 1,nhr
          hr_grid(igrid) = ctrvar%fit_winwav_lim(1) + &
               hrd * REAL((igrid-1), KIND=r8)
       END DO
    END IF
    
    ! -------------------------------------------------
    ! Set the number of wavelengths for the common mode
    ! It has to be the maximum between nwav_radref &
    ! nwav_rad
    ! -------------------------------------------------
    IF (ctrvar%yn_radiance_reference) THEN
       n_comm_wvl_out = MAXVAL(nwav_radref)
    ELSE
       n_comm_wvl_out = MAXVAL(n_ins_database_wvl)
    END IF
    
    ! --------------------------------
    ! Loop over cross-track positions. 
    ! --------------------------------
    XTrackWavCal: DO ipix = first_pix, last_pix
       
       IF (cross_track_skippix(ipix)) CYCLE
       ! -------------------
       ! Initialize database
       ! -------------------
       database = r8_missval
       
       ! ----------------------------------------------------
       ! Assign number of radiance and irradiance wavelengths
       ! ----------------------------------------------------
       n_irradwvl = nwav_irrad(ipix)
       
       ! ---------------------------------------------------------------
       ! Restore solar fitting variables for across-track reference in
       ! Earthshine fitting.
       ! ---------------------------------------------------------------
       sol_wav_avg = ins_sol_wav_avg(ipix)
       hw1e        = solcal_pars(hwe_idx,ipix)
       e_asym      = solcal_pars(asy_idx,ipix)
       g_shap      = solcal_pars(sha_idx,ipix)
       
       ! -----------------------------------------------------
       ! Assign (hopefully predetermined) "reference" weights.
       ! -----------------------------------------------------
       IF ( .NOT. ctrvar%yn_solar_comp ) THEN
          ref_wgt(1:n_irradwvl) = irradiance_wght(1:n_irradwvl,ipix)
       ELSE
          ref_wgt(1:n_irradwvl) = normweight
       END IF
       
       ! --------------------------------------------------------
       ! Save wavelength calibrated irradiance spectra
       ! to be used when creating the database of reference cross
       ! sections (below).
       ! ---------------------------------------------------------
       ! Solar irradiance
       n_ref_wvl = n_irradwvl
       ref_wvl(1:n_ref_wvl) = irradiance_wavl(1:n_ref_wvl,ipix)
       ref_spc(1:n_ref_wvl) = irradiance_spec(1:n_ref_wvl,ipix)
       
       ! ----------------------------------------------------
       ! Spline reference spectra to current wavelength grid.
       ! ----------------------------------------------------
       Call prepare_databases ( &
            ipix, n_ref_wvl, ref_wvl(1:n_ref_wvl), ref_spc(1:n_ref_wvl), &
            n_ref_wvl, ref_wvl(1:n_ref_wvl), n_max_rspec, locerrstat )
       
       IF ( locerrstat >= pge_errstat_error ) EXIT XTrackWavCal
       
       ! ---------------------------------------------------------
       ! Save DATABASE in OMI_DATABASE for radiance fitting loops.
       ! ---------------------------------------------------------
       ins_database(1:max_rs_idx,1:n_ref_wvl,ipix) = database(1:max_rs_idx,1:n_ref_wvl)
       ins_database_wvl(1:n_ref_wvl,ipix) = ref_wvl(1:n_ref_wvl)
       n_ins_database_wvl(ipix) = n_ref_wvl     
       
    END DO XTrackWavCal
    
    errstat = MAX ( errstat, locerrstat )
    
    RETURN
  END SUBROUTINE xtrack_prepare_database
  
  SUBROUTINE prepare_databases ( &
       xtrack_pix, n_sol_wvl, sol_wvl, sol_spc, n_rad_wvl, curr_rad_wvl, n_max_rspec, errstat )

    ! ===========================================
    !
    ! Called from XTRACK_RADIANCE_WVL_CALIBRATION
    !
    ! ===========================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: xtrack_pix, n_rad_wvl, n_sol_wvl, n_max_rspec
    REAL    (KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: curr_rad_wvl
    REAL    (KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: sol_wvl, sol_spc
    
    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    
    locerrstat = pge_errstat_ok
    
    ! ---------------------------------------------------------
    ! Spline external reference spectra to common radiance grid
    ! ---------------------------------------------------------
    CALL dataspline ( xtrack_pix, n_rad_wvl, curr_rad_wvl(1:n_rad_wvl), n_max_rspec, locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( errstat >= pge_errstat_error ) RETURN
    
    ! -----------------------------------
    ! Calculate the undersampled spectrum
    ! -----------------------------------
    IF ( ANY (ctrvar%have_undersampling) ) &
         CALL undersample ( xtrack_pix, n_rad_wvl, curr_rad_wvl(1:n_rad_wvl), &
         n_sol_wvl, sol_wvl(1:n_sol_wvl), hw1e, e_asym, g_shap, ctrvar%phase, locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( errstat >= pge_errstat_error ) RETURN
    
    ! -------------------------------------------------------------------------------------
    ! Calculate the splined fitting database. This will be saved into a nXtrack-dimensional
    ! array in a calling routine higher up (inside the fitting loop). 
    ! -------------------------------------------------------------------------------------
    CALL prepare_solar_refspec ( &
         n_sol_wvl, sol_wvl(1:n_sol_wvl), sol_spc(1:n_sol_wvl), n_rad_wvl, &
         curr_rad_wvl(1:n_rad_wvl), locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( errstat >= pge_errstat_error ) RETURN
    
    RETURN
  END SUBROUTINE prepare_databases

  SUBROUTINE prepare_solar_refspec ( &
       n_solpts, sol_wvl, sol_spc, n_radpts, curr_rad_wvl, errstat )
    
    ! ***********************************************************
    !
    !   Spline the solar measured solar spectrum to the radiance
    !   wavelength grid
    !
    !   Called from PREPARE_DATABASE
    !
    ! ***********************************************************  
    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                       INTENT (IN)    :: n_radpts, n_solpts
    REAL    (KIND=r8), DIMENSION (n_solpts), INTENT (IN)    :: sol_wvl, sol_spc
    REAL    (KIND=r8), DIMENSION (n_radpts), INTENT (IN)    :: curr_rad_wvl
    INTEGER (KIND=i4),                       INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL                                 :: yn_full_range
    INTEGER (KIND=i4)                       :: j, ll_rad, lu_rad, locerrstat, n_sol_tmp
    REAL    (KIND=r8), DIMENSION (n_radpts) :: spline_sun
    REAL    (KIND=r8), DIMENSION (n_solpts) :: tmp_sol_spec, tmp_sol_wvl
    
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'prepare_solar_refspec'
    
    
    locerrstat = pge_errstat_ok
    
    ! ---------------------------------------------
    ! Spline irradiance spectrum onto radiance grid
    ! ---------------------------------------------
    ! There are two issues here:
    ! * The radiance and irradiance grids are likely to be off-set, i.e., at one
    !   or both end points we will run into an extrapolation.
    ! * Any bad pixel in the solar spectrum must be avoided in the interpolation.
    
    ! ----------------------------------------------------------------------
    ! First create an array with solar spectrum and wavelengths that is free
    ! of any bad pixels. "Bad" pixels are set to -1.0_r8
    ! ----------------------------------------------------------------------
    n_sol_tmp = 0
    DO j = 1, n_solpts
       IF ( sol_spc(j) > 0.0_r8 ) THEN
          n_sol_tmp = n_sol_tmp + 1
          tmp_sol_wvl (n_sol_tmp) = sol_wvl(j)
          tmp_sol_spec(n_sol_tmp) = sol_spc(j)
          ! ------------------------------------------------
          ! Check whether wavelengths are strictly ascending
          ! ------------------------------------------------
          IF ( n_sol_tmp > 1 ) THEN
             IF ( tmp_sol_wvl (n_sol_tmp) <= tmp_sol_wvl (n_sol_tmp-1) ) n_sol_tmp = n_sol_tmp - 1
          END IF
       END IF
    END DO
    
    ! ----------------------------------------------------------------
    ! Interpolate (spectral overlap is taken care of in INTERPOLATION)
    ! ----------------------------------------------------------------
    CALL interpolation ( &
         n_sol_tmp, tmp_sol_wvl(1:n_sol_tmp), tmp_sol_spec(1:n_sol_tmp),             &
         n_radpts, curr_rad_wvl(1:n_radpts), spline_sun(1:n_radpts),                 &
         'endpoints', 0.0_r8, yn_full_range, locerrstat )
    CALL error_check ( &
         locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
         modulename, vb_lev_default, errstat )
    IF ( errstat >= pge_errstat_error ) RETURN
    IF ( .NOT. yn_full_range )   CALL error_check ( 0, 1, pge_errstat_warning, &
         OMSAO_W_INTERPOL_RANGE, modulename, vb_lev_develop, errstat )
    
    ! --------------------------------------------------------------------------
    ! Finally, we have to check whether we encountered any bad solar pixels. The
    ! adopted strategy is to exclude any such pixels from the fit, plus a window
    ! of +/- 0.1 nm around them.
    ! --------------------------------------------------------------------------
    DO j = 1, n_solpts
       IF ( sol_spc(j) < 0.0_r8 ) THEN
          WHERE ( ABS(curr_rad_wvl(1:n_radpts)-sol_wvl(j)) <= 0.1_r8)
             spline_sun(1:n_radpts) = -1.0_r8
          END WHERE
       END IF
    END DO
    
    ! --------------------------------------------------------
    ! Now we are ready to assign the splined solar spectrum to
    ! its final array - DATABASE(solar_idx,*).
    ! --------------------------------------------------------
    database(solar_idx,1:n_radpts) = spline_sun(1:n_radpts)
    
    ! =================================================================
    ! Note that the UNDERSAMPLING spectrum has already been assigned to
    ! DATABASE(us1/2_idx,*) in the UNDERSPEC routine.
    ! =================================================================
    
    ! ---------------------------------------------------------------------
    ! Set up Ring spectrum for high pass filtering of divided spectrum
    ! (afterward, add solar spectrum back in order to divide by solar
    ! spectrum with altered wavelength calibration in subroutine spectrum).
    ! ---------------------------------------------------------------------
    IF ( ctrvar%yn_doas ) database(ring_idx, 1:n_radpts) = &
         database(ring_idx, 1:n_radpts) / spline_sun(1:n_radpts)
    
    ! ---------------------------------------------------------
    ! For the DOAS case, high-pass filter the reference spectra
    ! ---------------------------------------------------------
    IF ( ctrvar%yn_doas ) THEN
       ll_rad = fit_winwav_idx(2) ; lu_rad = fit_winwav_idx(3)
       CALL subtract_cubic (curr_rad_wvl, n_radpts, ll_rad, lu_rad)
       database(ring_idx, 1:n_radpts) = &
            database(ring_idx, 1:n_radpts) * spline_sun(1:n_radpts)
    END IF
    
    IF ( ctrvar%yn_smooth ) database(1:max_rs_idx, 3:n_radpts-2) = &
         0.375_r8  * database(1:max_rs_idx, 3:n_radpts-2)  +  &
         0.25_r8   * (database(1:max_rs_idx, 4:n_radpts-1) + &
         database(1:max_rs_idx, 2:n_radpts-3)) + &
         0.0625_r8 * (database(1:max_rs_idx, 5:n_radpts)   + &
         database(1:max_rs_idx, 1:n_radpts-4))
    
    errstat = MAX ( errstat, locerrstat )
    
    RETURN
  END SUBROUTINE prepare_solar_refspec

END MODULE OMSAO_database_module
