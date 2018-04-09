MODULE OMSAO_database_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: normweight, r8_missval, &
       zerospec_string, solar_i0_scd, yn_i0_spc
  USE OMSAO_indices_module, ONLY: max_rs_idx, hwe_idx, asy_idx, &
       sha_idx, solar_idx, ring_idx, mxs_idx, max_calfit_idx, & 
       refspec_strings, comm_idx, us1_idx, us2_idx
  USE OMSAO_variables_module, ONLY: hw1e, e_asym, g_shap, &
       sol_wav_avg, database, ctrvar, fit_winwav_idx, refspecs_original, &
       common_mode_spec
  USE OMSAO_errstat_module, ONLY: omsao_w_interpol, omsao_e_interpol, &
       omsao_w_interpol_range, pge_errstat_error, pge_errstat_ok, &
       pge_errstat_warning, vb_lev_default, vb_lev_develop, error_check, &
       f_sep, omsao_w_interpol_range, omsao_e_interpol_refspec 

  IMPLICIT NONE

  ! Variables for high resolution convolved reference spectra
  INTEGER (KIND=i4) :: nhr
  REAL (KIND=r8) :: hrd
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:,:) :: ins_hr_database
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: hr_database
  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:) :: hr_grid
  

CONTAINS
  SUBROUTINE xtrack_prepare_database ( &
       first_pix, last_pix, n_max_rspec, n_comm_wvl_out, errstat )
    
    USE OMSAO_data_module, ONLY: cross_track_skippix, nwav_radref, &
         n_ins_database_wvl, solcal_pars, ins_sol_wav_avg,  &
         nwav_irrad, irradiance_wght, irradiance_wavl, &
         irradiance_spec, ins_database, ins_database_wvl, &
         n_irradwvl, nxtrack_rad
    
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
            hr_database(1:max_rs_idx,1:nhr), &
            ins_hr_database(1:max_rs_idx,1:nhr,1:nxtrack_rad), stat = locerrstat)
       ins_hr_database=r8_missval
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
       hr_database = r8_missval
       
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
       ins_hr_database(1:max_rs_idx,1:nhr,ipix) = &
            hr_database(1:max_rs_idx,1:nhr)
       ins_database(1:max_rs_idx,1:n_ref_wvl,ipix) = &
            database(1:max_rs_idx,1:n_ref_wvl)
       ins_database_wvl(1:n_ref_wvl,ipix) = ref_wvl(1:n_ref_wvl)
       n_ins_database_wvl(ipix) = n_ref_wvl     
       
    END DO XTrackWavCal
    
    errstat = MAX ( errstat, locerrstat )
    
    RETURN
  END SUBROUTINE xtrack_prepare_database
  
  SUBROUTINE prepare_databases ( xtrack_pix, n_sol_wvl, sol_wvl, sol_spc, &
       n_rad_wvl, curr_rad_wvl, n_max_rspec, errstat )

    ! ===========================================
    !
    ! Called from XTRACK_RADIANCE_WVL_CALIBRATION
    !
    ! ===========================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER(KIND=i4), INTENT (IN) :: xtrack_pix, n_rad_wvl, n_sol_wvl, &
         n_max_rspec
    REAL(KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: curr_rad_wvl
    REAL(KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: sol_wvl, sol_spc
    
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
    ! -------------------------------------------------------------------------------------
    ! Calculate the splined fitting database. This will be saved into a nXtrack-dimensional
    ! array in a calling routine higher up (inside the fitting loop). 
    ! -------------------------------------------------------------------------------------
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
    
    ! --------------------------------------------------
    ! Correct the reference solar spectra for bad pixels
    ! --------------------------------------------------
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

  SUBROUTINE dataspline ( xtrack_pix, n_radwvl, curr_rad_wvl, n_max_rspec, errstat )
    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    xtrack_pix ........... current cross-track pixel number
    !    n_radwvl ............. number of radiance wavelenghts
    !    curr_rad_wvl ......... radiance wavelengths
    !    n_max_rspec .......... maximum number of reference spectra points
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------
    USE OMSAO_data_module, ONLY : solcal_pars
    
    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                       INTENT (IN) :: xtrack_pix, n_radwvl, n_max_rspec
    REAL    (KIND=r8), DIMENSION (n_radwvl), INTENT (IN) :: curr_rad_wvl
    
    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL :: yn_full_range
    INTEGER (KIND=i4) :: idx, npts, locerrstat, iii, nsol, ios
    REAL (KIND=r8) :: DU_load
    REAL (KIND=r8), DIMENSION (n_max_rspec)    :: tmp_spec, tmp_wavl
    REAL (KIND=r8), DIMENSION (n_radwvl)       :: dbase_loc
    REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: solar_spc, solar_wvl, solar_conv, xsec_i0_spc
    
    CHARACTER (LEN=11), PARAMETER :: modulename = 'dataspline'
    
    locerrstat = pge_errstat_ok
    
    nsol = refspecs_original(solar_idx)%nPoints
    IF ( nsol > 0 ) THEN
       ALLOCATE ( solar_spc  (1:nsol), STAT=ios )
       ALLOCATE ( solar_wvl  (1:nsol), STAT=ios )
       ALLOCATE ( solar_conv (1:nsol), STAT=ios )
       ALLOCATE ( xsec_i0_spc(1:nsol), STAT=ios )
    END IF
    
    ! --------------------------------------
    ! Compute correction for Solar I0 effect
    ! --------------------------------------
    ! work out convolved high resolution
    ! solar spectrum
    ! ----------------------------------
    IF ( ctrvar%yn_solar_i0 ) THEN
       idx  = solar_idx
       solar_wvl(1:nsol) = refspecs_original(idx)%RefSpecWavs(1:nsol)
       solar_spc(1:nsol) = refspecs_original(idx)%RefSpecData(1:nsol)
       CALL convolve_data ( &
            xtrack_pix, nsol,  solar_wvl(1:nsol), solar_spc(1:nsol), ctrvar%yn_use_labslitfunc, &
            solcal_pars(hwe_idx,xtrack_pix), solcal_pars(asy_idx,xtrack_pix), &
            solcal_pars(sha_idx,xtrack_pix), solar_conv(1:nsol), errstat )
    END IF
    
    ! ---------------------------------------------------------------------
    ! Load results into the database array. The order of the spectra is
    ! determined by the molecule indices in OMSAO_indices_module. Only 
    ! those spectra that are requested for read-in are interpolated.
    ! Note that this might be different from the actual fitting parameters
    ! to be varied in the fit - a non-zero but constant fitting parameter
    ! will still require a spectrum to make a contribution. However, if all
    ! fitting parameters, including the upper and lower bounds, are ZERO,
    ! then there is no way that spectrum can contribute, and it therefore
    ! doesn't have to be interpolated.
    !
    ! If the original reference spectrum does not cover the wavelength
    ! range of the current radiance wavelength, we interpolate only the
    ! part that is covered and set the rest to Zero.
    ! ---------------------------------------------------------------------
    
    ! ------------------------
    ! Spline Reference Spectra
    ! ------------------------
    DO idx = 1, max_rs_idx
       
       ! --------------------------------------------------------------------------
       ! We don't spline the solar reference spectrum and the undersampling spectra
       ! (computed and assigned to DATABASE in the UNDERSAMPLE subroutine), and the
       ! common mode spectrum (which hasn't been defined yet anyway).
       ! --------------------------------------------------------------------------
       IF ( (idx == solar_idx) .OR. (idx == us1_idx) .OR. (idx == us2_idx) .OR. (idx == comm_idx) ) CYCLE
       
       ! -----------------------------------------------------------------------------
       ! Only carry on with the i0 correction, convolution and interpolation if we are
       ! using that particular cross-section
       ! -----------------------------------------------------------------------------
       iii = max_calfit_idx + (idx-1)*mxs_idx
       IF ( refspecs_original(idx)%nPoints /= 0                                          .AND. &
            INDEX (TRIM(ADJUSTL(refspecs_original(idx)%FileName)), zerospec_string) == 0 .AND. &
            ( ANY (ctrvar%fitvar_rad_init(iii+1:iii+mxs_idx) /= 0) .OR. &
            ANY (ctrvar%lo_radbnd      (iii+1:iii+mxs_idx) /= 0) .OR. &
            ANY (ctrvar%up_radbnd      (iii+1:iii+mxs_idx) /= 0)        ) ) THEN
          
          ! -----------------------------------------------
          ! Define short-hand for number of spectral points
          ! -----------------------------------------------
          npts = refspecs_original(idx)%nPoints
          tmp_wavl(1:npts) = refspecs_original(idx)%RefSpecWavs(1:npts)
          tmp_spec(1:npts) = refspecs_original(idx)%RefSpecData(1:npts)
          
          ! ---------------------------------------------------------
          ! Solar I0 correction, Yes or No?
          ! Check if it is a corrected species (set in params module)
          ! ---------------------------------------------------------
          IF ( ctrvar%yn_solar_i0 .AND. yn_i0_spc(idx)) THEN
             
             DU_load = solar_i0_scd(idx)
             ! --------------------------------------------------------------------
             ! 1: Interpolate cross sections to solar reference spectrum wavelength
             ! --------------------------------------------------------------------
             CALL interpolation ( &
                  npts, tmp_wavl(1:npts), tmp_spec(1:npts),                         &
                  nsol, solar_wvl(1:nsol), xsec_i0_spc(1:nsol),                     &
                  'fillvalue', 0.0_r8, yn_full_range, locerrstat )
             
             ! ---------------------
             ! 2: Undo normalization
             ! ---------------------
             xsec_i0_spc(1:nsol) = xsec_i0_spc(1:nsol) * refspecs_original(idx)%NormFactor
             
             ! ------------------------------------
             ! 3: Compute attenuated solar spectrum
             ! ------------------------------------
             tmp_spec(1:nsol) = solar_spc(1:nsol) * EXP(-xsec_i0_spc(1:nsol)*DU_load)
             
             ! ------------------------------
             ! 4: Convolve with slit function
             ! ------------------------------
             CALL convolve_data ( &
                  xtrack_pix, nsol, solar_wvl(1:nsol), tmp_spec(1:nsol), ctrvar%yn_use_labslitfunc, &
                  solcal_pars(hwe_idx,xtrack_pix), solcal_pars(asy_idx,xtrack_pix), &
                  solcal_pars(sha_idx,xtrack_pix), xsec_i0_spc(1:nsol), errstat )
             
             ! -----------------------------------
             ! 5: Compute corrected cross sections
             ! -----------------------------------
             WHERE ( solar_conv(1:nsol) > 0.0_r8 .AND. xsec_i0_spc(1:nsol) > 0.0_r8 )
                xsec_i0_spc(1:nsol) = -LOG( xsec_i0_spc(1:nsol) / solar_conv(1:nsol) ) / DU_load
             END WHERE
             
             ! ---------------------
             ! 6: Redo normalization
             ! ---------------------
             xsec_i0_spc(1:nsol) = xsec_i0_spc(1:nsol) / refspecs_original(idx)%NormFactor
             
             ! ------------------------------------------------------------
             ! Copy values to the variables used in the interpolation below
             ! ------------------------------------------------------------
             npts = nsol
             tmp_spec(1:npts) = xsec_i0_spc(1:npts)
             tmp_wavl(1:npts) = solar_wvl(1:npts)
             
          ELSE
             ! -------------------------------------------------------------------------
             ! All other cross sections just need to be convolved with the slit function,
             ! EXCEPT for the Common Mode spectra.
             ! -------------------------------------------------------------------------
             IF ( (idx == comm_idx) ) THEN
                tmp_spec(1:npts) = common_mode_spec%RefSpecData(xtrack_pix,1:npts)
             ELSE
                CALL convolve_data ( &
                     xtrack_pix, npts, tmp_wavl(1:npts), tmp_spec(1:npts), ctrvar%yn_use_labslitfunc, &
                     solcal_pars(hwe_idx,xtrack_pix), solcal_pars(asy_idx,xtrack_pix), &
                     solcal_pars(sha_idx,xtrack_pix), tmp_spec(1:npts), errstat )
             END IF
          END IF
          
          ! ----------------------------------------------------------------------------
          ! Call interpolation and check returned error status. WARNING status indicates
          ! missing parts of the interpolated spectrum, while ERROR status indicates a
          ! more serious condition that requires termination.
          ! ----------------------------------------------------------------------------
          CALL interpolation ( &
               npts, tmp_wavl(1:npts), tmp_spec(1:npts),                         &
               n_radwvl, curr_rad_wvl(1:n_radwvl), dbase_loc(1:n_radwvl),        &
               'fillvalue', 0.0_r8, yn_full_range, locerrstat )
          database(idx, 1:n_radwvl) = dbase_loc(1:n_radwvl)

          CALL error_check ( &
               locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL_REFSPEC, &
               modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_default, errstat )
          IF ( .NOT. yn_full_range ) THEN
             CALL error_check ( &
                  0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
                  modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_develop, errstat )
          END IF
          ! ----------------------------------------------------
          ! Save high resolution spectrum for each cross section
          ! They are saved in the hr_grid
          ! ----------------------------------------------------
          CALL interpolation ( &
               npts, tmp_wavl(1:npts), tmp_spec(1:npts),           &
               nhr, hr_grid(1:nhr), hr_database(idx,1:nhr),        &
               'fillvalue', 0.0_r8, yn_full_range, locerrstat )
          CALL error_check ( &
               locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL_REFSPEC, &
               modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_default, errstat )
          IF ( .NOT. yn_full_range ) THEN
             CALL error_check ( &
                  0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
                  modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_develop, errstat )
          END IF

       END IF
    END DO
    
    IF ( nsol > 0 ) THEN
       IF ( ALLOCATED ( solar_spc   ) )  DEALLOCATE ( solar_spc   )
       IF ( ALLOCATED ( solar_wvl   ) )  DEALLOCATE ( solar_wvl   )
       IF ( ALLOCATED ( solar_conv  ) )  DEALLOCATE ( solar_conv  )
       IF ( ALLOCATED ( xsec_i0_spc ) )  DEALLOCATE ( xsec_i0_spc )
    END IF
    
    RETURN
  END SUBROUTINE dataspline
  
  
  SUBROUTINE convolve_data (                                     &
       xtrack_pix, npts, wvl_in, spec_in, yn_labslit, hw1e, asy, &
       sha, spec_conv, errstat )
    
    USE OMSAO_slitfunction_module, ONLY: omi_slitfunc_convolve, &
         super_gaussian_sf
    
    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: xtrack_pix, npts
    REAL    (KIND=r8),                   INTENT (IN) :: hw1e, asy, sha
    LOGICAL,                             INTENT (IN) :: yn_labslit
    REAL    (KIND=r8), DIMENSION (npts), INTENT (IN) :: spec_in, wvl_in
    
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: spec_conv
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)           :: errstat
    CHARACTER (LEN=13), PARAMETER :: modulename = 'convolve_data'
    
    
    errstat = pge_errstat_ok
    
    ! -----------------------------------------------------------------
    ! Either laboratory slit function (tabulated) or Gaussian (fitted).
    ! -----------------------------------------------------------------
    IF ( yn_labslit ) THEN
       CALL omi_slitfunc_convolve (                     &
            xtrack_pix, npts, wvl_in(1:npts), spec_in(1:npts), spec_conv(1:npts), sha, errstat )
       CALL error_check ( &
            errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_INTERPOL, &
            modulename//f_sep//'Convolution', vb_lev_default, errstat )
    ELSE
       ! -----------------------------------------------------------------
       ! Here is the Gaussian branch. We need to make sure that we use the
       ! slit function information for the current pixel. There is also no
       ! need to save the convolved spectrum, because we have to convolve
       ! for each and every pixel due to the varying slit function.
       ! -----------------------------------------------------------------
       CALL super_gaussian_sf (                                           &
            npts, hw1e, asy, sha, wvl_in(1:npts), spec_in(1:npts), spec_conv(1:npts) )
    END IF
    
    RETURN
  END SUBROUTINE convolve_data

  SUBROUTINE deallocate_hr_database (errstat)

    IMPLICIT NONE

    ! Modified variable
    INTEGER(KIND=i4), INTENT(INOUT) :: errstat
    errstat = 0

    IF (ALLOCATED(hr_grid)) DEALLOCATE (hr_grid, hr_database, ins_hr_database, stat=errstat)

  END SUBROUTINE deallocate_hr_database

END MODULE OMSAO_database_module
