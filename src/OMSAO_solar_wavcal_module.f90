MODULE OMSAO_solar_wavcal_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: &
       wvl_idx, sig_idx, spc_idx, ccd_idx, max_calfit_idx, shi_idx, squ_idx, solcal_idx,  &
       hwe_idx, asy_idx, sha_idx, solar_idx, bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, &
       bl5_idx, sc0_idx, sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx, sin_idx
  USE OMSAO_parameters_module, ONLY: &
       i2_missval, i4_missval, r8_missval, normweight, downweight, maxchlen, max_spec_pts
  USE OMSAO_variables_module,  ONLY: &
       hw1e, e_asym, g_shap, curr_sol_spec, sol_wav_avg, fitvar_cal, fitvar_cal_saved,  &
       mask_fitvar_cal, n_fitvar_cal, lobnd, upbnd,&
       fitwavs, fitweights, currspec, refspecs_original,    &
       pcfvar, ctrvar
  USE OMSAO_data_module, ONLY: nxtrack_rad, nwav_irrad, solcal_itnum, solcal_xflag, &
       ins_sol_wav_avg, cross_track_skippix, solcal_pars, n_irradwvl, irradiance_wght, &
       irradiance_spec, irradiance_wavl, curr_xtrack_pixnum, &
       ccdpix_selection, ccdpix_exclusion
  USE OMSAO_errstat_module
  USE OMSAO_he5_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE xtrack_solar_calibration_loop ( first_pix, last_pix, errstat )

    IMPLICIT NONE
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i2) :: local_solcal_itnum
    INTEGER (KIND=i4) :: locerrstat, ipix, solcal_exval, i, ngood, nbad
    CHARACTER (LEN=maxchlen) :: addmsg
    REAL (KIND=r8) :: chisquav, sol_spec_avg, asum, ssum
    LOGICAL :: yn_skip_pix, yn_bad_pixel
    INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:) :: bad_idx
    REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: fitres_out, weightsum, &
         wvl_good, wvl_bad, spc_good, spc_bad

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=29), PARAMETER :: modulename = 'xtrack_solar_calibration_loop'

    fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_sol_init(1:max_calfit_idx)

    ! ---------------------------------------------------------------
    ! Loop for solar wavelength calibration and slit function fitting
    ! ---------------------------------------------------------------
    XtrackSolCal: DO ipix = first_pix, last_pix

       locerrstat = pge_errstat_ok
       curr_xtrack_pixnum = ipix

       n_irradwvl = nwav_irrad(ipix)
       IF ( n_irradwvl <= 0 ) CYCLE
       yn_skip_pix = .FALSE.
       ! ------------------------
       ! Allocate local variables
       ! ------------------------
       ALLOCATE(fitres_out(1:n_irradwvl), weightsum(1:n_irradwvl), wvl_good(1:n_irradwvl), &
            wvl_bad(1:n_irradwvl), spc_good(1:n_irradwvl), spc_bad(1:n_irradwvl), bad_idx(1:n_irradwvl) )

       ! -----------------------------
       ! Set up generic fitting arrays
       ! -----------------------------
       curr_sol_spec(wvl_idx,1:n_irradwvl) = irradiance_wavl(1:n_irradwvl,ipix)
       curr_sol_spec(spc_idx,1:n_irradwvl) = irradiance_spec(1:n_irradwvl,ipix)
       curr_sol_spec(sig_idx,1:n_irradwvl) = normweight

       ! -----------------------------------
       ! Make sure wavelengths are ascending
       ! -----------------------------------
       DO i = 2, n_irradwvl
          IF ( curr_sol_spec(wvl_idx,i) <= curr_sol_spec(wvl_idx,i-1) ) THEN
             curr_sol_spec(wvl_idx,i) = curr_sol_spec(wvl_idx,i-1) + 0.001_r8
             curr_sol_spec(sig_idx,i) = downweight
          END IF
       END DO

       ! ----------------------------------------
       ! No missing / negative values in spectrum
       ! ----------------------------------------
       WHERE ( curr_sol_spec(spc_idx,1:n_irradwvl) <= 0.0 )
          curr_sol_spec(sig_idx,1:n_irradwvl) = downweight
          curr_sol_spec(spc_idx,1:n_irradwvl) = 0.0_r8
       END WHERE

       ! -----------------------------------------------
       ! Compute normalization factor for solar spectrum
       ! -----------------------------------------------
       weightsum = 0.0_r8
       WHERE ( curr_sol_spec(sig_idx,1:n_irradwvl) /= downweight )
          weightsum = 1.0_r8
       END WHERE
       sol_spec_avg = SUM ( curr_sol_spec(spc_idx,1:n_irradwvl)*weightsum(1:n_irradwvl) ) / &
            MAX(1.0_r8, SUM(weightsum(1:n_irradwvl)))
       IF ( sol_spec_avg == 0.0_r8 ) sol_spec_avg = 1.0_r8

       ! -------------------------------------------------------------------------
       ! So far we have only taken care of/excluded any negative values in the
       ! spectrum, but there may abnormally high or low positive values also.
       ! Now we check for any values exceeding 100 times the average, which should
       ! be a large enough window to keep anything sensible and reject the real
       ! outliers.
       ! -------------------------------------------------------------------------
       WHERE ( weightsum(1:n_irradwvl) /= 0.0_r8 .AND. &
            ABS(curr_sol_spec(spc_idx,1:n_irradwvl)) >= 100.0_r8 * sol_spec_avg )
          weightsum(1:n_irradwvl) = 0.0_r8
          curr_sol_spec(sig_idx,1:n_irradwvl) = downweight
          curr_sol_spec(spc_idx,1:n_irradwvl) = 0.0_r8
       ENDWHERE

       ! ------------------------------------------------------------------
       ! Recompute the solar spectrum average, because it may have changed.
       ! ------------------------------------------------------------------
       sol_spec_avg = SUM ( curr_sol_spec(spc_idx,1:n_irradwvl)*weightsum(1:n_irradwvl) ) / &
            MAX(1.0_r8, SUM(weightsum(1:n_irradwvl)))
       IF ( sol_spec_avg <= 0.0_r8 ) THEN
          yn_skip_pix = .TRUE.
          sol_spec_avg = 1.0_r8
       END IF

       ! -----------------------
       ! Skip pixel if necessary
       ! -----------------------
       IF ( yn_skip_pix .OR. locerrstat >= pge_errstat_error ) THEN
          errstat = MAX ( errstat, locerrstat )
          cross_track_skippix (ipix) = .TRUE.
          addmsg = ''
          WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
          CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
               modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
               locerrstat )
          DEALLOCATE(fitres_out, weightsum, wvl_good, wvl_bad, spc_good, spc_bad, bad_idx )
          CYCLE
       END IF

       ! --------------------------------------
       ! Finally, normalize the solar spectrum.
       ! --------------------------------------
       IF ( ctrvar%yn_spectrum_norm ) &
            curr_sol_spec(spc_idx,1:n_irradwvl) = curr_sol_spec(spc_idx,1:n_irradwvl) / sol_spec_avg
       
       ! ---------------------------------------------
       ! Calculate SOL_WAV_AVG of measured solar spectra here,
       ! for use in calculated spectra.
       ! ---------------------------------------------
       asum = SUM ( curr_sol_spec(wvl_idx,1:n_irradwvl) * curr_sol_spec(sig_idx,1:n_irradwvl) )
       ssum = SUM ( 1.0_r8 * curr_sol_spec(sig_idx,1:n_irradwvl) )
       sol_wav_avg = asum / ssum

       ! -------------------------------------------------
       ! Count the good and the bad, and fill up bad pixel 
       ! with values interpolated from the good ones.
       ! -------------------------------------------------
       ngood = 0 ; nbad = 0 ; bad_idx = 0
       wvl_good = 0.0_r8 ; wvl_bad = 0.0_r8 ; spc_good = 0.0_r8 ; spc_bad = 0.0_r8
       DO i = 1, n_irradwvl
          IF ( curr_sol_spec(spc_idx,i) <= 0.0_r8 .OR. curr_sol_spec(spc_idx,i) == downweight ) THEN
             nbad          = nbad + 1
             bad_idx(nbad) = i
             wvl_bad(nbad) = curr_sol_spec(wvl_idx,i)
          ELSE
             ngood           = ngood + 1
             spc_good(ngood) = curr_sol_spec(spc_idx,i)
             wvl_good(ngood) = curr_sol_spec(wvl_idx,i)
          END IF
       END DO
       IF ( nbad > 0 ) THEN
          CALL ezspline_1d_interpolation (                      &
               ngood, wvl_good(1:ngood), spc_good(1:ngood),     &
               nbad, wvl_bad(1:nbad), spc_bad(1:nbad), locerrstat )
          DO i = 1, nbad
             curr_sol_spec(spc_idx,bad_idx(i)) = spc_bad(i)
             curr_sol_spec(sig_idx,bad_idx(i)) = downweight
          END DO
       END IF

       ! ------------------------------------------------------------
       ! Anything outside the fitting window will receive Zero weight
       ! ------------------------------------------------------------
       ! (the CCD indices are absolute positions, i.e., unlikely to be "1:n_sol_wvl")
       ! ----------------------------------------------------------------------------
       IF ( ccdpix_selection(ipix,2) > ccdpix_selection(ipix,1) ) &
            curr_sol_spec(sig_idx,1:ccdpix_selection(ipix,2)-ccdpix_selection(ipix,1)+1) = downweight
       IF ( ccdpix_selection(ipix,3) < ccdpix_selection(ipix,4) ) &
            curr_sol_spec(sig_idx,ccdpix_selection(ipix,3)-ccdpix_selection(ipix,1)+1:n_irradwvl) = downweight
              
       ! ------------------------------------------------------------------------
       ! Also any window excluded by the user (specified in fitting control file)
       ! ------------------------------------------------------------------------
       IF ( ccdpix_exclusion(ipix,1) >= 1 .AND. ccdpix_exclusion(ipix,2) >= ccdpix_exclusion(ipix,1) ) &
            curr_sol_spec(sig_idx, &
            ccdpix_exclusion(ipix,1)-ccdpix_selection(ipix,1)+1:&
            ccdpix_exclusion(ipix,2)-ccdpix_selection(ipix,1)+1) = downweight

       yn_bad_pixel   = .FALSE.
       CALL solar_fit ( & ! Solar wavelength calibration
            ctrvar%n_fitres_loop(solcal_idx), ctrvar%fitres_range(solcal_idx), n_irradwvl, &
            curr_sol_spec(wvl_idx:ccd_idx,1:n_irradwvl), hw1e, e_asym, solcal_exval, local_solcal_itnum, &
            chisquav, yn_bad_pixel, fitres_out(1:n_irradwvl), locerrstat )
       ! ------------------------------------------------------------------------------------------
       IF ( yn_bad_pixel .OR. locerrstat >= pge_errstat_error ) THEN
          errstat = MAX ( errstat, locerrstat )
          cross_track_skippix (ipix) = .TRUE.
          addmsg = ''
          WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
          CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
               modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
               locerrstat )
          DEALLOCATE(fitres_out, weightsum, wvl_good, wvl_bad, spc_good, spc_bad, bad_idx )
          CYCLE
       END IF
       
       ! -----------------------------------------------------------------------
       ! Save crucial variables for across-track reference in Earthshine fitting
       ! -----------------------------------------------------------------------
       ins_sol_wav_avg (ipix) = sol_wav_avg
       solcal_pars (1:max_calfit_idx,ipix) = fitvar_cal(1:max_calfit_idx)
       solcal_xflag(ipix) = INT (solcal_exval, KIND=i2)
       solcal_itnum(ipix) = INT (local_solcal_itnum, KIND=i2)
       
       ! ------------------------------------------------------------------------
       ! Save the processed solar spectrum in its original array. Note that the
       ! spectrum is now normalized, and that the wavelength array is calibrated.
       ! ------------------------------------------------------------------------
       nwav_irrad(ipix)                  = n_irradwvl
       irradiance_wavl(1:n_irradwvl,ipix) = curr_sol_spec(wvl_idx,1:n_irradwvl)
       irradiance_spec(1:n_irradwvl,ipix) = curr_sol_spec(spc_idx,1:n_irradwvl)
       irradiance_wght(1:n_irradwvl,ipix) = curr_sol_spec(sig_idx,1:n_irradwvl)

       addmsg = ''
       WRITE (addmsg, '(A,I2,5(A,1PE10.3),2(A,I5))') 'SOLAR FIT #', ipix, &
            ': hw 1/e = ', hw1e, '; e_asy = ', e_asym,  '; g_sha = ', g_shap, '; shift = ', &
            fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), '; exit val = ', &
            solcal_exval, '; iter num = ', local_solcal_itnum
       CALL error_check ( &
            0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
            vb_lev_omidebug, errstat )
       IF ( pcfvar%verb_thresh_lev >= vb_lev_screen  ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

       IF (ctrvar%yn_diagnostic_run) THEN
          CALL he5_write_solarwavcal (n_irradwvl, ipix, fitvar_cal(shi_idx), fitvar_cal(squ_idx), &
               fitres_out(1:n_irradwvl), locerrstat)
       ENDIF
       DEALLOCATE(fitres_out, weightsum, wvl_good, wvl_bad, spc_good, spc_bad, bad_idx )

    END DO XtrackSolCal
    
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_solar_calibration_loop


  SUBROUTINE solar_fit ( &
       n_fitres_loop, fitres_range, n_sol_wvl,                            &
       curr_sol_spec, hw1e, e_asym, solcal_exval, local_solcal_itnum, chisquav, &
       yn_bad_pixel, fitres_out, errstat )

    ! ***************************************************************
    !
    !   Perform solar wavelength calibration and slit width fitting
    !
    ! ***************************************************************
    USE OMSAO_parameters_module, ONLY: elsunc_less_is_noise, r8_missval
    USE OMSAO_variables_module,  ONLY: ctrvar
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: n_fitres_loop, n_sol_wvl, fitres_range

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: hw1e, e_asym, chisquav
    INTEGER (KIND=i4), INTENT (OUT) :: solcal_exval
    INTEGER (KIND=i2), INTENT (OUT) :: local_solcal_itnum

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL, INTENT (OUT) :: yn_bad_pixel
    REAL (KIND=r8), DIMENSION(n_sol_wvl), INTENT (OUT) :: fitres_out
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    REAL (KIND=r8), DIMENSION(wvl_idx:ccd_idx,1:n_sol_wvl), INTENT (INOUT) :: curr_sol_spec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, i, j, locitnum, n_nozero_wgt
    REAL (KIND=r8) :: rms, mean, sdev, loclim
    REAL (KIND=r8), DIMENSION (n_sol_wvl) :: fitres, fitspec
    REAL (KIND=r8), DIMENSION (max_calfit_idx) :: fitvar
    REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: covar


    EXTERNAL specfit_func_sol
    ! ----------------------------------------------------------------
    ! Initialize local error status variable; note that error handling
    ! is rudimentary in this subroutine - no error is reported.
    ! ----------------------------------------------------------------
    locerrstat = pge_errstat_ok
    solcal_exval = i4_missval
    local_solcal_itnum = i2_missval
    chisquav = r8_missval
    yn_bad_pixel = .FALSE.
    
    ! --------------------------------------------------------------
    ! Calculate and iterate on the irradiance spectrum.
    ! --------------------------------------------------------------
    fitwavs (1:n_sol_wvl) = curr_sol_spec(wvl_idx,1:n_sol_wvl)
    fitweights(1:n_sol_wvl) = curr_sol_spec(sig_idx,1:n_sol_wvl)
    currspec (1:n_sol_wvl) = curr_sol_spec(spc_idx,1:n_sol_wvl)

    ! -------------------------------------------------------------
    ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
    ! set to the initial values in the calling routine. outside the
    ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be 
    ! updated with current values from the previous fit if that fit
    ! has gone well.
    ! -------------------------------------------------------------
    fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)

    ! ---------------------------------------------------------
    ! Assign varied fitting variables to array passed to Elsunc
    ! ---------------------------------------------------------
    fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8
    n_fitvar_cal = 0
    DO i = 1, max_calfit_idx
       IF (ctrvar%lo_sunbnd(i) < ctrvar%up_sunbnd(i) ) THEN
          n_fitvar_cal  = n_fitvar_cal + 1
          mask_fitvar_cal(n_fitvar_cal) = i
          fitvar(n_fitvar_cal) = fitvar_cal(i)
          lobnd(n_fitvar_cal) = ctrvar%lo_sunbnd(i)
          upbnd(n_fitvar_cal) = ctrvar%up_sunbnd(i)
       END IF
    END DO

    ! --------------------------------------------------------
    ! Check whether we have enough spectral points to carry 
    ! out the fitting. If not, call it a bad pixel and return.
    ! --------------------------------------------------------
    IF ( n_fitvar_cal >= n_sol_wvl ) THEN
       yn_bad_pixel = .TRUE.  ;  RETURN
    END IF

    ALLOCATE ( covar(1:n_fitvar_cal,1:n_fitvar_cal) )
    CALL specfit (                                                    &
         n_fitvar_cal, fitvar(1:n_fitvar_cal), n_sol_wvl,             &
         lobnd(1:n_fitvar_cal), upbnd(1:n_fitvar_cal), ctrvar%max_itnum_sol, &
         covar(1:n_fitvar_cal,1:n_fitvar_cal), fitspec(1:n_sol_wvl),  &
         fitres(1:n_sol_wvl), solcal_exval, locitnum, specfit_func_sol )
    IF ( ALLOCATED (covar) ) DEALLOCATE (covar)

    ! ------------------------------------------
    ! Assign iteration number from the first fit
    ! ------------------------------------------
    local_solcal_itnum = INT ( locitnum, KIND=i2 )

    ! ---------------------------------------------------------------------
    ! Attempt to standardize the re-iteration with spectral points excluded
    ! that have fitting residuals larger than a pre-set window. Needs more
    ! thinking before it can replace a simple window determined empirically
    ! from fitting lots of spectra.
    ! ---------------------------------------------------------------------
    n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_sol_wvl)) ) )
    mean = SUM ( fitres(1:n_sol_wvl) ) / REAL(n_nozero_wgt, KIND=r8)
    sdev = SQRT ( SUM ( (fitres(1:n_sol_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
    loclim = REAL (fitres_range, KIND=r8)*sdev

    ! ----------------------
    ! Fitting RMS and CHI**2
    ! ----------------------
    IF ( n_nozero_wgt > 0 ) THEN
       rms = SQRT ( SUM ( fitres(1:n_sol_wvl)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
       ! ---------------------------------------------
       ! This gives the same CHI**2 as the NR routines
       ! ---------------------------------------------
      chisquav = SUM ( fitres(1:n_sol_wvl)**2 )
    ELSE
       rms = r8_missval
       chisquav = r8_missval
    END IF

    ! -----------------------------------------------------------------------
    ! Refit if any part of the fitting residual computed above is larger than
    ! the pre-set window given by FITRES_RANGE. N_FITRES_LOOP must be set > 0
    ! since it determines the maximum number of re-iterations (we don't want
    ! to fit forever!).
    ! -----------------------------------------------------------------------
    IF ( ( n_fitres_loop                    >  0             ) .AND. &
         ( loclim                           >  0.0_r8        ) .AND. &
         ( MAXVAL(ABS(fitres(1:n_sol_wvl))) >= loclim        ) .AND. &
         ( n_nozero_wgt                     >  n_fitvar_cal  )  ) THEN

       fitloop: DO j = 1, n_fitres_loop
          WHERE ( ABS(fitres(1:n_sol_wvl)) > loclim )
             fitweights(1:n_sol_wvl) = downweight
          END WHERE

          ALLOCATE ( covar(1:n_fitvar_cal,1:n_fitvar_cal) )
          CALL specfit ( &
               n_fitvar_cal, fitvar(1:n_fitvar_cal), n_sol_wvl, &
               lobnd(1:n_fitvar_cal), upbnd(1:n_fitvar_cal), ctrvar%max_itnum_sol, &
               covar(1:n_fitvar_cal,1:n_fitvar_cal), fitspec(1:n_sol_wvl), &
               fitres(1:n_sol_wvl), solcal_exval, locitnum, specfit_func_sol )
          IF ( ALLOCATED (covar) ) DEALLOCATE (covar)

          IF ( solcal_exval > 0 ) THEN
             fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
          ELSE
             fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_sol_init(1:max_calfit_idx)
          END IF

          ! ----------------------
          ! Fitting RMS and CHI**2
          ! ----------------------
          n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_sol_wvl)) ) )
          IF ( n_nozero_wgt > 0.0_r8 ) THEN
             rms     = SQRT ( SUM ( fitres(1:n_sol_wvl)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
             ! ---------------------------------------------
             ! This gives the same CHI**2 as the NR routines
             ! ---------------------------------------------
             chisquav = SUM  ( fitres(1:n_sol_wvl)**2 )
          ELSE
             rms      = r8_missval
             chisquav = r8_missval
          END IF

          ! -----------------------------
          ! Add any subsequent iterations
          ! -----------------------------
          local_solcal_itnum = local_solcal_itnum + INT ( locitnum, KIND=i2 )

          ! --------------------------------------------------------------
          ! Exit iteration loop if fitting residual is within contstraints
          ! --------------------------------------------------------------
          IF ( MAXVAL(ABS(fitres(1:n_sol_wvl))) <= loclim ) EXIT fitloop

       END DO fitloop
    END IF

    ! ---------------------------------------------------------------
    ! The following assignment makes sense only because FITVAR_CAL is
    ! updated with FITVAR (using the proper mask) in SPECTRUM_SOLAR.
    ! ---------------------------------------------------------------
    IF ( solcal_exval >= INT(elsunc_less_is_noise, KIND=i4) ) THEN
       fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
    ELSE
       fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_sol_init(1:max_calfit_idx)
    END IF
    
    ! ---------------------------------------------------------------
    ! Save shifted&squeezed wavelength array, and the fitting weights
    ! ---------------------------------------------------------------
    curr_sol_spec(wvl_idx,1:n_sol_wvl) = (fitwavs (1:n_sol_wvl) - fitvar_cal_saved(shi_idx) + &
         sol_wav_avg * fitvar_cal_saved(squ_idx)) / (1.0_r8 + fitvar_cal_saved(squ_idx))
    curr_sol_spec(sig_idx,1:n_sol_wvl) = fitweights (1:n_sol_wvl)

    ! --------------------------------------------
    ! Save the slit function parameters for output
    ! --------------------------------------------
    hw1e   = fitvar_cal(hwe_idx) ; e_asym = fitvar_cal(asy_idx) ; g_shap = fitvar_cal(sha_idx)
    fitres_out(1:n_sol_wvl) = fitres(1:n_sol_wvl)
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE solar_fit

END MODULE OMSAO_solar_wavcal_module
