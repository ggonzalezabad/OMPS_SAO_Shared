MODULE OMSAO_solar_wavcal_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: &
       wvl_idx, sig_idx, spc_idx, ccd_idx, max_calfit_idx, shi_idx, squ_idx, solcal_idx,  &
       hwe_idx, asy_idx, solar_idx, bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, bl5_idx, &
       sc0_idx, sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx, sin_idx
  USE OMSAO_parameters_module, ONLY: &
       i2_missval, i4_missval, r8_missval, normweight, downweight, maxchlen, max_spec_pts
  USE OMSAO_variables_module,  ONLY: &
       hw1e, e_asym, curr_sol_spec, sol_wav_avg, fitvar_cal, fitvar_cal_saved,  &
       mask_fitvar_cal, n_fitvar_cal, lobnd, upbnd,&
       fitwavs, fitweights, currspec, refspecs_original,    &
       pcfvar, ctrvar
  USE OMSAO_omidata_module
  USE OMSAO_errstat_module
  USE OMSAO_he5_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE xtrack_solar_calibration_loop ( first_pix, last_pix, errstat )

    USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze

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
    INTEGER   (KIND=i2)                     :: solcal_itnum
    INTEGER   (KIND=i4)                     :: locerrstat, ipix, solcal_exval, n_sol_wvl
    CHARACTER (LEN=maxchlen)                :: addmsg
    REAL      (KIND=r8)                     :: chisquav
    LOGICAL                                 :: yn_skip_pix, yn_bad_pixel
    INTEGER (KIND=i4), DIMENSION (4)        :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)        :: exclud_idx
    REAL   (KIND=r8), DIMENSION(nwavel_max) :: fitres_out

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=29), PARAMETER :: modulename = 'xtrack_solar_calibration_loop'

    omi_solcal_chisq = r8_missval

    fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_sol_init(1:max_calfit_idx)

    ! ---------------------------------------------------------------
    ! Loop for solar wavelength calibration and slit function fitting
    ! ---------------------------------------------------------------
    XtrackSolCal: DO ipix = first_pix, last_pix

       locerrstat = pge_errstat_ok

       curr_xtrack_pixnum = ipix

       n_omi_irradwvl = omi_nwav_irrad(ipix)

       IF ( n_omi_irradwvl <= 0 ) CYCLE

       saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

       ! -------------------------------------------------------------------------
       select_idx(1:4) = omi_ccdpix_selection(ipix,1:4)
       exclud_idx(1:2) = omi_ccdpix_exclusion(ipix,1:2)
       CALL omi_adjust_irradiance_data ( &           ! Set up generic fitting arrays
            select_idx(1:4), exclud_idx(1:2),             &
            n_omi_irradwvl,                               &
            omi_irradiance_wavl  (1:n_omi_irradwvl,ipix), &
            omi_irradiance_spec  (1:n_omi_irradwvl,ipix), &
            omi_irradiance_ccdpix(1:n_omi_irradwvl,ipix), &
            n_sol_wvl, curr_sol_spec(wvl_idx:ccd_idx,1:n_omi_irradwvl), &
            yn_skip_pix, locerrstat )

       ! -------------------------------------------------------------------------
       IF ( yn_skip_pix .OR. locerrstat >= pge_errstat_error ) THEN
          errstat = MAX ( errstat, locerrstat )
          omi_cross_track_skippix (ipix) = .TRUE.
          addmsg = ''
          WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
          CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
               modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
               locerrstat )
          CYCLE
       END IF

       yn_bad_pixel   = .FALSE.
       CALL solar_fit ( &                              ! Solar wavelength calibration
            ctrvar%n_fitres_loop(solcal_idx), ctrvar%fitres_range(solcal_idx), n_sol_wvl,                       &
            curr_sol_spec(wvl_idx:ccd_idx,1:n_sol_wvl), hw1e, e_asym, solcal_exval, solcal_itnum, &
            chisquav, yn_bad_pixel, fitres_out(1:n_sol_wvl), locerrstat )
       ! ------------------------------------------------------------------------------------------
       IF ( yn_bad_pixel .OR. locerrstat >= pge_errstat_error ) THEN
          errstat = MAX ( errstat, locerrstat )
          omi_cross_track_skippix (ipix) = .TRUE.
          addmsg = ''
          WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
          CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
               modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
               locerrstat )
          CYCLE
       END IF

       ! -----------------------------------------------------------------------
       ! Save crucial variables for across-track reference in Earthshine fitting
       ! -----------------------------------------------------------------------
       omi_sol_wav_avg (ipix)                     = sol_wav_avg
       omi_solcal_chisq(ipix)                     = chisquav
       omi_solcal_pars (1:max_calfit_idx,ipix)    = fitvar_cal(1:max_calfit_idx)
       omi_solcal_xflag(ipix)                     = INT (solcal_exval, KIND=i2)
       omi_solcal_itnum(ipix)                     = INT (solcal_itnum, KIND=i2)
       
       ! ------------------------------------------------------------------------
       ! Save the processed solar spectrum in its original array. Note that the
       ! spectrum is now normalized, has bad pixels set to -1, and that the
       ! wavelength array is calibrated.
       ! ------------------------------------------------------------------------
       omi_nwav_irrad(ipix)                  = n_sol_wvl
       omi_irradiance_wavl(1:n_sol_wvl,ipix) = curr_sol_spec(wvl_idx,1:n_sol_wvl)
       omi_irradiance_spec(1:n_sol_wvl,ipix) = curr_sol_spec(spc_idx,1:n_sol_wvl)
       omi_irradiance_wght(1:n_sol_wvl,ipix) = curr_sol_spec(sig_idx,1:n_sol_wvl)

       addmsg = ''
       WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'SOLAR FIT          #', ipix, &
            ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; shift = ', &
            fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), '; exit val = ', &
            solcal_exval, '; iter num = ', solcal_itnum
       CALL error_check ( &
            0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
            vb_lev_omidebug, errstat )
       IF ( pcfvar%verb_thresh_lev >= vb_lev_screen  ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

       CALL he5_write_solarwavcal (n_sol_wvl, ipix, fitvar_cal(shi_idx), fitres_out(1:n_sol_wvl), locerrstat)

    END DO XtrackSolCal

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_solar_calibration_loop


  SUBROUTINE solar_fit ( &
       n_fitres_loop, fitres_range, n_sol_wvl,                            &
       curr_sol_spec, hw1e, e_asym, solcal_exval, solcal_itnum, chisquav, &
       yn_bad_pixel, fitres_out, errstat )

    ! ***************************************************************
    !
    !   Perform solar wavelength calibration and slit width fitting
    !
    ! ***************************************************************
    USE OMSAO_parameters_module, ONLY: elsunc_less_is_noise, r8_missval
    USE OMSAO_variables_module,  ONLY: yn_newshift
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: n_fitres_loop, n_sol_wvl, fitres_range

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), INTENT (OUT)   :: hw1e, e_asym, chisquav
    INTEGER (KIND=i4), INTENT (OUT)   :: solcal_exval
    INTEGER (KIND=i2), INTENT (OUT)   :: solcal_itnum

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,                                                   INTENT (OUT)   :: yn_bad_pixel
    REAL    (KIND=r8), DIMENSION(n_sol_wvl),                   INTENT (OUT)   :: fitres_out
    INTEGER (KIND=i4),                                         INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION(wvl_idx:ccd_idx,1:n_sol_wvl), INTENT (INOUT) :: curr_sol_spec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locerrstat, i, j, locitnum, n_nozero_wgt
    REAL    (KIND=r8)  :: rms, mean, mdev, sdev, loclim
    REAL    (KIND=r8), DIMENSION (n_sol_wvl)         :: fitres, fitspec
    REAL    (KIND=r8), DIMENSION (max_calfit_idx)    :: fitvar
    REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE  :: covar


    EXTERNAL specfit_func_sol

    ! ----------------------------------------------------------------
    ! Initialize local error status variable; note that error handling
    ! is rudimentary in this subroutine - no error is reported.
    ! ----------------------------------------------------------------
    locerrstat = pge_errstat_ok

    solcal_exval = i4_missval
    solcal_itnum = i2_missval

    chisquav = r8_missval

    yn_bad_pixel = .FALSE.

    ! --------------------------------------------------------------
    ! Calculate and iterate on the irradiance spectrum.
    ! --------------------------------------------------------------
    fitwavs   (1:n_sol_wvl) = curr_sol_spec(wvl_idx,1:n_sol_wvl)
    fitweights(1:n_sol_wvl) = curr_sol_spec(sig_idx,1:n_sol_wvl)
    currspec  (1:n_sol_wvl) = curr_sol_spec(spc_idx,1:n_sol_wvl)

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
          lobnd (n_fitvar_cal) = ctrvar%lo_sunbnd(i)
          upbnd (n_fitvar_cal) = ctrvar%up_sunbnd(i)
       END IF
    END DO

    ! --------------------------------------------------------------------
    ! Check whether we enough spectral points to carry out the fitting. If
    ! not, call it a bad pixel and return.
    ! --------------------------------------------------------------------
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
    solcal_itnum = INT ( locitnum, KIND=i2 )

    ! ---------------------------------------------------------------------
    ! Attempt to standardize the re-iteration with spectral points excluded
    ! that have fitting residuals larger than a pre-set window. Needs more
    ! thinking before it can replace a simple window determined empirically
    ! from fitting lots of spectra.
    ! ---------------------------------------------------------------------
     n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_sol_wvl)) ) )
     mean         = SUM  ( fitres(1:n_sol_wvl) )                 / REAL(n_nozero_wgt,   KIND=r8)
     sdev         = SQRT ( SUM ( (fitres(1:n_sol_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
     mdev         = SUM  ( ABS(fitres(1:n_sol_wvl)-mean) )       / REAL(n_nozero_wgt,   KIND=r8)
     loclim       = REAL (fitres_range, KIND=r8)*sdev

    ! ----------------------
    ! Fitting RMS and CHI**2
    ! ----------------------
    IF ( n_nozero_wgt > 0 ) THEN
       rms     = SQRT ( SUM ( fitres(1:n_sol_wvl)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
       ! ---------------------------------------------
       ! This gives the same CHI**2 as the NR routines
       ! ---------------------------------------------
      chisquav = SUM  ( fitres(1:n_sol_wvl)**2 )
    ELSE
       rms      = r8_missval
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
          solcal_itnum = solcal_itnum + INT ( locitnum, KIND=i2 )

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
    IF (yn_newshift .EQV. .true.) THEN !gga
       curr_sol_spec(wvl_idx,1:n_sol_wvl) = &
            (fitwavs (1:n_sol_wvl) - fitvar_cal_saved(shi_idx) + &
             sol_wav_avg * fitvar_cal_saved(squ_idx)) /          &
             (1.0_r8 + fitvar_cal_saved(squ_idx))
    ELSE !gga
       curr_sol_spec(wvl_idx,1:n_sol_wvl) = fitwavs (1:n_sol_wvl)
    END IF
    curr_sol_spec(sig_idx,1:n_sol_wvl) = fitweights (1:n_sol_wvl)

    ! ------------------------------------------------
    !  Save the slit function parameters for later use 
    ! in the undersampling correction.
    ! ------------------------------------------------
    hw1e   = fitvar_cal(hwe_idx)  ;  e_asym = fitvar_cal(asy_idx)
    fitres_out(1:n_sol_wvl) = fitres(1:n_sol_wvl)
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE solar_fit

  SUBROUTINE xtrack_inverse_solar_calibration_loop( first_pix, last_pix, errstat )

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
    INTEGER (KIND=i4) :: locerrstat, ipix, ivar
    INTEGER (KIND=i4) :: nwav
    REAL    (KIND=r8) :: sol_wav_avg
    REAL    (KIND=r8), DIMENSION(nwavel_max) :: tmpwav, tmpspe, modspe, &
                                                del, delx

    locerrstat = 0_i4

    DO ipix = first_pix, last_pix

       nwav        = omi_nwav_irrad(ipix)
       sol_wav_avg = omi_sol_wav_avg(ipix)
       
       tmpwav(1:nwav) = omi_irradiance_wavl(1:nwav,ipix)
       modspe(1:nwav) = omi_irradiance_spec(1:nwav,ipix)
       del(1:nwav)    = tmpwav(1:nwav) - sol_wav_avg

       ! Add baseline parameters
       tmpspe(1:nwav) = 0.0_r8
       delx(1:nwav)   = 1.0_r8
       DO ivar = bl0_idx, bl5_idx
          IF (omi_solcal_pars(ivar,ipix) .NE. 0.0_r4) THEN
             tmpspe(1:nwav) = tmpspe(1:nwav) + omi_solcal_pars(ivar,ipix) * delx(1:nwav)
          END IF
          delx(1:nwav) = delx(1:nwav) * del(1:nwav)
       END DO
       modspe(1:nwav) = modspe(1:nwav) - tmpspe

       ! Add scaling parameters
       tmpspe(1:nwav) = 0.0_r8
       delx(1:nwav)   = 1.0_r8
       DO ivar = sc0_idx, sc5_idx
          IF (omi_solcal_pars(ivar,ipix) .NE. 0.0_r4) THEN
             tmpspe(1:nwav) = tmpspe(1:nwav) + omi_solcal_pars(ivar,ipix) * delx(1:nwav)
          END IF
          delx(1:nwav) = delx(1:nwav) * del(1:nwav)
       END DO
       modspe(1:nwav) = modspe(1:nwav) / tmpspe

       ! Add "albedo"
       modspe(1:nwav) = modspe(1:nwav) / omi_solcal_pars(sin_idx,ipix)

       ! Saving the inverse calibrated spectra in the omi_irradiance_spec
       omi_irradiance_spec(1:nwav, ipix) = modspe(1:nwav)

  

    END DO

    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE xtrack_inverse_solar_calibration_loop

END MODULE OMSAO_solar_wavcal_module
