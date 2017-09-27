SUBROUTINE xtrack_radiance_wvl_calibration (             &
     yn_radiance_reference, yn_solar_comp,               &
     first_pix, last_pix, n_max_rspec, n_comm_wvl_out, errstat )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, &
       max_calfit_idx, max_rs_idx, hwe_idx, asy_idx,  &
       shi_idx, squ_idx, solar_idx, ccd_idx, radcal_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, downweight, normweight
  USE OMSAO_variables_module, ONLY: hw1e, e_asym, &
       n_rad_wvl, curr_rad_spec, sol_wav_avg, database, fitvar_cal, &
       fitvar_cal_saved, pcfvar, ctrvar
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_data_module, ONLY: nwavel_max, nxtrack_max, &
       cross_track_skippix, nwav_radref, radcal_itnum, &
       radcal_xflag, radcal_chisq, n_ins_database_wvl, &
       solcal_pars, omi_sol_wav_avg,  &
       nwav_irrad, irradiance_wght, irradiance_wavl, &
       irradiance_spec, ins_database, ins_database_wvl, &
       radref_spec, radref_wavl, radref_qflg, radref_wght, &
       nwav_rad, radiance_spec, radiance_wavl, radiance_qflg, &
       ccdpix_selection, ccdpix_exclusion, radiance_ccdpix, &
       radcal_pars, curr_xtrack_pixnum, n_irradwvl, n_radwvl      
  USE OMSAO_errstat_module, ONLY: f_sep, omsao_s_progress, omsao_w_skippix, &
       pge_errstat_error,pge_errstat_ok, pge_errstat_warning, vb_lev_default, &
       vb_lev_omidebug, vb_lev_screen, error_check
  USE OMSAO_he5_module, ONLY : NrofScanLines
  USE EZspline_obj
  USE EZspline

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix, n_max_rspec
  LOGICAL,           INTENT (IN) :: yn_radiance_reference, yn_solar_comp

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
  INTEGER   (KIND=i2)      :: local_radcal_itnum
  INTEGER   (KIND=i4)      :: locerrstat, ipix, radcal_exval, i, imax, n_ref_wvl, cline
  REAL      (KIND=r8)      :: chisquav, rad_spec_avg
  LOGICAL                  :: yn_skip_pix, yn_bad_pixel, yn_full_range
  CHARACTER (LEN=maxchlen) :: addmsg
  REAL   (KIND=r8),  DIMENSION(nwavel_max) :: fitres_out
  REAL   (KIND=r8),  DIMENSION(nwavel_max) :: calibration_wavl
  REAL   (KIND=r8),  DIMENSION(nwavel_max) :: calibration_spec
  REAL   (KIND=r8),  DIMENSION(nwavel_max) :: calibration_qflg
  INTEGER (KIND=i4), DIMENSION (4)            :: select_idx
  INTEGER (KIND=i4), DIMENSION (2)            :: exclud_idx
  REAL    (KIND=r8), DIMENSION (n_max_rspec) :: ref_wvl, ref_spc, ref_wgt, rad_wvl

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=31), PARAMETER :: modulename = 'xtrack_radiance_wvl_calibration'


  locerrstat = pge_errstat_ok

  fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_rad_init(1:max_calfit_idx)

  ! -------------------------------------------------
  ! Set the number of wavelengths for the common mode
  ! -------------------------------------------------
  n_comm_wvl_out = MAXVAL ( nwav_radref(first_pix:last_pix) )
  IF ( MAXVAL(nwav_rad(first_pix:last_pix,0)) > n_comm_wvl_out ) &
  n_comm_wvl_out = MAXVAL(nwav_rad(first_pix:last_pix,0))

  ! ---------------------------------------------------
  ! Find a line close or nearby the center of the swath
  ! ---------------------------------------------------
  cline = INT(NrofScanLines / 2, KIND = i4)

  ! --------------------------------
  ! Loop over cross-track positions. 
  ! --------------------------------
  XTrackWavCal: DO ipix = first_pix, last_pix

     locerrstat = pge_errstat_ok

     curr_xtrack_pixnum = ipix

     ! ---------------------------------------------------------------------
     ! If we already determined that this cross track pixel position carries
     ! an error, we don't even have to start processing.
     ! ---------------------------------------------------------------------
     IF ( cross_track_skippix(ipix) ) CYCLE

     ! ---------------------------------------------------------------------------
     ! For each cross-track position we have to initialize the saved Shift&Squeeze
     ! ---------------------------------------------------------------------------
     saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

     ! ----------------------------------------------------
     ! Assign number of radiance and irradiance wavelengths
     ! ----------------------------------------------------
     n_irradwvl = nwav_irrad(ipix  )
     n_radwvl   = nwav_rad  (ipix,cline)

     ! -----------------------------------------------------------------
     ! tpk: Should the following be "> n_fitvar_rad"??? No, because that
     !      value is set only inside OMI_ADJUST_RADIANCE_DATA!!!
     ! -----------------------------------------------------------------
     IF ( n_irradwvl <= 0 .OR. n_radwvl <= 0 ) CYCLE

     ! ---------------------------------------------------------------
     ! Restore solar fitting variables for across-track reference in
     ! Earthshine fitting. Use the Radiance References if appropriate.
     ! ---------------------------------------------------------------
     sol_wav_avg = omi_sol_wav_avg(ipix)
     hw1e        = solcal_pars(hwe_idx,ipix)
     e_asym      = solcal_pars(asy_idx,ipix)

     ! -----------------------------------------------------
     ! Assign (hopefully predetermined) "reference" weights.
     ! -----------------------------------------------------
     IF ( .NOT. yn_solar_comp ) THEN
        n_irradwvl            = nwav_irrad(ipix)
        ref_wgt(1:n_irradwvl) = irradiance_wght(1:n_irradwvl,ipix)

       ! -----------------------------------------------------
       ! Catch the possibility that N_OMI_RADWVL > N_OMI_IRRADWVL
       ! -----------------------------------------------------
       IF ( n_radwvl > n_irradwvl ) THEN
          i = n_radwvl - n_irradwvl
          ref_wgt(n_irradwvl+1:n_irradwvl+i) = downweight
          n_irradwvl = n_radwvl
       END IF
     ELSE
        n_irradwvl          = n_radwvl
        ref_wgt(1:n_radwvl) = normweight
     END IF

     ! -----------------------------------------------------------------
     ! If a Radiance Reference is being used, then it must be calibrated
     ! rather than the swath line that has been read, otherwise we copy
     ! the a line close to the center of the swath into position 0 and
     ! performe the calibration with it.
     ! -----------------------------------------------------------------
     IF ( yn_radiance_reference ) THEN
        calibration_wavl(1:n_radwvl) = radref_wavl(1:n_radwvl,ipix)
        calibration_spec(1:n_radwvl) = radref_spec(1:n_radwvl,ipix)
        calibration_qflg(1:n_radwvl) = radref_qflg(1:n_radwvl,ipix)
     ELSE
        calibration_wavl(1:n_radwvl) = radiance_wavl(1:n_radwvl,ipix,cline)
        calibration_spec(1:n_radwvl) = radiance_spec(1:n_radwvl,ipix,cline)
        calibration_qflg(1:n_radwvl) = radiance_qflg(1:n_radwvl,ipix,cline)
     END IF

     ! ---------------------------------------------------------------------------
     ! Set up generic fitting arrays. Remember that OMI_RADIANCE_XXX arrays are
     ! 3-dim with the last dimension being the scan line numbers. For the radiance
     ! wavelength calibration we only have one scan line at index "0".
     ! ---------------------------------------------------------------------------
     select_idx(1:4) = ccdpix_selection(ipix,1:4)
     exclud_idx(1:2) = ccdpix_exclusion(ipix,1:2)

     CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2),            &
          n_radwvl,                                &
          calibration_wavl   (1:n_radwvl),         &
          calibration_spec   (1:n_radwvl),         &
          radiance_ccdpix(1:n_radwvl,ipix,0),  &
          n_irradwvl, ref_wgt(1:n_irradwvl),   &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_radwvl), rad_spec_avg, &
          yn_skip_pix )
     
     ! ------------------------------------------------------------------------------------
     IF ( yn_skip_pix .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
             modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
             locerrstat )
        CYCLE
     END IF

     ! -----------------------------------------------------
     ! Assign the solar average wavelength - the wavelength
     ! calibration will not converge without it!
     ! -----------------------------------------------------
     sol_wav_avg = &
          SUM ( curr_rad_spec(wvl_idx,1:n_radwvl) ) / REAL(n_radwvl,KIND=r8)
     yn_bad_pixel = .FALSE.
     CALL radiance_wavcal ( &                       ! Radiance wavelength calibration
          ipix, ctrvar%n_fitres_loop(radcal_idx), ctrvar%fitres_range(radcal_idx),       &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),           &
          radcal_exval, local_radcal_itnum, chisquav, yn_bad_pixel, fitres_out(1:n_rad_wvl), locerrstat )

     IF ( yn_bad_pixel .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
             modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
             locerrstat )
        CYCLE
     END IF
     ! ------------------------------------------------------------------------------------
           
     addmsg = ''
     WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'RADIANCE Wavcal    #', ipix, &
          ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; shift = ', &
          fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), &
         '; exit val = ', radcal_exval, '; iter num = ', local_radcal_itnum
     CALL error_check ( &
          0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
          vb_lev_omidebug, locerrstat )
     IF ( pcfvar%verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

     ! ---------------------------------
     ! Save crucial variables for output
     ! ---------------------------------
     radcal_pars (1:max_calfit_idx,ipix) = fitvar_cal(1:max_calfit_idx)
     radcal_xflag(ipix)                  = INT (radcal_exval, KIND=i2)
     radcal_itnum(ipix)                  = INT (local_radcal_itnum, KIND=i2)
     radcal_chisq(ipix)                  = chisquav

     ! -----------------------------------------------------------------------

     IF ( .NOT. (yn_radiance_reference) ) THEN
        n_ref_wvl = n_irradwvl
        ref_wvl(1:n_ref_wvl) = irradiance_wavl(1:n_ref_wvl,ipix)
        ref_spc(1:n_ref_wvl) = irradiance_spec(1:n_ref_wvl,ipix)
        ref_wgt(1:n_ref_wvl) = irradiance_wght(1:n_ref_wvl,ipix)
     ELSE
        n_ref_wvl = n_rad_wvl
        ref_wvl(1:n_ref_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
        ref_spc(1:n_ref_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl)
        ref_wgt(1:n_ref_wvl) = curr_rad_spec(sig_idx,1:n_rad_wvl)
     END IF
    
     nwav_radref(ipix)                        = n_ref_wvl
     radref_wavl(1:n_ref_wvl,ipix)            = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     radref_spec(1:n_ref_wvl,ipix)            = curr_rad_spec(spc_idx,1:n_rad_wvl)
     radref_wght(1:n_rad_wvl,ipix)            = curr_rad_spec(sig_idx,1:n_rad_wvl)
     radref_wght(n_rad_wvl+1:nwavel_max,ipix) = downweight

     CALL he5_write_radiancewavcal ( n_rad_wvl, ipix, fitvar_cal(shi_idx), fitres_out(1:n_rad_wvl), locerrstat)

     ! ----------------------------------------------------
     ! Spline reference spectra to current wavelength grid.
     ! ----------------------------------------------------
     rad_wvl(1:n_rad_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     Call prepare_databases ( &
          ipix, n_ref_wvl, ref_wvl(1:n_ref_wvl), ref_spc(1:n_ref_wvl), &
          n_rad_wvl, rad_wvl(1:n_rad_wvl), n_max_rspec, locerrstat )

     IF ( locerrstat >= pge_errstat_error ) EXIT XTrackWavCal

     ! ---------------------------------------------------------
     ! Save DATABASE in OMI_DATABASE for radiance fitting loops.
     ! ---------------------------------------------------------
     ins_database      (1:max_rs_idx,1:n_rad_wvl,ipix) = database     (1:max_rs_idx,1:n_rad_wvl)
     ins_database_wvl  (             1:n_rad_wvl,ipix) = curr_rad_spec(     wvl_idx,1:n_rad_wvl)
     n_ins_database_wvl(                         ipix) = n_rad_wvl
     
     ! ----------------------------------------------------------------------
     ! Update the radiance reference with the wavelength calibrated values.
     ! ----------------------------------------------------------------------
     IF ( yn_radiance_reference ) THEN

        ! --------------------------------------------------------
        ! Update the solar spectrum entry in OMI_DATABASE. First
        ! re-sample the solar reference spectrum to the OMI grid
        ! then assign to data base.
        !
        ! We need to keep the irradiance spectrum because we still
        ! have to fit the radiance reference, and we can't really
        ! do that against itself. In a later module the irradiance
        ! is replaced by the radiance reference.
        ! --------------------------------------------------------

        ! ------------------------------------------------------------------
        ! Prevent failure of interpolation by finding the maximum wavelength
        ! of the irradiance wavelength array.
        ! ------------------------------------------------------------------
        imax = MAXVAL ( MAXLOC ( irradiance_wavl(1:n_irradwvl,ipix) ) )

        CALL interpolation ( &
             imax, irradiance_wavl(1:imax,ipix),                     &
             irradiance_spec(1:imax,ipix),                           &
             n_rad_wvl, ins_database_wvl(1:n_rad_wvl,ipix),              &
             ins_database(solar_idx,1:n_rad_wvl,ipix),                   &
             'endpoints', 0.0_r8, yn_full_range, locerrstat )

        IF ( locerrstat >= pge_errstat_error ) THEN
           errstat = MAX ( errstat, locerrstat )
           cross_track_skippix (ipix) = .TRUE.
           addmsg = ''
           WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
           CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
                modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
                locerrstat )
           CYCLE
        END IF

     END IF

  END DO XTrackWavCal

  ! CCM Write splined/convolved databases if necessary
  IF( ctrvar%yn_diagnostic_run ) THEN
     
     ! ins_database maybe ins_database_wvl?
     CALL he5_write_ins_database(ins_database(1:max_rs_idx,1:n_rad_wvl,1:nxtrack_max), &
          ins_database_wvl(1:n_rad_wvl, 1:nxtrack_max), &
          max_rs_idx, n_rad_wvl, nxtrack_max, errstat)
     
  ENDIF

  errstat = MAX ( errstat, locerrstat )


  RETURN
END SUBROUTINE xtrack_radiance_wvl_calibration


SUBROUTINE xtrack_radiance_fitting_loop (                             &
     n_max_rspec, first_pix, last_pix, pge_idx, iloop,                &
     n_fitvar_rad, allfit_cols, allfit_errs, corr_matrix, &
     target_var, errstat, fitspc_out, fitspc_out_dim0                 )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, &
       o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, &
       pge_o3_idx, solar_idx, ccd_idx, radfit_idx
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval
  USE OMSAO_variables_module,  ONLY: database, curr_sol_spec, n_rad_wvl, &
       curr_rad_spec, sol_wav_avg, hw1e, e_asym, n_database_wvl, ctrvar
  USE OMSAO_radiance_ref_module, ONLY: yn_reference_fit
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_data_module, ONLY: nxtrack_max, n_comm_wvl, &
       column_uncert, column_amount, fit_rms, radfit_chisq, &
       itnum_flag, fitconv_flag, solcal_pars, omi_sol_wav_avg, &
       n_ins_database_wvl, nwav_rad, szenith, xtrackpix_no, &
       cross_track_skippix, n_radwvl, n_irradwvl, &
       curr_xtrack_pixnum, o3_uncert, o3_amount, radiance_wavl, &
       ccdpix_exclusion, ccdpix_selection, ins_database, &
       ins_database_wvl, max_rs_idx, radiance_spec, radiance_ccdpix, &
       radref_wght
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok
     
  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: &
       pge_idx, iloop, first_pix, last_pix, n_max_rspec, n_fitvar_rad, &
       fitspc_out_dim0

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  REAL    (KIND=r8), INTENT (OUT  ), DIMENSION (n_fitvar_rad,first_pix:last_pix) :: &
       allfit_cols, allfit_errs, corr_matrix

  ! ---------------------------------------------------------
  ! Optional output variable (fitted variable for target gas)
  ! ---------------------------------------------------------
  REAL (KIND=r8), DIMENSION(ctrvar%n_fincol_idx,first_pix:last_pix), INTENT (OUT) :: target_var

  ! CCM Output fit spectra
  REAL (KIND=r8), DIMENSION(fitspc_out_dim0,nxtrack_max,4), INTENT (OUT) :: fitspc_out

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum
  REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg  
  REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
  LOGICAL                                     :: yn_skip_pix, yn_cycle_this_pix
  LOGICAL                                     :: yn_bad_pixel
  INTEGER (KIND=i4), DIMENSION (4)            :: select_idx
  INTEGER (KIND=i4), DIMENSION (2)            :: exclud_idx
  INTEGER (KIND=i4)                           :: n_solar_pts
  REAL    (KIND=r8), DIMENSION (n_max_rspec)  :: solar_wvl

  ! CCM Array for holding fitted spectra
  REAL    (KIND=r8), DIMENSION (n_comm_wvl)   :: fitspc

  locerrstat = pge_errstat_ok

  XTrackPix: DO ipix = first_pix, last_pix

     curr_xtrack_pixnum = ipix

     ! ---------------------------------------------------------------------
     ! If we already determined that this cross track pixel position carries
     ! an error, we don't even have to start processing.
     ! ---------------------------------------------------------------------
     IF ( cross_track_skippix(ipix) .OR. ctrvar%szamax < szenith(ipix,iloop) ) CYCLE
    
     locerrstat = pge_errstat_ok

     n_database_wvl = n_ins_database_wvl(ipix)
     n_radwvl   = nwav_rad      (ipix,iloop)

     ! ---------------------------------------------------------------------------
     ! For each cross-track position we have to initialize the saved Shift&Squeeze
     ! ---------------------------------------------------------------------------
     saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

     ! ----------------------------------------------------------------------------
     ! Assign the solar wavelengths. Those should be current in the DATABASE array
     ! and can be taken from there no matter which case - YN_SOLAR_COMP and/or
     ! YN_RADIANCE_REFRENCE we are processing.
     ! ----------------------------------------------------------------------------
     n_solar_pts              = n_ins_database_wvl(ipix)
     if (n_solar_pts < 1) cycle  ! JED fix

     solar_wvl(1:n_solar_pts) = ins_database_wvl  (1:n_solar_pts, ipix)
     n_irradwvl           = n_solar_pts

     CALL check_wavelength_overlap ( &
          n_fitvar_rad,                                                &
          n_solar_pts,          solar_wvl (1:n_solar_pts),             &
          n_radwvl, radiance_wavl (1:n_radwvl,ipix,iloop), &
          yn_cycle_this_pix )

     IF ( yn_cycle_this_pix                .OR. &
          (n_database_wvl <= 0) .OR. (n_radwvl <= 0) ) CYCLE
          !(n_database_wvl <= n_fitvar_rad) .OR. (n_radwvl <= n_fitvar_rad) ) CYCLE

     ! ----------------------------------------------
     ! Restore DATABASE from OMI_DATABASE (see above)
     ! ----------------------------------------------
     database (1:max_rs_idx,1:n_database_wvl) = ins_database (1:max_rs_idx,1:n_database_wvl,ipix)
                 
     ! ---------------------------------------------------------------------------------
     ! Restore solar fitting variables for across-track reference in Earthshine fitting.
     ! Note that, for the YN_SOLAR_COMP case, some variables have been assigned already
     ! in the XTRACK_RADIANCE_WAVCAL loop.
     ! ---------------------------------------------------------------------------------
     sol_wav_avg                             = omi_sol_wav_avg(ipix)
     hw1e                                    = solcal_pars(hwe_idx,ipix)
     e_asym                                  = solcal_pars(asy_idx,ipix)
     curr_sol_spec(wvl_idx,1:n_database_wvl) = ins_database_wvl(1:n_database_wvl,ipix)
     curr_sol_spec(spc_idx,1:n_database_wvl) = ins_database    (solar_idx,1:n_database_wvl,ipix)
     ! --------------------------------------------------------------------------------

     xtrackpix_no = ipix

     ! -------------------------------------------------------------------------
     select_idx(1:4) = ccdpix_selection(ipix,1:4)
     exclud_idx(1:2) = ccdpix_exclusion(ipix,1:2)

     CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2),                        &
          n_radwvl,                                            &
          radiance_wavl  (1:n_radwvl,ipix,iloop),          &
          radiance_spec  (1:n_radwvl,ipix,iloop),          &
          radiance_ccdpix(1:n_radwvl,ipix,iloop),          &
          n_radwvl, radref_wght(1:n_radwvl,ipix),      &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_radwvl),&
          rad_spec_avg, yn_skip_pix )

     ! --------------------
     ! The radiance fitting
     ! --------------------
     fitcol       = r8_missval
     dfitcol      = r8_missval
     radfit_exval = INT(i2_missval, KIND=i4)
     radfit_itnum = INT(i2_missval, KIND=i4)
     rms          = r8_missval

     yn_reference_fit = .FALSE.
     IF ( MAXVAL(curr_rad_spec(spc_idx,1:n_rad_wvl)) > 0.0_r8 .AND.     &
          n_rad_wvl > n_fitvar_rad .AND. (.NOT. yn_skip_pix)          ) THEN

        yn_bad_pixel = .FALSE.

        CALL radiance_fit ( &
             pge_idx, ipix, ctrvar%n_fitres_loop(radfit_idx), ctrvar%fitres_range(radfit_idx),   &
             yn_reference_fit,                                                     &
             n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),                &
             fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,           &
             o3fit_cols, o3fit_dcols, target_var(1:ctrvar%n_fincol_idx,ipix),             &
             allfit_cols(1:n_fitvar_rad,ipix), allfit_errs(1:n_fitvar_rad,ipix),   &
             corr_matrix(1:n_fitvar_rad,ipix), yn_bad_pixel, fitspc(1:n_rad_wvl) )

        IF ( yn_bad_pixel ) CYCLE
     END IF

     ! -----------------------------------
     ! Assign pixel values to final arrays
     ! -----------------------------------
     fitconv_flag (ipix,iloop) = INT (radfit_exval, KIND=i2)
     itnum_flag   (ipix,iloop) = INT (radfit_itnum, KIND=i2)
     radfit_chisq (ipix,iloop) = chisquav
     fit_rms      (ipix,iloop) = rms
     column_amount(ipix,iloop) = fitcol
     column_uncert(ipix,iloop) = dfitcol

     ! CCM assign fit residual
     fitspc_out(1:n_rad_wvl,ipix,1) = fitspc(1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,2) = curr_rad_spec(spc_idx,1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,3) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,4) = curr_rad_spec(sig_idx,1:n_rad_wvl)

     IF ( pge_idx == pge_o3_idx ) THEN
        o3_amount(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_cols (o3_t1_idx:o3_t3_idx)
        o3_uncert(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_dcols(o3_t1_idx:o3_t3_idx)
     END IF

  END DO XTrackPix

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE xtrack_radiance_fitting_loop
