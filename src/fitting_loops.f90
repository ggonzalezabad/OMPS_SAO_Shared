SUBROUTINE xtrack_radiance_wvl_calibration ( &
     first_pix, last_pix, n_max_rspec, n_comm_wvl_out, errstat )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, &
       max_calfit_idx, max_rs_idx, hwe_idx, asy_idx, sha_idx, &
       shi_idx, squ_idx, ccd_idx, radcal_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, downweight, normweight, &
       i2_missval, i4_missval, r8_missval
  USE OMSAO_variables_module, ONLY: hw1e, e_asym, g_shap, &
       n_rad_wvl, curr_rad_spec, sol_wav_avg, database, fitvar_cal, &
       fitvar_cal_saved, pcfvar, ctrvar
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_data_module, ONLY: &
       cross_track_skippix, nwav_radref, radcal_itnum, &
       radcal_xflag, radcal_chisq, n_ins_database_wvl, &
       solcal_pars, ins_sol_wav_avg,  &
       nwav_irrad, irradiance_wght, irradiance_wavl, &
       irradiance_spec, ins_database, ins_database_wvl, &
       radref_spec, radref_wavl, radref_qflg, radref_wght, &
       nwav_rad, radiance_spec, radiance_wavl, radiance_qflg, &
       ccdpix_selection_rad, ccdpix_exclusion_rad, &
       ccdpix_selection, ccdpix_exclusion, &
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
  INTEGER (KIND=i2) :: local_radcal_itnum
  INTEGER (KIND=i4) :: locerrstat, ipix, radcal_exval, i, n_ref_wvl, cline
  REAL (KIND=r8) :: chisquav, rad_spec_avg
  LOGICAL :: yn_skip_pix, yn_bad_pixel
  CHARACTER (LEN=maxchlen) :: addmsg
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: fitres_out
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: calibration_wavl
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: calibration_spec
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: calibration_wght
  INTEGER (KIND=i2), ALLOCATABLE, DIMENSION(:) :: qflg_mask
  INTEGER (KIND=i4), DIMENSION (4) :: select_idx
  INTEGER (KIND=i4), DIMENSION (2) :: exclud_idx
  REAL (KIND=r8), DIMENSION (n_max_rspec) :: ref_wvl, ref_spc, ref_wgt, rad_wvl

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=31), PARAMETER :: modulename = 'xtrack_radiance_wvl_calibration'

  locerrstat = pge_errstat_ok

  ! Initialize some variables that are to be output
  radcal_itnum = i2_missval ; radcal_xflag = i2_missval
  ins_database = r8_missval ; ins_database_wvl = r8_missval ; n_ins_database_wvl = i4_missval
  
  fitvar_cal_saved(1:max_calfit_idx) = ctrvar%fitvar_rad_init(1:max_calfit_idx)

  ! -------------------------------------------------
  ! Set the number of wavelengths for the common mode
  ! It has to be the maximum between nwav_radref &
  ! nwav_rad
  ! -------------------------------------------------
  IF (ctrvar%yn_radiance_reference) THEN
     n_comm_wvl_out = MAXVAL(nwav_radref)
  ELSE
     n_comm_wvl_out = MAXVAL(nwav_rad)
  END IF

  ! ---------------------------------------------------
  ! Find a line close or nearby the center of the swath
  ! to perform wavelength calibration on that swath if
  ! not using radiance reference.
  ! ---------------------------------------------------
  cline = INT(NrofScanLines / 2, KIND = i4)
  
  ! --------------------------------
  ! Loop over cross-track positions. 
  ! --------------------------------
  XTrackWavCal: DO ipix = first_pix, last_pix

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

     ! -------------------
     ! Initialize database
     ! -------------------
     database = r8_missval

     ! ----------------------------------------------------
     ! Assign number of radiance and irradiance wavelengths
     ! ----------------------------------------------------
     n_irradwvl = nwav_irrad(ipix  )
     IF ( ctrvar%yn_radiance_reference) THEN
        n_radwvl = nwav_radref(ipix)
     ELSE
        n_radwvl = nwav_rad(ipix,cline)
     END IF

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
        n_irradwvl = n_radwvl
        ref_wgt(1:n_radwvl) = normweight
     END IF

     ! -----------------------------------------------------------------
     ! If a Radiance Reference is being used, then it must be calibrated
     ! rather than the swath line that has been read, otherwise we copy
     ! the a line close to the center of the swath into position 0 and
     ! performe the calibration with it.
     ! -----------------------------------------------------------------
     ! First allocate variables
     ! ------------------------
     ALLOCATE(calibration_wavl(1:n_radwvl),calibration_spec(1:n_radwvl), &
          calibration_wght(1:n_radwvl), qflg_mask(1:n_radwvl), &
          fitres_out(1:n_radwvl), &
          stat = locerrstat)
     IF ( ctrvar%yn_radiance_reference ) THEN
        calibration_wavl(1:n_radwvl) = radref_wavl(1:n_radwvl,ipix)
        calibration_spec(1:n_radwvl) = radref_spec(1:n_radwvl,ipix)
        calibration_wght(1:n_radwvl) = radref_wght(1:n_radwvl,ipix)
        select_idx(1:4) = ccdpix_selection(ipix,1:4)
        exclud_idx(1:2) = ccdpix_exclusion(ipix,1:2)
     ELSE
        calibration_wavl(1:n_radwvl) = radiance_wavl(1:n_radwvl,ipix,cline)
        calibration_spec(1:n_radwvl) = radiance_spec(1:n_radwvl,ipix,cline)
        calibration_wght(1:n_radwvl) = ref_wgt(1:n_radwvl)
        qflg_mask(1:n_radwvl) = radiance_qflg(1:n_radwvl,ipix,cline)
        WHERE ( qflg_mask(1:n_radwvl) > 0_i2)
           calibration_wght(1:n_radwvl) = downweight
        END WHERE
        select_idx(1:4) = ccdpix_selection_rad(ipix,cline,1:4)
        exclud_idx(1:2) = ccdpix_exclusion_rad(ipix,cline,1:2)
     END IF

     CALL omi_adjust_radiance_data ( & ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2), n_radwvl, &
          calibration_wavl(1:n_radwvl), calibration_spec(1:n_radwvl), &
          calibration_wght(1:n_radwvl), &
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

     ! ------------------------
     ! Proceed with calibration
     ! ------------------------
     yn_bad_pixel = .FALSE.
     CALL radiance_wavcal ( & ! Radiance wavelength calibration
          ctrvar%n_fitres_loop(radcal_idx), ctrvar%fitres_range(radcal_idx), &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl), &
          radcal_exval, local_radcal_itnum, chisquav, yn_bad_pixel, &
          fitres_out(1:n_rad_wvl), locerrstat )

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
     WRITE (addmsg, '(A,I2,5(A,1PE10.3),2(A,I5))') 'RADIANCE Wavcal    #', ipix, &
          ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; g_sha = ', g_shap, '; shift = ', &
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
     radcal_xflag(ipix) = INT (radcal_exval, KIND=i2)
     radcal_itnum(ipix) = INT (local_radcal_itnum, KIND=i2)
     radcal_chisq(ipix) = chisquav

     ! -------------------------------------------------------------
     ! Save covolved and wavelength calibrated irradiance spectra or
     ! radiance reference spectra to ref_XXX arrays to be used when
     ! creating the database of reference cross sections (below).
     ! -------------------------------------------------------------
     IF ( .NOT. (ctrvar%yn_radiance_reference) ) THEN
        ! Solar irradiance
        n_ref_wvl = n_irradwvl
        ref_wvl(1:n_ref_wvl) = irradiance_wavl(1:n_ref_wvl,ipix)
        ref_spc(1:n_ref_wvl) = irradiance_spec(1:n_ref_wvl,ipix)
        ref_wgt(1:n_ref_wvl) = irradiance_wght(1:n_ref_wvl,ipix)
     ELSE
        ! Radiance reference
        n_ref_wvl = n_rad_wvl
        ref_wvl(1:n_ref_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
        ref_spc(1:n_ref_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl)
        ref_wgt(1:n_ref_wvl) = curr_rad_spec(sig_idx,1:n_rad_wvl)
        ! Update radiance reference
        nwav_radref(ipix)  = n_ref_wvl
        radref_wavl(1:n_ref_wvl,ipix) = curr_rad_spec(wvl_idx,1:n_ref_wvl)
        radref_spec(1:n_ref_wvl,ipix) = curr_rad_spec(spc_idx,1:n_ref_wvl)
        radref_wght(1:n_ref_wvl,ipix) = curr_rad_spec(sig_idx,1:n_ref_wvl)
     END IF

     ! ---------------------------------------------------
     ! Write wavelength calibration results to output file
     ! ---------------------------------------------------
     IF (ctrvar%yn_diagnostic_run) &
          CALL he5_write_radiancewavcal ( n_rad_wvl, ipix, fitvar_cal(squ_idx), &
          fitvar_cal(shi_idx), fitres_out(1:n_rad_wvl), locerrstat)

     ! --------------------
     ! Deallocate variables
     ! --------------------
     DEALLOCATE(calibration_wavl,calibration_spec, calibration_wght, qflg_mask, &
          fitres_out, stat = locerrstat)

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
     ins_database(1:max_rs_idx,1:n_rad_wvl,ipix) = database(1:max_rs_idx,1:n_rad_wvl)
     ins_database_wvl(1:n_rad_wvl,ipix) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     n_ins_database_wvl(ipix) = n_rad_wvl     

  END DO XTrackWavCal

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE xtrack_radiance_wvl_calibration


SUBROUTINE xtrack_radiance_fitting_loop ( yn_common_fit, &
     n_max_rspec, first_pix, last_pix, iloop, &
     n_fitvar_rad, allfit_cols, allfit_errs, corr_matrix, &
     target_var, errstat, fitspc_out, fitspc_out_dim0                 )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, &
       hwe_idx, asy_idx, sha_idx, &
       solar_idx, ccd_idx, radfit_idx
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, downweight
  USE OMSAO_variables_module,  ONLY: database, curr_sol_spec, &
       curr_rad_spec, sol_wav_avg, hw1e, e_asym, g_shap, ctrvar, &
       n_database_wvl
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_data_module, ONLY: n_comm_wvl, &
       column_uncert, column_amount, fit_rms, radfit_chisq, &
       itnum_flag, fitconv_flag, solcal_pars, ins_sol_wav_avg, &
       n_ins_database_wvl, szenith, nwav_rad, &
       cross_track_skippix, radiance_qflg, &
       curr_xtrack_pixnum, radiance_wavl, &
       ccdpix_exclusion_rad, ccdpix_selection_rad, ins_database, &
       ins_database_wvl, max_rs_idx, radiance_spec, &
       radref_wght, irradiance_wght
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok
     
  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: &
       iloop, first_pix, last_pix, n_max_rspec, n_fitvar_rad, &
       fitspc_out_dim0
  LOGICAL :: yn_common_fit

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  REAL    (KIND=r8), INTENT (OUT), DIMENSION (n_fitvar_rad,first_pix:last_pix) :: &
       allfit_cols, allfit_errs, corr_matrix

  ! ---------------------------------------------------------
  ! Optional output variable (fitted variable for target gas)
  ! ---------------------------------------------------------
  REAL (KIND=r8), DIMENSION(ctrvar%n_fincol_idx,first_pix:last_pix), INTENT (OUT) :: target_var

  ! CCM Output fit spectra
  REAL (KIND=r8), DIMENSION(fitspc_out_dim0,first_pix:last_pix,4), INTENT (OUT) :: fitspc_out

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum
  REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg  
  LOGICAL :: yn_skip_pix, yn_cycle_this_pix
  LOGICAL :: yn_bad_pixel
  INTEGER (KIND=i4), DIMENSION (4) :: select_idx
  INTEGER (KIND=i4), DIMENSION (2) :: exclud_idx
  INTEGER (KIND=i4) :: n_solar_pts, n_l1b_wvl, n_l1b_adj_wvl
  REAL    (KIND=r8), DIMENSION (n_max_rspec) :: solar_wvl, ref_wgh
  INTEGER (KIND=i2), DIMENSION (n_max_rspec) :: qflg_mask

  ! CCM Array for holding fitted spectra
  REAL    (KIND=r8), DIMENSION (n_comm_wvl) :: fitspc

  locerrstat = pge_errstat_ok

  XTrackPix: DO ipix = first_pix, last_pix

     curr_xtrack_pixnum = ipix

     ! ---------------------------------------------------------------------
     ! If we already determined that this cross track pixel position carries
     ! an error, we don't even have to start processing.
     ! ---------------------------------------------------------------------
     IF (cross_track_skippix(ipix) .OR. ctrvar%szamax < szenith(ipix,iloop) ) CYCLE
    
     n_database_wvl = n_ins_database_wvl(ipix)
     n_l1b_wvl = nwav_rad(ipix,iloop)

     ! ---------------------------------------------------------------------------
     ! For each cross-track position we have to initialize the saved Shift&Squeeze
     ! ---------------------------------------------------------------------------
     saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

     ! ---------------------------------------------------------------------------
     ! Assign the solar wavelengths. Those should be current in the DATABASE array
     ! and can be taken from there no matter which case - YN_SOLAR_COMP and/or
     ! YN_RADIANCE_REFRENCE we are processing.
     ! ---------------------------------------------------------------------------
     n_solar_pts = n_ins_database_wvl(ipix)
     if (n_solar_pts < 1) cycle

     solar_wvl(1:n_solar_pts) = ins_database_wvl (1:n_solar_pts, ipix)

     CALL check_wavelength_overlap ( &
          n_fitvar_rad, &
          n_solar_pts, solar_wvl (1:n_solar_pts), &
          n_l1b_wvl, radiance_wavl (1:n_l1b_wvl,ipix,iloop), &
          yn_cycle_this_pix )
     IF ( yn_cycle_this_pix .OR. n_database_wvl <= 0 .OR. n_l1b_wvl <= 0) CYCLE

     ! ----------------------------------------------
     ! Restore DATABASE from OMI_DATABASE (see above)
     ! ----------------------------------------------
     database (1:max_rs_idx,1:n_database_wvl) = ins_database (1:max_rs_idx,1:n_database_wvl,ipix)
                 
     ! ---------------------------------------------------------------------------------
     ! Restore solar fitting variables for across-track reference in Earthshine fitting.
     ! ---------------------------------------------------------------------------------
     sol_wav_avg = ins_sol_wav_avg(ipix)
     hw1e = solcal_pars(hwe_idx,ipix)
     e_asym = solcal_pars(asy_idx,ipix)
     g_shap = solcal_pars(sha_idx,ipix)
     curr_sol_spec(wvl_idx,1:n_database_wvl) = ins_database_wvl(1:n_database_wvl,ipix)
     curr_sol_spec(spc_idx,1:n_database_wvl) = ins_database(solar_idx,1:n_database_wvl,ipix)

     ! --------------------------------------
     ! Exclude data from fitting if necessary
     ! --------------------------------------
     select_idx(1:4) = ccdpix_selection_rad(ipix,iloop,1:4)
     exclud_idx(1:2) = ccdpix_exclusion_rad(ipix,iloop,1:2)

     ! --------------------------------------------
     ! Depending on using radiance reference or not
     ! the weighting comes from irradiance_wght or
     ! from radref_wght
     ! --------------------------------------------
     IF ( .NOT. (ctrvar%yn_radiance_reference) ) THEN
        ! Solar irradiance
        ref_wgh(1:n_l1b_wvl) = irradiance_wght(1:n_l1b_wvl,ipix)
     ELSE
        ! Radiance reference
        ref_wgh(1:n_l1b_wvl) = radref_wght(1:n_l1b_wvl,ipix)
     END IF

     ! Bad radiance pixels are downweighted
     qflg_mask(1:n_l1b_wvl) = radiance_qflg(1:n_l1b_wvl,ipix,iloop)
     WHERE ( qflg_mask(1:n_l1b_wvl) > 0_i2)
        ref_wgh(1:n_l1b_wvl) = downweight
     END WHERE

     CALL omi_adjust_radiance_data ( & ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2), &
          n_l1b_wvl, &
          radiance_wavl(1:n_l1b_wvl,ipix,iloop), &
          radiance_spec(1:n_l1b_wvl,ipix,iloop), &
          ref_wgh(1:n_l1b_wvl), &
          n_l1b_adj_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_l1b_wvl),&
          rad_spec_avg, yn_skip_pix )

     ! --------------------
     ! The radiance fitting
     ! --------------------
     fitcol = r8_missval
     dfitcol = r8_missval
     radfit_exval = INT(i2_missval, KIND=i4)
     radfit_itnum = INT(i2_missval, KIND=i4)
     rms = r8_missval

     IF ( MAXVAL(curr_rad_spec(spc_idx,1:n_l1b_adj_wvl)) > 0.0_r8 .AND.     &
          n_l1b_adj_wvl > n_fitvar_rad .AND. (.NOT. yn_skip_pix) ) THEN

        yn_bad_pixel = .FALSE.

        CALL radiance_fit ( &
             ipix, ctrvar%n_fitres_loop(radfit_idx), ctrvar%fitres_range(radfit_idx),   &
             yn_common_fit,                                                     &
             n_l1b_adj_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_l1b_adj_wvl),                &
             fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,           &
             target_var(1:ctrvar%n_fincol_idx,ipix),             &
             allfit_cols(1:n_fitvar_rad,ipix), allfit_errs(1:n_fitvar_rad,ipix),   &
             corr_matrix(1:n_fitvar_rad,ipix), yn_bad_pixel, fitspc(1:n_l1b_adj_wvl) )
        IF ( yn_bad_pixel ) CYCLE

     END IF
!!$     write(*,'(3E13.4,I6)') fitcol, dfitcol, rms, radfit_exval
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
     fitspc_out(1:n_l1b_adj_wvl,ipix,1) = fitspc(1:n_l1b_adj_wvl)
     fitspc_out(1:n_l1b_adj_wvl,ipix,2) = curr_rad_spec(spc_idx,1:n_l1b_adj_wvl)
     fitspc_out(1:n_l1b_adj_wvl,ipix,3) = curr_rad_spec(wvl_idx,1:n_l1b_adj_wvl)
     fitspc_out(1:n_l1b_adj_wvl,ipix,4) = curr_rad_spec(sig_idx,1:n_l1b_adj_wvl)

  END DO XTrackPix

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE xtrack_radiance_fitting_loop
