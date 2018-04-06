SUBROUTINE xtrack_radiance_fitting_loop ( yn_common_fit, &
     n_max_rspec, first_pix, last_pix, iloop, &
     n_fitvar_rad, allfit_cols, allfit_errs, corr_matrix, &
     target_var, errstat, fitspc_out, fitspc_out_dim0                 )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, &
       hwe_idx, asy_idx, sha_idx, &
       solar_idx, ccd_idx, radfit_idx
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval!, downweight
  USE OMSAO_variables_module,  ONLY: database, curr_sol_spec, &
       curr_rad_spec, sol_wav_avg, hw1e, e_asym, g_shap, ctrvar, &
       n_database_wvl
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_data_module, ONLY: n_comm_wvl, &
       column_uncert, column_amount, fit_rms, radfit_chisq, &
       itnum_flag, fitconv_flag, solcal_pars, ins_sol_wav_avg, &
       n_ins_database_wvl, szenith, nwav_rad, &
       cross_track_skippix, & !radiance_qflg, &
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
!!$  INTEGER (KIND=i2), DIMENSION (n_max_rspec) :: qflg_mask

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
!     qflg_mask(1:n_l1b_wvl) = radiance_qflg(1:n_l1b_wvl,ipix,iloop)
!     WHERE ( qflg_mask(1:n_l1b_wvl) > 0_i2)
!        ref_wgh(1:n_l1b_wvl) = downweight
!     END WHERE

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
