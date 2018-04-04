MODULE OMSAO_radiance_ref_module

  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: i2_missval, r4_missval, &
       r8_missval, downweight, normweight, maxchlen
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, omsao_s_progress, &
       vb_lev_screen, vb_lev_omidebug, error_check, pge_errstat_fatal, &
       pge_errstat_error


  IMPLICIT NONE

  ! -------------------------------------------------------
  ! Variables connected with  a radiance reference spectrum
  ! -------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_reference_lnums

CONTAINS

  SUBROUTINE xtrack_radiance_reference_loop (     &
       yn_sun, yn_remove_target, nx, nw, fpix, lpix, errstat )

    USE OMSAO_indices_module, ONLY: solar_idx, wvl_idx, spc_idx, sig_idx, &
         hwe_idx, asy_idx, shi_idx, sha_idx, squ_idx, &
         ccd_idx, radref_idx
    USE OMSAO_variables_module,  ONLY: database, curr_sol_spec, n_rad_wvl, &
         curr_rad_spec, hw1e, e_asym, g_shap, n_fitvar_rad, &
         fitvar_rad_saved, n_database_wvl, fitvar_rad, &
         pcfvar, ctrvar
    USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_data_module, ONLY: irradiance_wght, irradiance_wavl, &
         irradiance_spec, nwav_irrad, cross_track_skippix, &
         curr_xtrack_pixnum, n_radwvl, max_rs_idx, ins_database, &
         ins_database_wvl, n_ins_database_wvl, ins_sol_wav_avg, &
         solcal_pars, radref_wavl, radref_spec, ccdpix_selection, &
         nwav_radref, ccdpix_exclusion, radref_wght, radref_pars, &
         max_calfit_idx, radref_xflag, radref_itnum, radref_chisq, &
         radref_col, radref_rms, radref_dcol, radref_xtrcol

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nx, nw, fpix, lpix
    LOGICAL, INTENT (IN) :: yn_sun, yn_remove_target

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum, i
    REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg, rad_wav_avg
    REAL    (KIND=r8), DIMENSION (n_fitvar_rad)        :: corr_matrix_tmp, allfit_cols_tmp, allfit_errs_tmp
    LOGICAL                  :: yn_skip_pix
    CHARACTER (LEN=maxchlen) :: addmsg
    LOGICAL :: yn_bad_pixel
    INTEGER (KIND=i4), DIMENSION (4) :: select_idx
    INTEGER (KIND=i4), DIMENSION (2) :: exclud_idx
    INTEGER (KIND=i4) :: n_solar_pts, npix
    REAL    (KIND=r8), DIMENSION (1:nw) :: solar_wgt
    REAL    (KIND=r8), DIMENSION (ctrvar%n_fincol_idx,1:nx) :: target_var 

    ! CCM fitted spectrum now returned from radiance_fit.f90
    REAL    (KIND=r8), DIMENSION (1:nw)              :: fitspctmp

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    locerrstat = pge_errstat_ok
    target_var = r8_missval
    fitvar_rad_saved = ctrvar%fitvar_rad_init

    ! --------------------------------------------------
    ! Note that this initialization will overwrite valid
    ! results on any second call to this subroutine.
    ! --------------------------------------------------
    radref_pars (1:max_calfit_idx,1:nx) = r8_missval
    radref_xflag (1:nx) = i2_missval
    radref_itnum (1:nx) = i2_missval
    radref_chisq (1:nx) = r8_missval
    radref_col (1:nx) = r8_missval
    radref_dcol (1:nx) = r8_missval
    radref_rms (1:nx) = r8_missval
    npix = lpix-fpix+1 ! For output

    XTrackPix: DO ipix = fpix, lpix

       locerrstat = pge_errstat_ok

       ! -----------------------------------------------------------
       ! The current cross-track pixel number is required further on
       ! when the slit function is computed: The CCD position based
       ! hyper-parameterization requires the knowledge of the row#.
       ! -----------------------------------------------------------
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

       ! ------------------------------------------------------------------
       ! Assign number of irradiance wavelengths and the fitting weights
       ! from the solar wavelength calibration if retrieving against the
       ! Sun spectra, otherwise use the ins_database info
       ! ------------------------------------------------------------------
       n_solar_pts = n_ins_database_wvl(ipix)
       n_database_wvl = n_solar_pts
       IF (yn_sun) THEN
          solar_wgt(1:n_solar_pts) = irradiance_wght(1:n_solar_pts,ipix)
       ELSE
          solar_wgt(1:n_solar_pts) = radref_wght(1:n_solar_pts,ipix)
       END IF

       ! ------------------------
       ! Number of wavelengths in
       ! the radiance reference
       ! ------------------------
       n_radwvl = nwav_radref(ipix)

       ! -------------------------------------------------
       ! Catch the possibility that N_RADWVL > N_SOLAR_PTS
       ! Not a very good fix!!!
       ! -------------------------------------------------
       IF ( n_radwvl > n_solar_pts ) THEN
          i = n_radwvl - n_solar_pts
          solar_wgt(n_solar_pts+1:n_solar_pts+i) = downweight
          n_solar_pts = n_radwvl
       END IF

       IF ( n_solar_pts > 0 .AND. n_radwvl > 0 ) THEN
          ! ----------------------------------
          ! Restore DATABASE from OMI_DATABASE
          ! ----------------------------------
          database (1:max_rs_idx,1:n_database_wvl) = ins_database (1:max_rs_idx,1:n_database_wvl,ipix)

          ! --------------------------------------------------------------------------------
          ! Restore solar fitting variables for across-track reference in Earthshine fitting
          ! We can not retrieve the BrO column in the radiance reference spectrum against it
          ! self. If yn_sun we do it against the solar irradiance, if yn_sun = .FALSE. then
          ! we do it against the radiance reference with the target column substracted in a
          ! first pass and saved in ins_database.
          ! --------------------------------------------------------------------------------
          hw1e   = solcal_pars(hwe_idx,ipix)
          e_asym = solcal_pars(asy_idx,ipix)
          g_shap = solcal_pars(sha_idx,ipix)
          IF (yn_sun) THEN
             curr_sol_spec(wvl_idx,1:n_solar_pts) = irradiance_wavl(1:n_solar_pts,ipix)
             curr_sol_spec(spc_idx,1:n_solar_pts) = irradiance_spec(1:n_solar_pts,ipix)
          ELSE
             curr_sol_spec(wvl_idx,1:n_solar_pts) = ins_database_wvl(1:n_solar_pts,ipix)
             curr_sol_spec(spc_idx,1:n_solar_pts) = ins_database(solar_idx,1:n_solar_pts,ipix)
          END IF
          
          ! --------------------------------------
          ! Prepare radiance reference for fitting
          ! --------------------------------------
          select_idx(1:4) = ccdpix_selection(ipix,1:4)
          exclud_idx(1:2) = ccdpix_exclusion(ipix,1:2)

          CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
               select_idx(1:4), exclud_idx(1:2), &
               n_radwvl, &
               radref_wavl(1:n_radwvl,ipix), &
               radref_spec(1:n_radwvl,ipix), &
               solar_wgt(1:n_radwvl),  &
               n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_radwvl), rad_spec_avg, &
               yn_skip_pix )

          ! --------------------------------------------------------------------
          ! Update the weights for the Reference/Wavelength Calibration Radiance
          ! and the mean fitting window wavelength value.
          ! --------------------------------------------------------------------
          radref_wght(1:n_radwvl,ipix) = curr_rad_spec(sig_idx,1:n_radwvl)
          rad_wav_avg = (SUM ( curr_rad_spec(wvl_idx,1:n_rad_wvl) * curr_rad_spec(sig_idx,1:n_rad_wvl) * &
               curr_rad_spec(sig_idx,1:n_rad_wvl) )  / SUM (curr_rad_spec(sig_idx,1:n_rad_wvl) * &
               curr_rad_spec(sig_idx,1:n_rad_wvl) ) )

          ! --------------------
          ! The radiance fitting
          ! --------------------
          fitcol       = r8_missval
          dfitcol      = r8_missval
          radfit_exval = INT(i2_missval, KIND=i4)
          radfit_itnum = INT(i2_missval, KIND=i4)
          rms          = r8_missval
          addmsg = ''

          IF ( MAXVAL(curr_rad_spec(spc_idx,1:n_rad_wvl)) > 0.0_r8 .AND.     &
               n_rad_wvl > n_fitvar_rad .AND. (.NOT. yn_skip_pix)              ) THEN
             yn_bad_pixel     = .FALSE.
             CALL radiance_fit ( &
                  ipix, ctrvar%n_fitres_loop(radref_idx), ctrvar%fitres_range(radref_idx), &
                  .FALSE., & !Common mode fit logical
                  n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),                    &
                  fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,               &
                  target_var(1:ctrvar%n_fincol_idx,ipix),                 &
                  allfit_cols_tmp(1:n_fitvar_rad), allfit_errs_tmp(1:n_fitvar_rad),         &
                  corr_matrix_tmp(1:n_fitvar_rad), yn_bad_pixel, fitspctmp(1:n_rad_wvl) )
 
             IF ( yn_bad_pixel ) THEN
                 cross_track_skippix(ipix) = .TRUE.
                 CYCLE
              END IF

             WRITE (addmsg, '(A,I2,5(A,1PE10.3),2(A,I5))') 'RADIANCE Reference #', ipix, &
                  ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; g_sha = ', g_shap, '; shift = ', &
                  fitvar_rad(shi_idx), '; squeeze = ', fitvar_rad(squ_idx),&
                  '; exit val = ', radfit_exval, '; iter num = ', radfit_itnum
          ELSE
             WRITE (addmsg, '(A,I2,A)') 'RADIANCE Reference #', ipix, ': Skipped!'
          END IF
          
          ! ------------------
          ! Report on progress
          ! ------------------
          CALL error_check ( &
               0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), vb_lev_omidebug, errstat )
          IF ( pcfvar%verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))
          
          ! -----------------------------------
          ! Assign pixel values to final arrays
          ! -----------------------------------
          radref_pars (1:max_calfit_idx,ipix) = fitvar_rad(1:max_calfit_idx)
          radref_xflag(ipix)                  = INT (radfit_exval, KIND=i2)
          radref_itnum(ipix)                  = INT (radfit_itnum, KIND=i2)
          radref_chisq(ipix)                  = chisquav
          radref_col  (ipix)                  = fitcol
          radref_dcol (ipix)                  = dfitcol
          radref_rms  (ipix)                  = rms

          ! -------------------------------------------------------------------------
          ! Remember weights for the reference radiance, to be used as starting point
          ! in the regular radiance fitting
          ! -------------------------------------------------------------------------
          radref_wght(1:n_rad_wvl,ipix) = curr_rad_spec(sig_idx,1:n_rad_wvl)

          ! ---------------------------------------------
          ! Update wavelength grid with shift and squeeze
          ! ---------------------------------------------
          print*, rad_wav_avg
          radref_wavl(1:n_rad_wvl,ipix) = curr_rad_spec(wvl_idx,1:n_rad_wvl) - fitvar_rad(shi_idx) + &
               (rad_wav_avg * fitvar_rad(squ_idx) / ( 1.0_r8 + fitvar_rad(squ_idx) ) )

       END IF
    END DO XTrackPix

    ! -----------------------------------------
    ! Remove target gas from radiance reference
    ! Only if we are fitting against the Sun
    ! -----------------------------------------
    IF ( yn_remove_target .AND. yn_sun ) THEN
       ! ----------------------------------------------------------------
       ! Removing the target gas from the radiance reference will alter
       ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
       ! subroutine via MODULE use rather than through the argument list.
       ! ----------------------------------------------------------------
       radref_xtrcol(1:nx) = r8_missval
       CALL remove_target_from_radiance (                              &
            fpix, lpix, ctrvar%n_fincol_idx, ctrvar%fincol_idx(1:2,1:ctrvar%n_fincol_idx),  &
            ctrvar%target_npol, target_var(1:ctrvar%n_fincol_idx,fpix:lpix), radref_xtrcol(fpix:lpix) )
    END IF

    ! --------------------------------------------------------
    ! Update values of the solar spectra. We need that for the 
    ! second call to xtrack_prepare_database since we will be 
    ! doing the fitting against the radiance reference. 
    ! This will only have an effect the first time we call 
    ! xtrack_radiance_reference_loop since the database
    ! is not recalculated after the second call
    ! --------------------------------------------------------
    DO ipix = fpix, lpix
       ! ---------------------------------------------------------------------
       ! If we already determined that this cross track pixel position carries
       ! an error, we skip this step
       ! ---------------------------------------------------------------------
       IF ( cross_track_skippix(ipix) ) CYCLE

       ! Calculate average wavelength
       rad_wav_avg = (SUM ( radref_wavl(1:n_rad_wvl,ipix) * radref_wght(1:n_rad_wvl,ipix) * &
            radref_wght(1:n_rad_wvl,ipix) )  / SUM (radref_wght(1:n_rad_wvl,ipix) * &
            radref_wght(1:n_rad_wvl,ipix) ) )
       print*, rad_wav_avg
       nwav_irrad(ipix) = n_rad_wvl
       ins_sol_wav_avg(ipix) = rad_wav_avg
       irradiance_wght(1:nwav_irrad(ipix),ipix) = radref_wght(1:n_rad_wvl,ipix)
       irradiance_wavl(1:nwav_irrad(ipix),ipix) = radref_wavl(1:n_rad_wvl,ipix)
       irradiance_spec(1:nwav_irrad(ipix),ipix) = radref_spec(1:n_rad_wvl,ipix)
    END DO

    ! ----------------------------------------
    ! Write radiance reference fitting results
    ! If we go over it twice results be over
    ! written.
    ! ----------------------------------------
    IF (ctrvar%yn_diagnostic_run) CALL he5_write_radrefcal(npix, fpix, lpix, errstat)
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_reference_loop


  SUBROUTINE remove_target_from_radiance (   &
       ipix, jpix, n_fincol_idx, fincol_idx, &
       target_npol, target_var, target_fit   )

    USE OMSAO_data_module, ONLY: radref_spec, nwav_radref, nwavel_max, ins_database
    USE OMSAO_variables_module, ONLY: refspecs_original
    USE OMSAO_median_module, ONLY: median

    IMPLICIT NONE


    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4),                                     INTENT (IN) :: &
         ipix, jpix, n_fincol_idx, target_npol
    INTEGER (KIND=i4), DIMENSION (2,n_fincol_idx),         INTENT (IN) :: fincol_idx

    ! ------------------
    ! Modified Variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (n_fincol_idx,ipix:jpix), INTENT (INOUT) :: target_var
    REAL (KIND=r8), DIMENSION              (ipix:jpix), INTENT (OUT)   :: target_fit

    ! ------------------------------
    ! Local Variables and Parameters
    ! ------------------------------
    INTEGER (KIND=i4)                         :: i, j, k, l, nwvl
    REAL    (KIND=r8)                         :: yfloc
    REAL    (KIND=r8), DIMENSION (nwavel_max) :: tmpexp

    ! ----------------
    ! DPOLFt variables
    ! ----------------
    INTEGER (KIND=i4)            :: ndeg, ierr, nx, npol, nfit
    REAL    (KIND=r8)            :: eps
    REAL    (KIND=r8), DIMENSION (jpix-ipix+1) :: x, y, yf, w
    REAL    (KIND=r8), DIMENSION (3*((jpix-ipix+1)+target_npol+1)) :: a

    target_fit = 0.0_r8

    npol = target_npol
    nx   = jpix-ipix+1
    DO j = 1, n_fincol_idx
       k = fincol_idx(2,j)
     
       ! ----------------------------------------------------------------------
       ! If we can/have to, fit a cross-track polynomial to the fitted columns,
       ! we do this individually for each FINCOL_IDX and remove the smoothed 
       ! column loading rather than the originally fitted one. In any case, YF
       ! will contain the column values to be removed. Hence the outer loop 
       ! over N_FINCOL_IDX rather than cross-track position.
       ! ----------------------------------------------------------------------
       nfit = 0
       DO i = ipix, jpix
          IF ( target_var(j,i) > r8_missval ) THEN
             nfit    = nfit + 1
             y(nfit) = target_var(j,i)
             w(nfit) = 1.0_r8
          END IF
       END DO

       IF ( nfit /= (jpix-ipix+1) .AND. nfit > 0 ) THEN
          WHERE ( target_var(j,ipix:jpix) <= r8_missval )
             target_var(j,ipix:jpix) = SUM(y(1:nfit))/REAL(nfit,KIND=r8)
          ENDWHERE
       END IF

       ! ----------------------------------------------------
       ! We either fit a polynomial or simply use the Median.
       ! The distinction is made depending on
       !
       ! (a) the order of the cross-track polynomial, and
       ! (b) the number of cross-track points we can fit.
       ! ----------------------------------------------------
       IF ( npol > 0 .AND. jpix-ipix+1 > npol ) THEN
          eps =  0.0_r8  ! Fit the complete NPOL polynomial
          ndeg = npol
          x(1:nx) = (/ ( REAL(i-nx/2, KIND=r8), i = 1, nx ) /) / REAL(nx/2, KIND=r8)
          y(1:nx) = target_var(j,ipix:jpix)
          WHERE ( y(1:nx) > r8_missval )
             w(1:nx) = 1.0_r8
          ELSEWHERE
             w(1:nx) = downweight
          END WHERE
          CALL dpolft (&
               nx, x(1:nx), y(1:nx), w(1:nx), npol, ndeg, eps, yf(1:nx), ierr, a )
       ELSE IF (npol .EQ. 0) THEN
          ! --------------------------------------
          ! Just assing each cross track possition
          ! --------------------------------------
          yf(1:nx) = target_var(j,ipix:jpix)
       ELSE
          ! -----------------------------------------------------------------
          ! The Median is a better choice than the Mean, since the former
          ! is less sensitive to outliers. The Mean may be skewed towards
          ! abnormally high values at the edges of the swath.
          ! -----------------------------------------------------------------
          nfit = nx
          yf(1:nx) = median(nx, target_var(j,ipix:jpix))
       END IF

       DO i = ipix, jpix

          l = i - ipix + 1

          nwvl = nwav_radref(i)

          IF ( yf(l) > r8_missval ) THEN
             yfloc = yf(l)
          ELSE
             yfloc = 0.0_r8
          END IF

          tmpexp(1:nwvl) = yfloc * ins_database(k,1:nwvl,i)
          WHERE ( tmpexp >= MAXEXPONENT(1.0_r8) )
             tmpexp = MAXEXPONENT(1.0_r8) - 1.0_r8
          ENDWHERE
          WHERE ( tmpexp <= MINEXPONENT(1.0_r8) )
             tmpexp = MINEXPONENT(1.0_r8) + 1.0_r8
          ENDWHERE

          ! -------------------------------------
          ! Update the radiance reference with
          ! the target column density substracted
          ! in the reference spectra database.
          ! -------------------------------------
          radref_spec(1:nwvl,i) = radref_spec(1:nwvl,i) * &
               EXP(+tmpexp(1:nwvl))

          target_fit(i) = target_fit(i) + yfloc/refspecs_original(k)%NormFactor

       END DO
    END DO

    RETURN
  END SUBROUTINE remove_target_from_radiance

SUBROUTINE create_radiance_reference (omps_data, nt, nx, nw, locerrstat)

  USE OMSAO_data_module, ONLY: &
       ccdpix_selection, nwav_radref, radref_spec, radref_wavl,     &
       radref_qflg, radref_sza, radref_vza, radref_wght,            &
       ccdpix_exclusion
  USE OMSAO_variables_module, ONLY : ctrvar
  USE OMSAO_omps_reader, ONLY: omps_nmev_type

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  TYPE(omps_nmev_type), INTENT(IN) :: omps_data
  INTEGER (KIND=i4), INTENT(IN) :: nt, nx, nw

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4) :: locerrstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2), DIMENSION(nw)           :: qflg_mask
  REAL    (KIND=r4)                          :: lat_midpt
  REAL    (KIND=r4), DIMENSION(nx,0:nt-1)    :: latr4, tmp_szenith, tmp_vzenith
  REAL    (KIND=r4), DIMENSION(nx)           :: szacount
  REAL    (KIND=r8)                          :: specsum
  REAL    (KIND=r8), DIMENSION(nw)           :: cntr8
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: local_radref_spec, local_radref_wavl
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: allcount
  REAL    (KIND=r8), DIMENSION(nw,nx,0:nt-1) :: tmp_radiance_spec, tmp_radiance_wavl
  INTEGER (KIND=i4), DIMENSION(0:nt-1,2)     :: xtrange
  INTEGER (KIND=i4)                          :: fpix, lpix, midpt_line, iline, ix, &
                                                icnt, imin, imax, j1
  LOGICAL                                    :: yn_have_scanline
  LOGICAL,           DIMENSION (2)           :: yn_have_limits

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = pge_errstat_ok

  latr4(1:nx,0:nt-1)                  = omps_data%Latitude(1:nx,1:nt)
  tmp_szenith(1:nx,0:nt-1)            = omps_data%SolarZenithAngle(1:nx,1:nt)
  tmp_vzenith(1:nx,0:nt-1)            = omps_data%SatelliteZenithAngle(1:nx,1:nt)
  tmp_radiance_spec(1:nw,1:nx,0:nt-1) = omps_data%Radiance(1:nw,1:nx,1:nt)
  tmp_radiance_wavl(1:nw,1:nx,0:nt-1) = omps_data%BandCenterWavelengths(1:nw,1:nx,1:nt)
  
  ! ------------------------------
  ! Initialize some some variables
  ! ------------------------------
  radiance_reference_lnums = -1  ! This will be written to file, hence needs a value
  lat_midpt = SUM ( ctrvar%radref_latrange ) / 2.0_r4

  CALL omi_set_xtrpix_range ( &
       nt, nx, ctrvar%pixnum_lim(3:4), &
       xtrange(0:nt-1,1:2), fpix, lpix, locerrstat )
  IF (locerrstat .GT. pge_errstat_ok) RETURN

  ! ----------------------------------------------------------------------
  ! Locate the swath line numbers corresponding the center of the latitude
  ! range to average into radiance reference spectrum.
  ! ----------------------------------------------------------------------
  CALL find_swathline_by_latitude ( &
       nx, 0, nt-1, latr4(1:nx,0:nt-1), lat_midpt, &
       xtrange(0:nt-1,1:2), midpt_line, yn_have_scanline )
  IF (.NOT. yn_have_scanline) THEN
     WRITE(*,*) 'Radiance reference file does not cover selected latitude range.'
     locerrstat = pge_errstat_fatal
     RETURN
  ENDIF

  ! --------------------------------------------------------------------
  ! If lower and upper bounds of the radiance reference block to average
  ! are identical, then we keep the midpoint line number as the only
  ! reference. Else locate the corresponding swath line numbers.
  ! --------------------------------------------------------------------
  IF ( ctrvar%radref_latrange(1) == ctrvar%radref_latrange(2) ) THEN
     radiance_reference_lnums(1:2) = midpt_line
     yn_have_limits(1:2)           = .TRUE.
  ELSE
     CALL find_swathline_by_latitude ( &
          nx, 0, midpt_line, latr4(1:nx,0:midpt_line), ctrvar%radref_latrange(1), &
          xtrange(0:midpt_line,1:2), radiance_reference_lnums(1), yn_have_limits(1)   )
     CALL find_swathline_by_latitude ( &
          nx, midpt_line, nt-1, latr4(1:nx,midpt_line:nt-1), ctrvar%radref_latrange(2), &
          xtrange(midpt_line:nt-1,1:2), radiance_reference_lnums(2), yn_have_limits(2) )
  END IF

  ! ---------------------------------------------------------
  ! Now we can average the spectra and the wavelength arrays.
  ! ---------------------------------------------------------
  allcount = 0.0_r8; szacount = 0.0_r4
  local_radref_wavl = 0.0_r8; local_radref_spec = 0.0_r8
  radref_sza = 0.0_r4; radref_vza = 0.0_r4

  DO iline = radiance_reference_lnums(1), radiance_reference_lnums(2)

     ! ------------------------------------------------------
     ! Skip this cross-track position if there isn't any data
     ! ------------------------------------------------------
     fpix = xtrange(iline,1)
     lpix = xtrange(iline,2)

     DO ix = 1, nx
        
        IF ( (ix < fpix) .OR. (ix > lpix) ) CYCLE

        ! ----------------------------------------------------------------
        ! For now with the quality flags only check for zero values (good)
        ! ----------------------------------------------------------------
        qflg_mask(1:nw) = omps_data%PixelQualityFlags(1:nw,nx,nt)
        cntr8(1:nw) = 1.0_r8
        WHERE ( qflg_mask(1:nw) > 0_i2 )
           cntr8(1:nw) = 0.0_r8
        END WHERE

        ! ------------------------------------
        ! Only proceed if we have a good value
        ! ------------------------------------
        IF ( ANY ( cntr8(1:nw) > 0.0_r8 ) ) THEN

           tmp_radiance_spec(1:nw,ix,iline) = tmp_radiance_spec(1:nw,ix,iline) * cntr8(1:nw)
           
           specsum = SUM ( tmp_radiance_spec(1:nw,ix,iline) ) / SUM ( cntr8(1:nw) )
           IF ( specsum == 0.0_r8 ) specsum = 1.0_r8
           
           local_radref_spec(ix,1:nw) = &
                local_radref_spec(ix,1:nw) + tmp_radiance_spec(1:nw,ix,iline)/specsum
           local_radref_wavl(ix,1:nw) = &
                local_radref_wavl(ix,1:nw) + tmp_radiance_wavl(1:nw,ix,iline)
           allcount(ix,1:nw) = allcount(ix,1:nw) + cntr8(1:nw)
           
           IF ( tmp_szenith(ix,iline) /= r4_missval .AND. &
                tmp_vzenith(ix,iline) /= r4_missval         ) THEN
              radref_sza(ix) = radref_sza(ix) + tmp_szenith(ix,iline)
              radref_vza(ix) = radref_vza(ix) + tmp_vzenith(ix,iline)
              szacount(ix) = szacount(ix) + 1.0_r4
           END IF
           
        END IF
        
     END DO ! xtrack loop
     
  END DO ! line loop
  
  ! -----------------------------------------------
  ! Do the averaging and assignment of final arrays
  ! -----------------------------------------------
  DO ix = 1, nx

     ! -----------------------------------
     ! Average the wavelengths and spectra
     ! -----------------------------------
     WHERE ( allcount(ix,1:nw) /= 0.0_r8 )
        local_radref_spec(ix,1:nw) = local_radref_spec(ix,1:nw) / allcount(ix,1:nw)
        local_radref_wavl(ix,1:nw) = local_radref_wavl(ix,1:nw) / allcount(ix,1:nw)
     END WHERE

     ! -------------------------------------------
     ! Average the Solar and Viewing Zenith Angles
     ! -------------------------------------------
     IF ( szacount(ix) > 0.0_r4 ) THEN
        radref_sza(ix) = radref_sza(ix) / szacount(ix)
        radref_vza(ix) = radref_vza(ix) / szacount(ix)
     ELSE
        radref_sza(ix) = r4_missval
        radref_vza(ix) = r4_missval             
     END IF

     ! -------------------------------------------------------------------------------
     ! Determine the CCD pixel numbers based on the selected wavelength fitting window
     ! -------------------------------------------------------------------------------     
     DO j1 = 1, 3, 2
        CALL array_locate_r8 ( &
             nw, local_radref_wavl(ix,1:nw), REAL(ctrvar%fit_winwav_lim(j1  ),KIND=r8), 'LE', &
             ccdpix_selection(ix,j1  ) )
        CALL array_locate_r8 ( &
             nw, local_radref_wavl(ix,1:nw), REAL(ctrvar%fit_winwav_lim(j1+1),KIND=r8), 'GE', &
             ccdpix_selection(ix,j1+1) )
     END DO

     imin = ccdpix_selection(ix,1)
     imax = ccdpix_selection(ix,4)
     icnt = imax - imin + 1

     radref_wavl(1:icnt,ix) = local_radref_wavl(ix,imin:imax)
     radref_spec(1:icnt,ix) = local_radref_spec(ix,imin:imax)
     radref_qflg(1:icnt,ix) = 0_i2
     radref_wght(1:icnt,ix) = normweight
     nwav_radref(ix) = icnt

     ! ------------------------------------------------------------------
     ! Set weights and quality flags to "bad" for missing spectral points
     ! ------------------------------------------------------------------
     allcount(ix,1:icnt) = allcount(ix,imin:imax)
     WHERE ( allcount(ix,1:icnt) == 0.0_r8 )
        radref_qflg(1:icnt,ix) = 1_i2
        radref_wght(1:icnt,ix) = downweight
     END WHERE

     ! ------------------------------------------------------------------------------
     ! If any window is excluded, find the corresponding indices. This has to be done
     ! after the array assignements above because we need to know which indices to
     ! exclude from the final arrays, not the complete ones read from the HE4 file.
     ! ------------------------------------------------------------------------------
     ccdpix_exclusion(ix,1:2) = -1
     IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
             nw, local_radref_wavl(ix,1:nw), REAL(ctrvar%fit_winexc_lim(1),KIND=r8), 'GE', &
             ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nw, local_radref_wavl(ix,1:nw), REAL(ctrvar%fit_winexc_lim(2),KIND=r8), 'LE', &
             ccdpix_exclusion(ix,2) )
     END IF
     
  END DO

END SUBROUTINE create_radiance_reference

END MODULE OMSAO_radiance_ref_module
