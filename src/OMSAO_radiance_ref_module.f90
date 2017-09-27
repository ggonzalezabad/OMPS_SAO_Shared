MODULE OMSAO_radiance_ref_module

  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: i2_missval, r4_missval, &
       r8_missval, downweight, normweight, maxchlen
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, omsao_s_progress, &
       vb_lev_screen, vb_lev_omidebug, error_check


  IMPLICIT NONE

  ! -------------------------------------------------------
  ! Variables connected with  a radiance reference spectrum
  ! -------------------------------------------------------
  LOGICAL                          :: yn_reference_fit
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_reference_lnums

CONTAINS

  SUBROUTINE xtrack_radiance_reference_loop (     &
       yn_radiance_reference, yn_remove_target, nx, nw, fpix, lpix, pge_idx, errstat )

    USE OMSAO_indices_module, ONLY: solar_idx, wvl_idx, spc_idx, sig_idx, &
         o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, shi_idx, squ_idx, &
         ccd_idx, radref_idx
    USE OMSAO_variables_module,  ONLY: database, curr_sol_spec, n_rad_wvl, &
         curr_rad_spec, sol_wav_avg, hw1e, e_asym, n_fitvar_rad, &
         fitvar_rad_saved, n_database_wvl, fitvar_rad, &
         pcfvar, ctrvar
    USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_omidata_module, ONLY: nwav_irrad, irradiance_wght, &
         nwav_rad, n_ins_database_wvl, cross_track_skippix, &
         curr_xtrack_pixnum, n_omi_radwvl, max_rs_idx, ins_database, &
         ins_database_wvl, omi_sol_wav_avg, solcal_pars, radref_wavl, &
         radref_spec, ccdpix_selection, radiance_ccdpix, &
         ccdpix_exclusion, xtrackpix_no, radref_wght, &
         radref_pars, max_calfit_idx, radref_xflag, radref_itnum, &
         radref_chisq, radref_col, radref_rms, radref_dcol, &
         radref_xtrcol

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nw, fpix, lpix
    LOGICAL,           INTENT (IN) :: yn_radiance_reference, yn_remove_target

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, jpix, radfit_exval, radfit_itnum, i
    REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg
    REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
    REAL    (KIND=r8), DIMENSION (n_fitvar_rad)        :: corr_matrix_tmp, allfit_cols_tmp, allfit_errs_tmp
    LOGICAL                  :: yn_skip_pix
    CHARACTER (LEN=maxchlen) :: addmsg
    LOGICAL                                          :: yn_bad_pixel
    INTEGER (KIND=i4), DIMENSION (4)                 :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)                 :: exclud_idx
    INTEGER (KIND=i4)                                :: n_solar_pts
    REAL    (KIND=r8), DIMENSION (1:nw)              :: solar_wgt
    REAL    (KIND=r8), DIMENSION (ctrvar%n_fincol_idx,1:nx) :: target_var 

    ! CCM fitted spectrum now returned from radiance_fit.f90
    REAL    (KIND=r8), DIMENSION (1:nw)              :: fitspctmp

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    locerrstat          = pge_errstat_ok
    target_var          = r8_missval
    fitvar_rad_saved    = ctrvar%fitvar_rad_init


    ! ---------------------------------------------------
    ! Note that this initialization will overwrite valid
    ! results on any second call to this subroutine. This
    ! happens, for example, when YN_RADIANCE_REFERENCE
    ! and YN_REMOVE_TARGET are selected simultaneously.
    ! In that case, however, we write the results to file
    ! before the second call.
    ! ---------------------------------------------------
    radref_pars  (1:max_calfit_idx,1:nx) = r8_missval
    radref_xflag (1:nx)                  = i2_missval
    radref_itnum (1:nx)                  = i2_missval
    radref_chisq (1:nx)                  = r8_missval
    radref_col   (1:nx)                  = r8_missval
    radref_dcol  (1:nx)                  = r8_missval
    radref_rms   (1:nx)                  = r8_missval
    radref_xtrcol(1:nx)                  = r8_missval

    XTrackPix: DO jpix = fpix, lpix

       locerrstat = pge_errstat_ok

       ! -----------------------------------------------------------
       ! The current cross-track pixel number is required further on
       ! when the slit function is computed: The CCD position based
       ! hyper-parameterization requires the knowledge of the row#.
       ! -----------------------------------------------------------
       ipix = jpix
       curr_xtrack_pixnum = ipix

       ! ---------------------------------------------------------------------
       ! If we already determined that this cross track pixel position carries
       ! an error, we don't even have to start processing.
       ! ---------------------------------------------------------------------
       IF ( cross_track_skippix(ipix) ) CYCLE

       n_database_wvl = n_ins_database_wvl(ipix)
       n_omi_radwvl   = nwav_rad      (ipix,0)

       ! ---------------------------------------------------------------------------
       ! For each cross-track position we have to initialize the saved Shift&Squeeze
       ! ---------------------------------------------------------------------------
       saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

       ! ------------------------------------------------------------------
       ! Assign number of irradiance wavelengths and the fitting weights
       ! from the solar wavelength calibration. Why? gga
       ! ------------------------------------------------------------------
       n_solar_pts              = nwav_irrad(ipix)
       solar_wgt(1:n_solar_pts) = irradiance_wght(1:n_solar_pts,ipix)

       ! -----------------------------------------------------
       ! Catch the possibility that N_OMI_RADWVL > N_SOLAR_PTS
       ! -----------------------------------------------------
       IF ( n_omi_radwvl > n_solar_pts ) THEN
          i = n_omi_radwvl - n_solar_pts
          solar_wgt(n_solar_pts+1:n_solar_pts+i) = downweight
          n_solar_pts = n_omi_radwvl
       END IF

       IF ( n_database_wvl > 0 .AND. n_omi_radwvl > 0 ) THEN

          ! ----------------------------------------------
          ! Restore DATABASE from OMI_DATABASE (see above)
          ! ----------------------------------------------
          database (1:max_rs_idx,1:n_database_wvl) = ins_database (1:max_rs_idx,1:n_database_wvl,ipix)

          ! -----------------------------------------------------------------------
          ! Restore solar fitting variables for across-track reference in Earthshine fitting
          ! --------------------------------------------------------------------------------
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
               select_idx(1:4), exclud_idx(1:2),            &
               n_omi_radwvl,                                &
               radref_wavl(1:n_omi_radwvl,ipix),        &
               radref_spec(1:n_omi_radwvl,ipix),        &
               radiance_ccdpix (1:n_omi_radwvl,ipix,0), &
               n_solar_pts, solar_wgt(1:n_solar_pts),       &
               n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_omi_radwvl), rad_spec_avg, &
               yn_skip_pix )

          ! --------------------------------------------------------------------
          ! Update the weights for the Reference/Wavelength Calibration Radiance
          ! --------------------------------------------------------------------
          radref_wght(1:n_omi_radwvl,ipix) = curr_rad_spec(sig_idx,1:n_omi_radwvl)

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
                  pge_idx, ipix, ctrvar%n_fitres_loop(radref_idx), ctrvar%fitres_range(radref_idx), &
                  yn_reference_fit,                                                         &
                  n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),                    &
                  fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,               &
                  o3fit_cols, o3fit_dcols, target_var(1:ctrvar%n_fincol_idx,ipix),                 &
                  allfit_cols_tmp(1:n_fitvar_rad), allfit_errs_tmp(1:n_fitvar_rad),         &
                  corr_matrix_tmp(1:n_fitvar_rad), yn_bad_pixel, fitspctmp(1:n_rad_wvl) )

             yn_reference_fit = .FALSE.

             IF ( yn_bad_pixel ) CYCLE

             WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'RADIANCE Reference #', ipix, &
                  ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; shift = ', &
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

       END IF
      
    END DO XTrackPix

    ! -----------------------------------------
    ! Remove target gas from radiance reference
    ! -----------------------------------------
    IF ( yn_remove_target ) THEN
       ! ----------------------------------------------------------------
       ! Removing the target gas from the radiance reference will alter
       ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
       ! subroutine via MODULE use rather than through the argument list.
       ! ----------------------------------------------------------------
       CALL remove_target_from_radiance (                              &
            fpix, lpix, ctrvar%n_fincol_idx, ctrvar%fincol_idx(1:2,1:ctrvar%n_fincol_idx),  &
            ctrvar%target_npol, target_var(1:ctrvar%n_fincol_idx,fpix:lpix), radref_xtrcol(fpix:lpix) )
       
    END IF

    ! -----------------------------------------------
    ! Update the solar spectrum entry in OMI_DATABASE
    ! -----------------------------------------------
    IF ( yn_radiance_reference ) &
         ins_database (solar_idx,1:n_rad_wvl,fpix:lpix) = radref_spec(1:n_rad_wvl,fpix:lpix)



    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_reference_loop


  SUBROUTINE remove_target_from_radiance (   &
       ipix, jpix, n_fincol_idx, fincol_idx, &
       target_npol, target_var, target_fit   )

    USE OMSAO_omidata_module, ONLY: radref_spec, ins_database, nwav_radref, nwavel_max
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

          !IF ( npol >=0 .AND. nfit > npol ) THEN
          !eps  = -0.1_r8  ! Chose the best-fitting order
           
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

          !CALL dpolft (&
          !     nfit, x(1:fit), y(1:nfit), w(1:nfit), npol, ndeg, eps, yf(1:nfit), ierr, a )

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

          radref_spec(1:nwvl,i) = radref_spec(1:nwvl,i) * &
               EXP(+tmpexp(1:nwvl))

          target_fit(i) = target_fit(i) + yfloc/refspecs_original(k)%NormFactor

       END DO
    END DO

    RETURN
  END SUBROUTINE remove_target_from_radiance

SUBROUTINE create_radiance_reference (nt, nx, nw, locerrstat)

  USE OMSAO_omidata_module, ONLY: &
       ccdpix_selection, nwav_radref, radref_spec, radref_wavl,     &
       radref_qflg, radref_sza, radref_vza, radref_wght,            &
       ccdpix_exclusion, omi_sol_wav_avg 
  USE OMSAO_variables_module, ONLY : ctrvar, pcfvar
  USE OMSAO_omps_reader, ONLY: omps_nmev_type, omps_nmev_reader

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT(IN) :: nt, nx, nw

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4) :: locerrstat

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE(omps_nmev_type) :: omps_data
  INTEGER (KIND=i2)                          :: omps_reader_status_reference
  INTEGER (KIND=i2), DIMENSION(nw)           :: qflg_mask
  REAL    (KIND=r4)                          :: lat_midpt
  REAL    (KIND=r4), DIMENSION(nx,0:nt-1)    :: latr4, tmp_szenith, tmp_vzenith
  REAL    (KIND=r4), DIMENSION(nx)           :: szacount
  REAL    (KIND=r8)                          :: specsum
  REAL    (KIND=r8), DIMENSION(nw)           :: cntr8, radref_wavl_ix
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: local_radref_spec, local_radref_wavl
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: allcount, dumcount
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

  ! ------------------------------------
  ! Read OMPS radiance reference granule
  ! ------------------------------------
  omps_reader_status_reference = omps_nmev_reader(omps_data, &
                                                    pcfvar%l1b_radref_fname)

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

  ! ----------------------------------------------------------------------
  ! Locate the swath line numbers corresponding the center of the latitude
  ! range to average into radiance reference spectrum.
  ! ----------------------------------------------------------------------
  CALL find_swathline_by_latitude ( &
       nx, 0, nt-1, latr4(1:nx,0:nt-1), lat_midpt, &
       xtrange(0:nt-1,1:2), midpt_line, yn_have_scanline )

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

  ! --------------------------------------------------------------------
  ! Now we can average the spectra and the wavelength arrays. Loop over
  ! the block of swath lines in multiples of NLINES_MAX (100 by default)
  ! --------------------------------------------------------------------
  allcount       = 0.0_r8  ; dumcount       = 0.0_r8 ;  szacount = 0.0_r4
  local_radref_wavl    = 0.0_r8  ; local_radref_spec    = 0.0_r8
  radref_sza = 0.0_r4  ; radref_vza = 0.0_r4

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
           
           specsum = 1.0_r8
           
           local_radref_spec(ix,1:nw) = &
                local_radref_spec(ix,1:nw) + tmp_radiance_spec(1:nw,ix,iline)/specsum
           local_radref_wavl(ix,1:nw) = &
                local_radref_wavl(ix,1:nw) + tmp_radiance_wavl(1:nw,ix,iline)
           allcount(ix,1:nw) = allcount(ix,1:nw) + cntr8(1:nw)
           dumcount(ix,1:nw) = dumcount(ix,1:nw) + 1.0_r8
           
           IF ( tmp_szenith(ix,iline) /= r4_missval .AND. &
                tmp_vzenith(ix,iline) /= r4_missval         ) THEN
              radref_sza(ix) = radref_sza(ix) + tmp_szenith(ix,iline)
              radref_vza(ix) = radref_vza(ix) + tmp_vzenith(ix,iline)
              szacount      (ix) = szacount(ix) + 1.0_r4
           END IF
           
        END IF
        
     END DO ! xtrack loop
     
  END DO ! line loop
  
  ! -----------------------------------------------------------
  ! Now for the actual averaging and assignment of final arrays
  ! -----------------------------------------------------------
  DO ix = 1, nx
     ! -----------------------------------
     ! Average the wavelengths and spectra
     ! -----------------------------------
     WHERE ( allcount(ix,1:nw) /= 0.0_r8 )
        local_radref_spec(ix,1:nw) = local_radref_spec(ix,1:nw) / allcount(ix,1:nw)
     END WHERE
     WHERE ( dumcount(ix,1:nw) /= 0.0_r8 )
        local_radref_wavl(ix,1:nw) = local_radref_wavl(ix,1:nw) / dumcount(ix,1:nw)
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
     radref_wavl_ix = local_radref_wavl(ix,1:nw)
     DO j1 = 1, 3, 2
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winwav_lim(j1  ),KIND=r8), 'LE', &
             ccdpix_selection(ix,j1  ) )
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winwav_lim(j1+1),KIND=r8), 'GE', &
             ccdpix_selection(ix,j1+1) )
     END DO

     imin = ccdpix_selection(ix,1)
     imax = ccdpix_selection(ix,4)

     icnt = imax - imin + 1
     nwav_radref(       ix) = icnt
     radref_spec(1:icnt,ix) = local_radref_spec(ix,imin:imax)
     radref_wavl(1:icnt,ix) = radref_wavl_ix(imin:imax)
     radref_qflg(1:icnt,ix) = 0_i2
     radref_wght(1:icnt,ix) = normweight

     ! -----------------------------------------------------------------
     ! Re-assign the average solar wavelength variable, sinfe from here
     ! on we are concerned with radiances.
     ! -----------------------------------------------------------------
     omi_sol_wav_avg(ix) =  SUM( radref_wavl(1:icnt,ix) ) / REAL(icnt, KIND=r8)

     ! ------------------------------------------------------------------
     ! Set weights and quality flags to "bad" for missing spectral points
     ! ------------------------------------------------------------------
     allcount(ix,1:icnt) = allcount(ix,imin:imax)
     WHERE ( allcount(ix,1:icnt) == 0.0_r8 )
        radref_qflg(1:icnt,ix) = 7_i2
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
             nw, radref_wavl_ix, REAL(ctrvar%fit_winexc_lim(1),KIND=r8), 'GE', &
             ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winexc_lim(2),KIND=r8), 'LE', &
             ccdpix_exclusion(ix,2) )
     END IF
     
  END DO

  locerrstat = omps_reader_status_reference  

END SUBROUTINE create_radiance_reference

END MODULE OMSAO_radiance_ref_module
