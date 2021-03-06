SUBROUTINE omi_adjust_radiance_data (                                   &
     ccdpix_idx, ccdpix_exc, n_radwvl, omi_rad_wvl, omi_rad_spc, &
     curr_sol_weight, n_rad_wvl, curr_rad_spec,   &
     rad_spec_avg, yn_skip_pix )

  USE OMSAO_precision_module, ONLY: i4, r8
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx, ccd_idx
  USE OMSAO_parameters_module, ONLY: downweight, r4_missval
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_solcomp_module, ONLY: solarcomp_pars
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: n_radwvl
  INTEGER (KIND=i4), DIMENSION (2), INTENT (IN) :: ccdpix_exc
  INTEGER (KIND=i4), DIMENSION (4), INTENT (IN) :: ccdpix_idx
  REAL (KIND=r8), DIMENSION (n_radwvl), INTENT (IN) :: omi_rad_wvl, omi_rad_spc
  REAL (KIND=r8), DIMENSION (n_radwvl), INTENT (IN) :: curr_sol_weight

  ! ----------------
  ! Output variables
  ! ----------------
  LOGICAL, INTENT (OUT) :: yn_skip_pix
  INTEGER (KIND=i4), INTENT (OUT) :: n_rad_wvl
  REAL (KIND=r8), INTENT (OUT) :: rad_spec_avg
  REAL (KIND=r8), DIMENSION (ccd_idx,1:n_radwvl), INTENT (OUT) :: curr_rad_spec

  ! ------------------
  ! Modified variables
  ! ------------------
  !INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, locerrstat, imin1, imax1, imin2, imax2, j1, j2
  LOGICAL :: have_good_window
  REAL (KIND=r8), DIMENSION (n_radwvl) :: weightsum

  locerrstat  = pge_errstat_ok
  yn_skip_pix = .FALSE.

  imin1 = ccdpix_idx(1) ; imax1 = ccdpix_idx(4)  ! The total window
  imin2 = ccdpix_idx(2) ; imax2 = ccdpix_idx(3)  ! The fitting window

  ! -------------------------------------------------------------
  ! Assign radiance spectrum to generic variables that are passed
  ! through the fitting routines down to the spectrum function.
  ! -------------------------------------------------------------
  n_rad_wvl                          = n_radwvl
  curr_rad_spec(wvl_idx,1:n_rad_wvl) = omi_rad_wvl(1:n_rad_wvl)
  curr_rad_spec(spc_idx,1:n_rad_wvl) = omi_rad_spc(1:n_rad_wvl)

  ! ---------------------------------------------------------
  ! Compute the weights. This is a bit tedious, as we have to
  ! check for a number of things that can go wrong. We start
  ! out assuming "all is well" and exclude/modify only those
  ! entries that are expected to give us trouble.
  ! ---------------------------------------------------------
  ! EXCEPT, of course, that we also want to exclude any CCD
  ! pixel already excluded from the solar fit. So instead of
  ! starting out with "all is well" we start out with 
  ! "anything that is well in the solar spectrum". 
  !
  ! We assume in this that the number of spectral points is 
  ! equal in both spectra. A reminder to ourselves: Check 
  ! that this indeed the case!
  ! ---------------------------------------------------------
  curr_rad_spec(sig_idx,1:n_rad_wvl) = curr_sol_weight(1:n_rad_wvl)

  ! ---------------------------------------
  ! (1) Anything outside the fitting window
  ! ---------------------------------------
  ! (the CCD indices are absolute positions, i.e., unlikely to be "1:n_sol_wvl")
  ! ----------------------------------------------------------------------------
  IF ( imin2 > imin1 ) curr_rad_spec(sig_idx,1:imin2-imin1+1)         = downweight
  IF ( imax2 < imax1 ) curr_rad_spec(sig_idx,imax2-imin1+1:n_rad_wvl) = downweight
 
  ! ----------------------------------------------------------------------
  ! (2) Any window excluded by the user (specified in fitting control file
  ! ----------------------------------------------------------------------
  IF ( ALL ( ccdpix_exc(1:2) > 0 ) ) THEN
     j1 = ccdpix_exc(1) - imin1 + 1 ; j2 = ccdpix_exc(2) - imin1 + 1
     IF ( j1 >= 1 .AND. j2 <= n_rad_wvl ) curr_rad_spec(sig_idx,j1:j2) = downweight
  END IF
  
  ! ----------------------------------
  ! (3) Wavelengths in ascending order
  ! ----------------------------------
  DO i = 2, n_rad_wvl
     IF ( curr_rad_spec(wvl_idx,i) <= curr_rad_spec(wvl_idx,i-1) ) THEN
        curr_rad_spec(wvl_idx,i) = curr_rad_spec(wvl_idx,i-1) + 0.001_r8
        curr_rad_spec(sig_idx,i) = downweight
     END IF
  END DO

  ! ---------------------------------
  ! (4) No missing values in spectrum
  ! ---------------------------------
  WHERE ( curr_rad_spec(spc_idx,1:n_rad_wvl) <= REAL( r4_missval, KIND=r8 ) )
     curr_rad_spec(sig_idx,1:n_rad_wvl) = downweight
     curr_rad_spec(spc_idx,1:n_rad_wvl) = 0.0_r8
  END WHERE

  ! -------------------------------------------------------------------------------
  ! Translate window limit wavelenghts into indices; making sure that Shift&Squeeze
  ! doesn't shift the wavelength array off the limits we set at the beginning. Do 
  ! this iteratively until we have found good window limits. In the best case, all
  ! is well the first time around, but we might have to adjust the window margins
  ! if we fail to read all the data.
  ! -------------------------------------------------------------------------------
  have_good_window = .TRUE.

  ! --------------------------------------------------
  ! Compute normalization factor for radiance spectrum
  ! --------------------------------------------------
  weightsum = 0.0_r8
  WHERE ( curr_rad_spec(sig_idx,1:n_rad_wvl) /= downweight )
     weightsum = 1.0_r8
  END WHERE

  rad_spec_avg = SUM ( &
       ABS(curr_rad_spec(spc_idx,1:n_rad_wvl))*weightsum(1:n_rad_wvl) ) / &
       MAX(1.0_r8, SUM(weightsum(1:n_rad_wvl)))

  IF ( rad_spec_avg == 0.0_r8 ) rad_spec_avg = 1.0_r8

  ! -------------------------------------------------------------------------
  ! So far we have only taken care of/excluded any negative values in the
  ! spectrum, but there may abnormally high or low positive values also.
  ! Now we check for any values exceeding 100 times the average, which should
  ! be a large enough window to keep anything sensible and reject the real
  ! outliers.
  ! -------------------------------------------------------------------------
  WHERE ( weightsum(1:n_rad_wvl) /= 0.0_r8 .AND. &
       ABS(curr_rad_spec(spc_idx,1:n_rad_wvl)) >= 100.0_r8 * rad_spec_avg )
     weightsum(1:n_rad_wvl) = 0.0_r8
     curr_rad_spec(sig_idx,1:n_rad_wvl) = downweight
     curr_rad_spec(spc_idx,1:n_rad_wvl) = 0.0_r8
  ENDWHERE

  ! ---------------------------------------------------------------------
  ! Recompute the radiance spectrum average, because it may have changed.
  ! ---------------------------------------------------------------------
  rad_spec_avg = SUM ( &
       curr_rad_spec(spc_idx,1:n_rad_wvl)*weightsum(1:n_rad_wvl) ) / &
       MAX(1.0_r8, SUM(weightsum(1:n_rad_wvl)))
  IF ( rad_spec_avg <= 0.0_r8 ) THEN
     yn_skip_pix = .TRUE.
     rad_spec_avg = 1.0_r8
  ELSE
     ! -----------------------------------------
     ! Finally, normalize the radiance spectrum.
     ! -----------------------------------------
     ! There are three possibilities:
     ! (1) If YN_SPECTRUM_NORM = .TRUE. the spectrum will be normalized to 1.
     ! (2) If YN_SPECTRUM_NORM = .FALSE. but YN_SOLAR_COMP = .TRUE. and 
     !     YN_RADIANCE_REFERENCE = .FALSE. then the radiance spectrum will be
     !      normalized with the composite solar norm.
     ! (3) If both of the above are .FALSE. do nothing.
     !
     ! "(2)" assures that radiance and irradiance retain their relative
     ! magnitudes, which means that the Solar Intensity parameter will be more
     ! closely associated with (but not necessarily identical to) the scene albedo.
     ! -----------------------------------------------------------------------------
     IF ( .NOT. ctrvar%yn_spectrum_norm ) THEN                               ! branch for "(2)" or "(3)"
        IF ( ctrvar%yn_solar_comp .AND. (.NOT. ctrvar%yn_radiance_reference) ) THEN ! "(2)"
           rad_spec_avg = solarcomp_pars%SolarNorm
        ELSE                                                          ! "(3)"
           rad_spec_avg = 1.0_r8
        END IF
     END IF
     curr_rad_spec(spc_idx,1:n_rad_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl) / rad_spec_avg
  END IF

  RETURN
END SUBROUTINE omi_adjust_radiance_data


SUBROUTINE omi_create_solcomp_irradiance ( nxt )

  ! ------------------------------------------------------------------
  ! Compute an initial spectrum for an equidistant wavelength array.
  ! This will be used in the solar wavelength calibration to determine
  ! the shift of the composite solar spectrum.
  ! ------------------------------------------------------------------

  USE OMSAO_parameters_module, ONLY: i2, i4, r8
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_solcomp_module, ONLY: soco_compute
  USE OMSAO_data_module, ONLY: nwavel_max, irradiance_spec, irradiance_qflg, &
       irradiance_prec, irradiance_wavl, nwav_irrad, ccdpix_selection, &
       ccdpix_exclusion, ins_sol_wav_avg
  

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nxt

  ! --------------------------------------------
  ! Spacing of the wavelength array to be set up 
  ! --------------------------------------------
  REAL    (KIND=r8), PARAMETER :: dwvl = 0.1_r8

  INTEGER (KIND=i4)                         :: j, ix, nwvl
  REAL    (KIND=r8)                         :: swvl, ewvl
  REAL    (KIND=r8), DIMENSION (nwavel_max) :: tmpwvl


  ! -----------------------------------------------------------
  ! Compute number of wavelengths and assign to temporary array
  ! -----------------------------------------------------------
  swvl = ctrvar%fit_winwav_lim(1) ; ewvl = ctrvar%fit_winwav_lim(4)
  nwvl           = INT ( (ewvl-swvl) / dwvl, KIND=i4 ) + 1
  !tmpwvl(1:nwvl) = swvl + (/ (REAL(j, KIND=r8), j = 0, nwvl) /) * dwvl
  tmpwvl(1:nwvl) = swvl + (/ (REAL(j, KIND=r8), j = 0, nwvl-1) /) * dwvl  ! JED fix

  DO ix = 1, nxt

     nwav_irrad (ix) = nwvl
     irradiance_wavl(1:nwvl,ix) = tmpwvl(1:nwvl)

     ! ---------------------------------------------------------------
     ! Compute the solar spectrum. Note that we are not requesting the
     ! normalized spectrum here, even in cases where we DO want to use
     ! one. Rather, we are keeping the Solar Composite branch as close
     ! as possible to the regular L1b irradiance branch, which at this
     ! point is not normalized. This will be done in a later routine.
     ! ---------------------------------------------------------------
     CALL soco_compute ( &
          .FALSE., ix, nwvl, &
          irradiance_wavl(1:nwvl,ix), irradiance_spec(1:nwvl,ix) )

     irradiance_prec(1:nwvl,ix) = 0.0_r8
     irradiance_qflg(1:nwvl,ix) = 0_i2

     ins_sol_wav_avg(ix) = irradiance_wavl(nwvl/2,ix)


     ! ------------------------------------------------------------------------------
     ! Determine indices included and excluded from the fit. We need to make sure 
     ! that further down the line of the fitting this doesn't screw up things by
     ! introducing incompatible indices into the radiance CCD positions.
     ! ------------------------------------------------------------------------------
     ccdpix_selection(ix,1:4) = -1
     ccdpix_exclusion(ix,1:2) = -1
     DO j = 1, 3, 2
        CALL array_locate_r8 ( &
             nwvl, tmpwvl(1:nwvl), ctrvar%fit_winwav_lim(j  ), 'LE', ccdpix_selection(ix,j  ) )
        CALL array_locate_r8 ( &
             nwvl, tmpwvl(1:nwvl), ctrvar%fit_winwav_lim(j+1), 'GE', ccdpix_selection(ix,j+1) )
     END DO

     IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
             nwvl, tmpwvl(1:nwvl), ctrvar%fit_winexc_lim(1), 'GE', ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nwvl, tmpwvl(1:nwvl), ctrvar%fit_winexc_lim(2), 'LE', ccdpix_exclusion(ix,2) )
     END IF

  END DO


  RETURN
END SUBROUTINE omi_create_solcomp_irradiance


SUBROUTINE compute_fitting_statistics ( &
     ntimes, nxtrack, xtrange, saocol, saodco, saorms, saofcf, saomqf, output, errstat )

  USE OMSAO_precision_module,  ONLY: i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: &
       i2_missval, r8_missval, elsunc_maxiter_eval, elsunc_usrstop_eval, &
       elsunc_less_is_noise, main_qa_good, main_qa_suspect, main_qa_bad
  USE OMSAO_he5_module,       ONLY:  &
       NrOfInputSamples, NrofGoodOutputSamples, NrofSuspectOutputSamples,        &
       NrofBadOutputSamples, NrofConvergedSamples, NrofFailedConvergenceSamples, &
       NrofExceededIterationsSamples, NrofOutofBoundsSamples, NrofMissingSamples, &
       NrofGoodInputSamples, NrofSuspectOutputSamples, NrofBadOutputSamples,      &
       NrofConvergedSamples, NrofFailedConvergenceSamples, &
       PercentGoodOutputSamples, PercentSuspectOutputSamples, &
       PercentBadOutputSamples, QAPercentMissingData, QAPercentOutofBoundsData, &
       AbsolutePercentMissingSamples
  USE OMSAO_errstat_module,   ONLY: vb_lev_screen, pge_errstat_ok
  USE OMSAO_variables_module, ONLY: pcfvar, ctrvar


  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
  INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
  REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco, saorms
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf
  LOGICAL                                          , INTENT (IN) :: output

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saomqf

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Local variables
  ! ----------------
  INTEGER (KIND=i4) :: locerrstat, ix, it, spix, epix, Nrofr8missvalUncertSamples
  REAL    (KIND=r4) :: PercentOutofBoundsSamples
  REAL    (KIND=r8) :: fitcol_avg, rms_avg, dfitcol_avg, nfitcol
  REAL    (KIND=r8) :: col2sig, col3sig


  locerrstat = pge_errstat_ok

  ! ------------------------------------------------------------------
  ! Compute all other fitting statistics variables over two nice loops
  ! ------------------------------------------------------------------
  saomqf                        = i2_missval
  NrofInputSamples              = 0_i4
  NrofGoodInputSamples          = 0_i4
  NrofGoodOutputSamples         = 0_i4
  NrofSuspectOutputSamples      = 0_i4
  NrofBadOutputSamples          = 0_i4
  NrofOutOfBoundsSamples        = 0_i4
  NrofConvergedSamples          = 0_i4
  NrofFailedConvergenceSamples  = 0_i4
  NrofExceededIterationsSamples = 0_i4
  NrofMissingSamples            = 0_i4
  Nrofr8missvalUncertSamples    = 0_i4

  nfitcol    = 0.0_r8
  fitcol_avg = 0.0_r8 ; rms_avg = 0.0_r8 ; dfitcol_avg = 0.0_r8
  DO it = 0, ntimes-1 

     spix = xtrange(it,1) ; epix = xtrange(it,2)
     DO ix = spix, epix

        col2sig = saocol(ix,it)+2.0_r8*saodco(ix,it)
        col3sig = saocol(ix,it)+3.0_r8*saodco(ix,it)

        ! ---------------------------------------------------------
        ! The total number of input samples is simply the number of
        ! pixels in the granule
        ! ---------------------------------------------------------
        NrofInputSamples = NrofInputSamples + 1

        ! ------------------------------------------------------
        ! The Good: Columns are postive within two sigma fitting
        !           uncertainty and the fitting has converged. 
        !           For this "sweet spot" we compute the average
        !           fitting statistics.
        ! ------------------------------------------------------
        IF ( (saofcf(ix,it)      >= INT(elsunc_less_is_noise,KIND=i4)) .AND. &
             (saocol(ix,it)      >  r8_missval                       ) .AND. &
             (saodco(ix,it)      >  r8_missval                       ) .AND. &
             (ABS(saocol(ix,it)) <= ctrvar%max_good_col            ) ) THEN 

           saomqf(ix,it) = main_qa_good

           NrofGoodInputSamples  = NrofGoodInputSamples  + 1
           NrofConvergedSamples  = NrofConvergedSamples  + 1
           NrofGoodOutputSamples = NrofGoodOutputSamples + 1

           fitcol_avg  = fitcol_avg  + saocol(ix,it)
           dfitcol_avg = dfitcol_avg + saodco(ix,it)
           rms_avg     = rms_avg     + saorms(ix,it)
           nfitcol     = nfitcol     + 1.0_r8

           CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Bad: Fitting hasn't converged or columns are negative
        !          within three sigma fitting uncertainty. Note that
        !          pixels can count towards both the number of out-
        !          of bounds and the failed convergence samples.
        ! ----------------------------------------------------------
        IF ( ( saofcf(ix,it) < 0_i2 ) .OR. &
             ( saodco(ix,it) .EQ. r8_missval) ) THEN

           saomqf(ix,it) = main_qa_bad

           NrofBadOutputSamples = NrofBadOutputSamples + 1
           IF ( saofcf(ix,it) < 0 .AND. saofcf(ix,it) >= elsunc_usrstop_eval ) &
                NrofFailedConvergenceSamples  = NrofFailedConvergenceSamples  + 1
           IF ( saofcf(ix,it) == elsunc_maxiter_eval )                      &
                NrofExceededIterationsSamples = NrofExceededIterationsSamples + 1
           IF ( saodco(ix,it) == r8_missval ) &
                Nrofr8missvalUncertSamples = Nrofr8missvalUncertSamples + 1
           CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Ugly: Whatever is left (outside plain missing columns)
        ! ----------------------------------------------------------
        IF ( saocol(ix,it) > r8_missval) THEN

           IF ( (saofcf(ix,it) >= 0_i2 .AND. saofcf(ix,it) < elsunc_less_is_noise) .OR. &
                (ABS(saocol(ix,it)) > ctrvar%max_good_col                        ) .OR. &
                (col3sig < 0) )  THEN

              saomqf(ix,it) = main_qa_suspect
              NrofSuspectOutputSamples = NrofSuspectOutputSamples + 1

              CYCLE
           END IF

        ELSE

           ! ----------------------------------------------------------
           ! The Missing: Not processed because of either missing input
           !              or restrictions on lat, lon, sza, etc.
           ! ----------------------------------------------------------
           NrofMissingSamples = NrofMissingSamples + 1

        END IF

     END DO
  END DO

  ! --------------------------------------------
  ! Now we can compute averages and percentages, 
  ! and write out the final statistics
  ! --------------------------------------------

  IF ( nfitcol >= 1.0_r8 ) THEN
     fitcol_avg  = fitcol_avg  / nfitcol
     rms_avg     = rms_avg     / nfitcol
     dfitcol_avg = dfitcol_avg / nfitcol
  END IF

  PercentGoodOutputSamples      = 100_r4    * &
       REAL(NrofGoodOutputSamples, KIND=r4) / &
       MAX ( 1.0_r4, REAL(NrofInputSamples,  KIND=r4) )

  PercentBadOutputSamples       = 100_r4        * &
       REAL(NrofBadOutputSamples, KIND=r4) / &
       MAX ( 1.0_r4, REAL(NrofInputSamples, KIND=r4) )

  PercentSuspectOutputSamples   =  100.0_r4         * &
       REAL(NrofSuspectOutputSamples, KIND=r4) / &
       MAX ( 1.0_r4, REAL(NrofInputSamples, KIND=r4) )

  PercentOutofBoundsSamples     =  100.0_r4         * &
       REAL(NrofOutofBoundsSamples, KIND=r4) / &
       MAX ( 1.0_r4, REAL(NrofInputSamples, KIND=r4) )

  AbsolutePercentMissingSamples = 100_r4 * &
       REAL(NrofMissingSamples, KIND=r4) / &
       MAX ( 1.0_4, REAL(NrofInputSamples, KIND=r4) )

  QAPercentMissingData     = ANINT ( AbsolutePercentMissingSamples, KIND=i4 )
  QAPercentOutofBoundsData = ANINT ( PercentOutofBoundsSamples,     KIND=i4 )

  ! ------------------------------------------------------------------------
  ! With the above information we can easily determine the Automatic QA Flag
  ! ------------------------------------------------------------------------
  ! CALL set_automatic_quality_flag ( PercentGoodOutputSamples )

  IF ( pcfvar%verb_thresh_lev >= vb_lev_screen ) THEN
     WRITE (*, '(A, 3(1PE15.5))')          'Col-DCol-RMS: ', fitcol_avg, dfitcol_avg, rms_avg
     WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
          MAX(NrofGoodOutputSamples,0), ' of ', MAX(NrofInputSamples,0), ' converged - ', &
          MAX(PercentGoodOutputSamples, 0.0), '%'
     WRITE (*, '(A, I7,A,I7,A)')  'Statistics:   ', &
          MAX(NrofFailedConvergenceSamples,0), ' of ', MAX(NrofInputSamples,0), ' did not converge'
     WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
          MAX(NrofExceededIterationsSamples,0), ' of ', MAX(NrofInputSamples,0), ' exceded maximum number of iterations'
     WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
          MAX(Nrofr8missvalUncertSamples,0), ' of ', MAX(NrofInputSamples,0), ' have r8_missval uncertainty'
     WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
          MAX(NrofOutofBoundsSamples,0), ' of ', MAX(NrofInputSamples,0), ' are ugly'
  END IF

  IF (output) THEN
     CALL he5_write_fitting_statistics ( &
          nxtrack, ntimes, saomqf(1:nxtrack,0:ntimes-1), &
          fitcol_avg, dfitcol_avg, rms_avg, locerrstat)
     errstat = MAX ( locerrstat, errstat )
  ENDIF

  RETURN

END SUBROUTINE compute_fitting_statistics

SUBROUTINE find_endstring ( slen, cstring, istart, iend)
  IMPLICIT NONE

  INTEGER,              INTENT (IN) :: slen, istart
  CHARACTER (LEN=slen), INTENT (IN) :: cstring

  INTEGER, INTENT (OUT) :: iend

  iend = INDEX ( cstring(istart:slen), "," )
  RETURN
END SUBROUTINE find_endstring

SUBROUTINE find_dimension ( slen, cstring, istart, iend, dimstring )
  IMPLICIT NONE

  INTEGER,              INTENT (IN ) :: slen, istart, iend
  CHARACTER (LEN=slen), INTENT (IN) :: cstring

  CHARACTER (LEN=*),    INTENT (OUT) :: dimstring


  dimstring = cstring(istart:iend)

  RETURN
END SUBROUTINE find_dimension

SUBROUTINE slice_string ( slen, instring, istart, ostring )
  IMPLICIT NONE

  INTEGER,              INTENT (IN) :: slen, istart
  CHARACTER (LEN=slen), INTENT (IN) :: instring

  CHARACTER (LEN=*),    INTENT (OUT) :: ostring


  ostring = TRIM(ADJUSTL( instring(istart:slen) ))

  RETURN
END SUBROUTINE slice_string

SUBROUTINE convert_1byte_to_8bits ( nbits, ndim, byte_num, bit_num )

  ! ==========================================================
  ! Takes an NDIM dimensional 1Byte integer BYTE_NUM and
  ! converts it into an NDIM x 8 dimensional integer BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module, ONLY: i1, i4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
  INTEGER (KIND=i1), DIMENSION (ndim), INTENT (IN) :: byte_num

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i1), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i1)                   :: i, powval
  INTEGER (KIND=i1), DIMENSION (ndim) :: tmp_byte

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  bit_num = 0_i1

  ! ------------------------------------------------
  ! Save input variable in TMP_BYTE for modification
  ! ------------------------------------------------
  tmp_byte(1:ndim) = byte_num(1:ndim)

  ! -------------------------------------------------------------------
  ! Starting with the highest power NBYTES-1, subtract powers of 2 and
  ! assign 1 whereever the power fits in the flag number. At the end we
  ! arrive at a 16 BIT binary representation from which we can extract
  ! the surface information.
  ! -------------------------------------------------------------------
  DO i = nbits-1, 0, -1
     powval = 2_i1**i
     IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
           bit_num (1:ndim,i) = 1_i1
           tmp_byte(1:ndim  ) = tmp_byte(1:ndim) - powval
        ENDWHERE
     END IF
  END DO

  RETURN
END SUBROUTINE convert_1byte_to_8bits

SUBROUTINE convert_2bytes_to_16bits ( nbits, ndim, byte_num, bit_num )

  ! ==========================================================
  ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
  ! converts it into an NDIM x 16 dimensional interger BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module, ONLY: i4, i2
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
  INTEGER (KIND=i2), DIMENSION (ndim), INTENT (IN) :: byte_num

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2)                   :: i, powval
  INTEGER (KIND=i2), DIMENSION (ndim) :: tmp_byte

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  bit_num = 0_i2

  ! ------------------------------------------------
  ! Save input variable in TMP_BYTE for modification
  ! ------------------------------------------------
  tmp_byte(1:ndim) = byte_num(1:ndim)

  ! -------------------------------------------------------------------
  ! Starting with the highest power NBYTES-1, subtract powers of 2 and
  ! assign 1 whereever the power fits in the flag number. At the end we
  ! arrive at a 16 BIT binary representation from which we can extract
  ! the surface information.
  ! -------------------------------------------------------------------
  DO i = nbits-1, 0, -1
     powval = 2_i2**i
     IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
           bit_num (1:ndim,i) = 1_i2
           tmp_byte(1:ndim  ) = tmp_byte(1:ndim) - powval
        ENDWHERE
     END IF
  END DO

  RETURN
END SUBROUTINE convert_2bytes_to_16bits

SUBROUTINE radiance_wvl_smoothing ( nxt, nwl, radiance_wavl )

  USE OMSAO_precision_module, ONLY: i4, r4, r8

  ! ---------------
  ! Input Variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nxt, nwl

  ! ------------------
  ! Modified Variables
  ! ------------------
  REAL (KIND=r4), DIMENSION(nwl,nxt), INTENT (INOUT) :: radiance_wavl

  ! ------------------------------
  ! Local Variables and Parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER       :: max_deg = 2
  INTEGER (KIND=i4)                  :: i, ierr, ndeg
  REAL    (KIND=r8)                  :: eps
  REAL    (KIND=r8), DIMENSION (nxt) :: x, y, w, yf
  REAL    (KIND=r8), DIMENSION (3*(nxt+max_deg+1)) :: a

  x(1:nxt) = (/ (REAL(i,KIND=KIND(r8)), i = 1, nxt) /)
  w(1:nxt) = -1.0_r8
  eps      = 0.0_r8

  DO i = 1, nwl
     y(1:nxt) = radiance_wavl(i,1:nxt)
     CALL dpolft (nxt, x, y, w, max_deg, ndeg, eps, yf, ierr, a)
     radiance_wavl(i,1:nxt) = yf(1:nxt)
  END DO

  RETURN
END SUBROUTINE radiance_wvl_smoothing

SUBROUTINE compact_fitting_spectrum ( n_fit_wvl, curr_fit_spec )

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_indices_module,    ONLY: wvl_idx, sig_idx, ccd_idx
  USE OMSAO_parameters_module, ONLY: normweight

  IMPLICIT NONE

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4),                                        INTENT (INOUT) :: n_fit_wvl
  REAL    (KIND=r8), DIMENSION (wvl_idx:ccd_idx,n_fit_wvl), INTENT (INOUT) :: curr_fit_spec


  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)                                        :: i, n_loc_wvl
  REAL    (KIND=r8), DIMENSION (wvl_idx:ccd_idx,n_fit_wvl) :: loc_fit_spec

  n_loc_wvl = 0_i4  ;  loc_fit_spec = 0.0_r8
  DO i = 1, n_fit_wvl
     IF ( curr_fit_spec(sig_idx,i) == normweight ) THEN
        n_loc_wvl = n_loc_wvl + 1
        loc_fit_spec(wvl_idx:ccd_idx,n_loc_wvl) = curr_fit_spec(wvl_idx:ccd_idx,i)
     END IF
  END DO

  ! ---------------------------
  ! Assign the output variables
  ! ---------------------------
  n_fit_wvl     = n_loc_wvl
  curr_fit_spec = loc_fit_spec

  RETURN
END SUBROUTINE compact_fitting_spectrum

SUBROUTINE check_wavelength_overlap ( &
     n_fitvar_rad, n_sol_wvl, irradiance_wvl, n_rad_wvl, radiance_wvl, &
     yn_cycle_this_pix )

  USE OMSAO_precision_module,  ONLY: i4, r8

  IMPLICIT NONE

  ! Input variables
  INTEGER (KIND=i4),                        INTENT (IN) :: n_sol_wvl, n_rad_wvl, n_fitvar_rad
  REAL    (KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: irradiance_wvl
  REAL    (KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: radiance_wvl

  ! Output variable
  LOGICAL, INTENT (OUT) :: yn_cycle_this_pix

  ! Local variables
  INTEGER (KIND=i4) :: j, n_overlap1, n_overlap2

  yn_cycle_this_pix = .FALSE.


  IF ( radiance_wvl(1)         >= irradiance_wvl(n_sol_wvl) .OR. &
       radiance_wvl(n_rad_wvl) <= irradiance_wvl(1)                 ) THEN
     yn_cycle_this_pix = .TRUE.
     RETURN
  END IF

  n_overlap1 = 0
  DO j = 1, n_rad_wvl
     IF ( radiance_wvl(j) >= irradiance_wvl(1)         .AND. &
          radiance_wvl(j) <= irradiance_wvl(n_sol_wvl)         ) n_overlap1 = n_overlap1 + 1
  END DO
  IF ( n_overlap1 < n_fitvar_rad ) THEN
     yn_cycle_this_pix = .TRUE.
     RETURN
  END IF

  n_overlap2 = 0
  DO j = 1, n_sol_wvl
     IF ( irradiance_wvl(j) >= radiance_wvl(1)         .AND. &
          irradiance_wvl(j) <= radiance_wvl(n_rad_wvl)         ) n_overlap2 = n_overlap2 + 1
  END DO
  IF ( n_overlap2 < n_fitvar_rad ) THEN
     yn_cycle_this_pix = .TRUE.
     RETURN
  END IF

  RETURN
END SUBROUTINE check_wavelength_overlap


SUBROUTINE omi_set_xtrpix_range ( &
     nTimes, nXtrack, pixnum_limits, &
     omi_xtrpix_range, first_wc_pix, last_wc_pix, &
     errstat )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_fatal

  IMPLICIT NONE

  ! =======================================================================
  !
  ! Purpose of this routine:
  !
  !    Set first and last cross track pixel position to process
  !
  !
  ! The approach:
  !
  !    For each swath line we define first and last cross track pixel
  !    based 
  ! =======================================================================

  ! ===================================================================
  !
  ! Explanation of subroutine arguments
  !
  !  nTimes .............. number of swath lines
  !  nXtrack ............. number of cross track ("XT") positions
  !  pixnum_limits ....... constraints on XT pixels
  !
  !  first_wc_pix ........ first XT position for wavelengh calibration
  !  last_wc_pix ......... last XT position for wavelengh calibration
  !  omi_xtrpix_range .... (nTimes,2) array of first and last XT for
  !                        radiance fitting
  ! 
  ! ====================================================================

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                        INTENT (IN) :: nTimes, nXtrack
  INTEGER (KIND=i4), DIMENSION(2),          INTENT (IN) :: pixnum_limits

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                          INTENT (OUT) :: first_wc_pix, last_wc_pix
  INTEGER (KIND=i4), DIMENSION(0:nTimes-1,2), INTENT (OUT) :: omi_xtrpix_range

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, locerrstat

  ! ---------------------------
  ! Initialize return variables
  ! ---------------------------
  omi_xtrpix_range(0:nTimes-1,1:2) = -1
  first_wc_pix                     = -1
  last_wc_pix                      = -1

  locerrstat = pge_errstat_ok

  ! ------------------------------------------
  ! Find the range of XT pixels to process
  ! ------------------------------------------
  IF (pixnum_limits(1) < 1 .OR. pixnum_limits(2) < 1) THEN
     first_wc_pix = 1
     last_wc_pix = nXtrack
  ELSE
     first_wc_pix = MAX (                                  1, pixnum_limits(1) )
     last_wc_pix  = MAX ( MIN ( nXtrack,  pixnum_limits(2) ), first_wc_pix     )
  END IF

  ! --------------------------------------------------------------------
  ! We go through the pixels one-by-one. Lots of redundancy, and no
  ! guarantee that we are actually doing the right thing. Proceed with
  ! fingers crossed.
  ! --------------------------------------------------------------------
  DO i = 0, nTimes -1 
     omi_xtrpix_range(i,1) = first_wc_pix
     omi_xtrpix_range(i,2) = last_wc_pix
  END DO

  IF ( ANY(omi_xtrpix_range .LT. 1) .OR. ANY(omi_xtrpix_range .GT. nXtrack)) locerrstat = pge_errstat_fatal
  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE omi_set_xtrpix_range


SUBROUTINE convert_tai_to_utc ( nUTCdim, time_tai, time_utc )

  USE OMSAO_precision_module, ONLY: i2, i4, r8

  IMPLICIT NONE

!  INCLUDE 'PGS_SMF.f'
!  INCLUDE 'PGS_TD_3.f'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER, EXTERNAL :: PGS_TD_TAItoUTC

  ! ---------------
  ! Input Variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nUTCdim
  REAL    (KIND=r8), INTENT (IN) :: time_tai

  ! ---------------
  ! Output Variable
  ! ---------------
  INTEGER (KIND=i2), DIMENSION (nUTCdim), INTENT (OUT) :: time_utc

  ! --------------
  ! Local Variable
  ! --------------
  INTEGER   (KIND=i4 ) :: locerrstat
  CHARACTER (LEN=27)   :: utc_string

  ! -----------------------
  ! Convert TAI to UTC time
  ! -----------------------
  locerrstat = PGS_TD_TAItoUTC ( time_tai, utc_string )

  ! ---------------------------------------------------------------
  ! Now we convert the UTC string to INTEGER values. We had rather
  ! write this to file as a CHAR but at least HE5 v1.6.4 has issues
  ! with that. Hence this kludge.
  ! ---------------------------------------------------------------
  ! The format of the UTC string is: YYYY-MM-DDThh:mm:ss.ddddddZ
  ! where 
  !         YYYY = year             (4 characters)
  !         MM   = month            (2)
  !         DD   = day              (2)
  !         T    = "T"              (separator)
  !         hh   = hour             (2)
  !         mm   = minutes          (2)
  !         ss   = seconds          (2)
  !         d    = decimal fraction (6)
  !         Z    = "Z" (terminator)
  ! ------------------------------------------------------------------
  ! The conversion below is sort-of Low-Brow, but there doesn't seem
  ! to be a simple routine that returns "interger value of a numerical
  ! character". Hence the recourse to old F77 practices.
  ! ------------------------------------------------------------------
  READ ( utc_string( 1: 4), * ) time_utc(1)   ! YYYY
  READ ( utc_string( 6: 7), * ) time_utc(2)   ! MM
  READ ( utc_string( 9:10), * ) time_utc(3)   ! DD
  READ ( utc_string(12:13), * ) time_utc(4)   ! hh
  READ ( utc_string(15:16), * ) time_utc(5)   ! mm
  READ ( utc_string(18:19), * ) time_utc(6)   ! ss

  RETURN
END SUBROUTINE convert_tai_to_utc


SUBROUTINE common_mode_save_vars (                                   &
     frrr, frrr_sav, frrf, frrf_sav, nflr, nflr_sav, nflf, nflf_sav, &
     fl, fl_sav, ll, ll_sav )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! -----------------------------------------------------------------
  ! The argument list is admittedly cryptic, but the purpose of this 
  ! subroutine is decidedly simple: Copy the values of any variable
  ! XXX to XXX_SAVE.
  !
  ! We have to do this for a range of settings when dealing with 
  ! on-line computation of the common mode spectrum.
  ! -----------------------------------------------------------------

  INTEGER (KIND=i4), INTENT (IN)  :: frrr, frrf, nflr, nflf, fl, ll
  INTEGER (KIND=i4), INTENT (OUT) :: frrr_sav, frrf_sav, nflr_sav, nflf_sav, fl_sav, ll_sav

  frrr_sav = frrr
  frrf_sav = frrf
  nflr_sav = nflr
  nflf_sav = nflf
  fl_sav   = fl
  ll_sav   = ll


  RETURN
END SUBROUTINE common_mode_save_vars


SUBROUTINE compute_common_mode ( &
     yn_common_fit, xti, nwvl, fitwvl, fitres, yn_final_call )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_indices_module,   ONLY: max_calfit_idx, comm_idx, mxs_idx
  USE OMSAO_variables_module, ONLY:                                           &
       common_mode_spec,                &
       refspecs_original, ctrvar
  USE OMSAO_data_module,   ONLY:                                           &
       common_spc, common_wvl, common_cnt, n_ins_database_wvl, ins_database,  &
       ccdpix_selection, n_comm_wvl, ins_database_wvl

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  LOGICAL,                             INTENT (IN) :: yn_common_fit, yn_final_call
  INTEGER (KIND=i4),                   INTENT (IN) :: xti, nwvl
  REAL    (KIND=r8), DIMENSION (nwvl), INTENT (IN) :: fitwvl, fitres

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER (KIND=i4) :: i, j, k, errstat
  REAL (KIND=r8) :: comnorm
  LOGICAL :: yn_full_range
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:) :: database_common_mode

  IF ( yn_final_call ) THEN
     ! ---------------------------------------------------
     ! Set the index value of the Common Mode spectrum and
     ! assign values to the fitting parameter arrays
     ! ---------------------------------------------------
     i = max_calfit_idx + (comm_idx-1)*mxs_idx + ctrvar%common_fitpos
     ctrvar%fitvar_rad_init(i) = ctrvar%common_fitvar(1)
     ctrvar%lo_radbnd      (i) = ctrvar%common_fitvar(2)
     ctrvar%up_radbnd      (i) = ctrvar%common_fitvar(3)
     DO i = 1, xti  ! NOTE: "xti == nxtrack" for this call
        j = n_ins_database_wvl(i)
        ! Cycle if we are not using this xtrack
        IF (j < 0) CYCLE
        ALLOCATE(database_common_mode(1:j))
        ! ------------------------------------------
        ! Average the wavelength and spectrum arrays
        ! ------------------------------------------
        k = MAX(1,common_mode_spec%RefSpecCount(i))
        common_mode_spec%RefSpecWavs(i,1:j)  = &
             common_mode_spec%RefSpecWavs(i,1:j) / REAL(k, KIND=r8)
        common_mode_spec%RefSpecData(i,1:j)  = &
             common_mode_spec%RefSpecData(i,1:j) / REAL(k, KIND=r8)

        ! ---------------------------------------
        ! Normalize the Common Mode Spectrum to 1
        ! ---------------------------------------
        ! Skip this for now until we bother with excluding the the
        ! low weights, which otherwise skew the norm.
        ! --------------------------------------------------------
        comnorm = 1.0_r8
        !comnorm = SUM(common_mode_spec%RefSpecData(i,1:j)) / REAL(k, KIND=r8)
        !IF ( comnorm == 0.0_r8 ) comnorm = 1.0_r8
        !common_mode_spec%RefSpecData(i,1:j) = &
        !     common_mode_spec%RefSpecData(i,1:j) / comnorm

        ! -------------------------------------------------
        ! Assign the common mode to the OMI data base array
        ! -------------------------------------------------
        ! First interpolate to database wavelength
        ! ----------------------------------------------------------------------------
        CALL interpolation ( &
             j, common_mode_spec%RefSpecWavs(i,1:j), common_mode_spec%RefSpecData(i,1:j), &
             j, ins_database_wvl(1:j,i), database_common_mode(1:j), &
             'endpoints', 0.0_r8, yn_full_range, errstat )
        ins_database(comm_idx,1:j,i) = database_common_mode(1:j)
        DEALLOCATE(database_common_mode)

        ! --------------------------------------------------------------
        ! Now assign a normalization factor to the original data base of
        ! reference spectra. This is needed in the computation of the
        ! columns and uncertainties of all fitting parameters.
        ! --------------------------------------------------------------
        refspecs_original(comm_idx)%NormFactor = 1.0_r8
     END DO

     ! ------------------------------------
     ! That is all we do for the final call
     ! ------------------------------------
     RETURN
  END IF

  IF ( .NOT. yn_common_fit ) THEN
     ! -------------------------------------------------------------
     ! The Radiance Reference Fit branch saves the wavelength values
     ! and initializes the count and spectrum arrays for the current
     ! cross track position. For goot measure, we also set the 
     ! common mode spectrum fitting parameters to Zero.
     ! -------------------------------------------------------------

     ! ---------------------------------------------------
     ! Set the index value of the Common Mode spectrum and
     ! assign values to the fitting parameter arrays
     ! ---------------------------------------------------
     i = max_calfit_idx + (comm_idx-1)*mxs_idx + ctrvar%common_fitpos
     ctrvar%fitvar_rad_init(i) = 0.0_r8
     ctrvar%lo_radbnd      (i) = 0.0_r8
     ctrvar%up_radbnd      (i) = 0.0_r8

     ! -------------------------------------------------------------
     ! Note that we set COMMON_CNT and COMMON_SPC to 0 for _all_
     ! positions _all_ the time. This is a safety measure just in 
     ! case one of the x-track positions is skipped during the
     ! reference fit due to non-convergence. The Common Mode can
     ! function without wavelengths (actually, we still have to
     ! create a good rationale for them), but not without both
     ! spectrum and cound arrays starting from Zero
     ! -------------------------------------------------------------
     common_wvl(xti, 1:nwvl) = fitwvl(1:nwvl)
     common_cnt              = 0_i4
     common_spc              = 0.0_r8

     common_mode_spec%nPoints      = n_comm_wvl
     common_mode_spec%RefSpecWavs  = 0.0_r8
     common_mode_spec%RefSpecData  = 0.0_r8
     common_mode_spec%RefSpecCount = 0

     common_mode_spec%CCDPixel(xti,1) = INT(ccdpix_selection(xti,1), KIND=i2)
     common_mode_spec%CCDPixel(xti,2) = INT(ccdpix_selection(xti,4), KIND=i2)
  ELSE
     ! --------------------------------------------------------
     ! The Reguar Fitting branch updates the spectrum and count
     ! --------------------------------------------------------
     comnorm = 1.0_r8
     IF ( nwvl > 0 ) THEN
        comnorm = SUM(ABS(fitres(1:nwvl)))/REAL(nwvl, KIND=r8)
        IF ( comnorm == 0.0_r8 ) comnorm = 1.0_r8
     END IF
     
     common_cnt(xti)        = common_cnt(xti) + 1
     common_spc(xti,1:nwvl) = common_spc(xti,1:nwvl) + fitres(1:nwvl)/comnorm
     
     common_mode_spec%RefSpecWavs(xti,1:nwvl)  = &
          common_mode_spec%RefSpecWavs(xti,1:nwvl) + fitwvl(1:nwvl)
     common_mode_spec%RefSpecData(xti,1:nwvl)  = &
          common_mode_spec%RefSpecData(xti,1:nwvl) + fitres(1:nwvl)/comnorm
     common_mode_spec%RefSpecCount(xti)        = &
          common_mode_spec%RefSpecCount(xti) + 1
  END IF

  RETURN
END SUBROUTINE compute_common_mode


SUBROUTINE find_swathline_range ( &
     nt, nx, l1blats, latrange, yn_in_range, errstat )

  USE OMSAO_precision_module, ONLY: i4, r4
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4),                          INTENT (IN) :: nt, nx
  REAL      (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: l1blats
  REAL      (KIND=r4), DIMENSION (2),           INTENT (IN) :: latrange

  ! -----------------------------
  ! Output and Modified variables
  ! -----------------------------
  LOGICAL, DIMENSION (0:nt-1), INTENT (INOUT) :: yn_in_range
  INTEGER (KIND=i4),           INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)                          :: estat, fpix, lpix, midnum, locerrstat
  REAL    (KIND=r4)                          :: midlat
  INTEGER (KIND=i4), DIMENSION (0:nt-1, 1:2) :: xtrange

  locerrstat = pge_errstat_ok
  estat      = pge_errstat_ok

  ! ----------------------------------------------------------------
  ! Read preparatory arrays for determining the range of swath lines
  ! that fall within the desired latitude interval.
  ! ----------------------------------------------------------------
  CALL omi_set_xtrpix_range ( &
       nt, nx, ctrvar%pixnum_lim(3:4), &
       xtrange(0:nt-1,1:2), fpix, lpix, estat    )

  ! ----------------------------------------------------------------------
  ! Determine the range of swath line numbers that go into the radiance
  ! reference spectrum. This is based either on a finite latitude interval
  ! or on a single latitude.
  ! ----------------------------------------------------------------------
  IF ( latrange(1) /= latrange(2) ) THEN
     midlat = SUM(latrange(1:2)) / 2.0_r4
  ELSE
     midlat = latrange(1)
  END IF

  CALL find_swathrange_by_latitude (                      &
       nt, nx, latrange(1), latrange(2),                  &
       l1blats(1:nx,0:nt-1), xtrange(0:nt-1,1:2), midlat, &
       midnum, yn_in_range(0:nt-1)                        )

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE find_swathline_range

SUBROUTINE find_swathline_by_latitude ( &
     nxrr, sline, eline, latr4, lat, xtrange, lnum, yn_found )

  USE OMSAO_precision_module, ONLY: i4, r4
  USE OMSAO_parameters_module, ONLY: r4_missval
  USE OMSAO_errstat_module, ONLY: pge_errstat_error

  IMPLICIT NONE

  ! --------------------------------------------------------------------------
  ! This subroutine returns a swath line number from an OMI swath
  ! based on a given latitude. Subroutine arguments:
  !
  ! nxrr ............. Number of cross-track entries in OMI swath
  ! sline ............ Lowest swath line number
  ! eline ............ Highest swath line number
  ! latr4 ............ Latitude array for whole swath
  ! lat .............. Latitude to locate
  ! xtrange .......... Number of valid cross-track positions in swath
  !                    (possibly smaller than nx)
  ! yn_found ......... TRUE if line number has been found, FALSE otherwise
  ! --------------------------------------------------------------------------

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                                 INTENT (IN) :: nxrr, sline, eline
  REAL    (KIND=r4),                                 INTENT (IN) :: lat
  INTEGER (KIND=i4), DIMENSION (sline:eline,2),      INTENT (IN) :: xtrange
  REAL    (KIND=r4), DIMENSION (1:nxrr,sline:eline), INTENT (IN) :: latr4

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (OUT) :: lnum
  LOGICAL,           INTENT (OUT) :: yn_found

  ! ---------------
  ! Local variables
  ! ---------------
  REAL    (KIND=r4)                   :: diff, mindiff
  REAL    (KIND=r4), DIMENSION (nxrr) :: cntr4, latdiff, loclat
  INTEGER (KIND=i4)                   :: &
       j1, j2, icnt, iline, fpix, lpix, locerr

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  lnum = -1
  yn_found = .FALSE.

  ! --------------------------------------------------------------------------
  ! First, start a bisection of the [0, NLINES-1] interval to find the closest 
  ! match in latitude to the mipoint of the latitude regime to average.
  ! --------------------------------------------------------------------------
  j1 = sline ; j2 = eline-1  ;  icnt = 0
  FindLine: DO WHILE ( .NOT. yn_found )
     icnt  = icnt + 1
     iline = (j1 + j2) / 2

     IF ( iline < sline .OR. iline > eline ) THEN
        locerr = pge_errstat_error
        EXIT FindLine
     END IF

     ! -----------------------------------------------------------------------
     ! Get first and last pixel.
     ! -----------------------------------------------------------------------
     fpix = xtrange(iline,1)
     lpix = xtrange(iline,2)

     loclat(fpix:lpix) = latr4(fpix:lpix,iline)
     cntr4 (fpix:lpix) = 1.0_r4

     WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
        loclat(fpix:lpix) = r4_missval
        cntr4 (fpix:lpix) = 0.0_r4
     END WHERE

     IF ( MAXVAL(loclat(fpix:lpix)) > r4_missval ) THEN
        latdiff(fpix:lpix) = ( loclat(fpix:lpix) - lat ) * cntr4(fpix:lpix)
        diff = SUM(latdiff(fpix:lpix))/SUM(cntr4(fpix:lpix))
        IF ( diff < 0.0_r4 ) THEN
           j1 = iline
        ELSE
           j2 = iline
        END IF
        IF ( ABS(j1 - j2) <= 2 ) THEN
           lnum = (j1 + j2) / 2
           yn_found     = .TRUE.
           EXIT FindLine
        END IF
     ELSE
        ! -------------------------------------------------------------
        ! Reaching here spells trouble. The only thing we can think of
        ! to do here is to alternately increase and decrease the top 
        ! boundary in the hope that we hit upon something. ICNT is 
        ! increased by one for each iteration, so "-1**ICNT" has 
        ! alternating sign between two iterations. The line below first
        ! subtracts 1, then adds 2, then subtracts 3, asf.
        ! -------------------------------------------------------------
        j2 = j2 + icnt * (-1)**icnt
     END IF
  END DO FindLine

  ! ----------------------------------------------------------------
  ! Now we fine-tune the retrieved scan line number by checking +/-2
  ! scan lines on either side.
  ! ----------------------------------------------------------------
  IF ( yn_found ) THEN

     ! -----------------------------------------------------------------------------
     ! MINDIFF will contain the smallest difference found; set to large value first.
     ! -----------------------------------------------------------------------------
     mindiff = REAL((lpix-fpix+1),KIND=r4)*ABS(r4_missval) ; j1 = lnum
     DO iline = lnum-2, lnum+2

        IF ( iline < 0 .OR. iline > eline ) CYCLE

        cntr4(fpix:lpix) = 1.0_r4
        fpix = xtrange(iline,1)
        lpix = xtrange(iline,2)
        loclat(fpix:lpix) = latr4(fpix:lpix,iline)

        WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
           loclat(fpix:lpix) = r4_missval
           cntr4 (fpix:lpix) = 0.0_r4
        END WHERE

        IF ( MAXVAL(loclat(fpix:lpix)) > r4_missval ) THEN
           latdiff(fpix:lpix) = ( loclat(fpix:lpix) - lat ) * cntr4(fpix:lpix)
           diff = SUM(latdiff(fpix:lpix)*latdiff(fpix:lpix))/SUM(cntr4(fpix:lpix))
           IF ( diff < mindiff ) THEN
              j1      = iline
              mindiff = diff
           END IF
        END IF

     END DO
     lnum = j1
  END IF

  RETURN
END SUBROUTINE find_swathline_by_latitude


SUBROUTINE find_swathrange_by_latitude ( &
     nt, nx, latlow, latupp, latr4, xtrange, latmid, latnum, yn_in_range )

  USE OMSAO_precision_module, ONLY: i4, r4
  USE OMSAO_parameters_module, ONLY: r4_missval

  IMPLICIT NONE

  ! ---------------------------------------------------------------------------------
  ! This subroutine returns a swath line number from an OMI swath
  ! based on a given latitude. Subroutine arguments:
  !
  ! nt .................. Number of swath lines
  ! nx .................. Number of cross-track entries in swath
  ! latlow .............. Lowest  latitude to include
  ! latupp .............. Highest latitude to include
  ! latr4 ............... Latitude array for whole swath
  ! xtrange ............. Number of valid cross-track positions in swath
  !                       (possibly smaller than nx)
  ! latmid .............. "Midpoint" latitude closest to average of latlow and latupp
  ! latnum .............. Swath line number of latmid
  ! yn_in_range ......... TRUE if swath line latitude falls between latlow and latupp
  ! ---------------------------------------------------------------------------------

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
  REAL    (KIND=r4),                          INTENT (IN) :: latlow, latupp, latmid
  INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
  REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: latr4

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),           INTENT (OUT)   :: latnum
  LOGICAL, DIMENSION (0:nt-1), INTENT (INOUT) :: yn_in_range

  ! ---------------
  ! Local variables
  ! ---------------
  REAL    (KIND=r4), PARAMETER      :: dlat = 10.0_r4
  INTEGER (KIND=i4)                 :: iline
  REAL    (KIND=r4), DIMENSION (nx) :: cntr4, latdiff, loclat
  REAL    (KIND=r4)                 :: diff, mindiff
  INTEGER (KIND=i4)                 :: fpix, lpix
  LOGICAL                           :: yn_single_lat

  ! -------------------------------------
  ! Initialize output and local variables
  ! -------------------------------------
  latnum = -1 ; mindiff = ABS(r4_missval)

  ! ------------------------------------------------------------
  ! Check whether we are working with a finite latitude interval
  ! or with a single latitude
  ! ------------------------------------------------------------
  yn_single_lat = .TRUE.
  IF ( latlow /= latupp ) yn_single_lat = .FALSE.

  ! --------------------------------------------------------------------------
  ! Owing to the discontiguous nature of NRT L1b storage, we can't assume that
  ! we are working with monotonously increasing latitues. The only possibility
  ! is to go throught the latitude array from start to finish and flag any
  ! swath line that falls within the desired range. So, out with the speedy
  ! bisection and in with the brute-force plowing through.
  ! --------------------------------------------------------------------------
  GetRange: DO iline = 0, nt-1

     ! -----------------------------------------------------------------------
     ! Get first and last pixel.
     ! -----------------------------------------------------------------------
     fpix = xtrange(iline,1)
     lpix = xtrange(iline,2)

     ! --------------------------------------
     ! Store the pixel slice in a local array
     ! --------------------------------------
     loclat(fpix:lpix) = latr4(fpix:lpix,iline)
     cntr4 (fpix:lpix) = 1.0_r4

     ! ---------------------------------------------------
     ! Replace and out-of-bounds entries by missing values
     ! ---------------------------------------------------
     WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
        loclat(fpix:lpix) = r4_missval
        cntr4 (fpix:lpix) = 0.0_r4
     END WHERE


     ! -----------------------------------------------------------------
     ! The trick is to handle both an extended latitude interval and a
     ! single latitude value to locate. Since we have 60 cross-track 
     ! postions in OMI with a somewhat slanted swath, we can't simply
     ! check for equality with a single latitude value.
     !
     ! Solution: If we are looking for a single latitude only, reject
     !           any swath lines that fall outside a 10deg window of
     !           the target. For the rest, we compute the RMS difference
     !           of the swath latitudes to the target latitude.
     ! -----------------------------------------------------------------

     IF ( yn_single_lat ) THEN
        ! -----------------------------------------------------------------------
        ! Check 1: Single latitude. Skip if nothing is within the
        !          [LATMID-DLAT, LAMID+DLAT] interval
        ! -----------------------------------------------------------------------
        IF ( .NOT. (                                        &
             ANY ( loclat(fpix:lpix) >= latmid-dlat ) .AND. &
             ANY ( loclat(fpix:lpix) <= latmid+dlat )      )  ) CYCLE

        ! ----------------------------------------------------------
        ! Update search for the swath line that is closest to LATMID
        ! ----------------------------------------------------------
        latdiff(fpix:lpix) = ( loclat(fpix:lpix) - latmid ) * cntr4(fpix:lpix)
        diff = SQRT( SUM(latdiff(fpix:lpix)*latdiff(fpix:lpix)) ) / SUM(cntr4(fpix:lpix))

        IF ( diff < mindiff ) THEN
           mindiff = diff
           latnum  = iline
        END IF
     ELSE
        ! -----------------------------------------------------------------------
        ! Check 2: Finite latitude interval. Skip if nothing is within the
        !          [LATLOW, LATUPP] interval
        ! -----------------------------------------------------------------------
        IF ( .NOT. (                                    &
             ANY ( loclat(fpix:lpix) >= latlow ) .AND.  &
             ANY ( loclat(fpix:lpix) <= latupp )       )  ) CYCLE

        ! -------------------------
        ! Set yn_in_range to .TRUE.
        ! -------------------------
        yn_in_range(iline) = .TRUE.
     END IF

  END DO GetRange

  ! ------------------------------------------
  ! If working with a single latitude, set the 
  ! corresponding line logical to .TRUE.
  ! ------------------------------------------
  IF ( yn_single_lat .AND. latnum >= 0 .AND. latnum <= nt-1 ) &
       yn_in_range(latnum) = .TRUE.


  RETURN
END SUBROUTINE find_swathrange_by_latitude
