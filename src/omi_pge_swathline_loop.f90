SUBROUTINE omi_pge_swathline_loop ( &
     nt, nx, n_max_rspec, yn_process,  &
     xtrange, yn_commit, yn_common_fit, & 
     errstat)


  USE OMSAO_precision_module, ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, maxchlen
  USE OMSAO_indices_module, ONLY: n_max_fitpars
  USE OMSAO_variables_module, ONLY: n_fitvar_rad, fitvar_rad_saved, &
       pcfvar, ctrvar
  USE OMSAO_data_module, ONLY:  n_comm_wvl, &
       itnum_flag, fitconv_flag, column_amount, &
       column_uncert, fit_rms
  USE OMSAO_errstat_module, ONLY: omsao_s_progress, pge_errstat_ok, &
       vb_lev_omidebug, vb_lev_screen, omi_smf_setmsg

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT(IN) :: nx, nt, n_max_rspec
  INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange
  LOGICAL, DIMENSION(0:nt-1), INTENT(IN) :: yn_process
  LOGICAL, INTENT(IN) :: yn_commit, yn_common_fit

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: iline, fpix, lpix, ipix, estat, locerrstat, npix
  CHARACTER (LEN=maxchlen) :: addmsg

  ! ---------------------------------------------------------------
  ! Variables to remove target gas from radiance reference spectrum
  ! ---------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (ctrvar%n_fincol_idx,1:nx) :: target_var

  ! ---------------------------------------------------------------------------------
  ! CCM Array to hold (1) Fitted Spec (2) Observed Spec (3) Spec Pos (4) Weight flags
  ! ---------------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (1:n_comm_wvl,1:nx,1:4) :: fitspc_tmp

  ! -------------------------------------
  ! Correlations with main output product 
  ! -------------------------------------
  REAL (KIND=r8), DIMENSION (n_fitvar_rad,nx,0:nt-1) :: &
       all_fitted_columns, all_fitted_errors, correlation_columns

  locerrstat = pge_errstat_ok
  
  ! --------------------------------
  ! Initialize fitting output arrays
  ! --------------------------------
  all_fitted_columns =  r8_missval
  all_fitted_errors  =  r8_missval
  correlation_columns = r8_missval
  fitspc_tmp = r8_missval

  ! ------------------------------------------
  ! Initialize output fields with MissingValue
  ! ------------------------------------------
  itnum_flag(1:nx, 0:nt-1) = i2_missval
  fitconv_flag(1:nx, 0:nt-1) = i2_missval
  column_amount(1:nx, 0:nt-1) = r8_missval
  column_uncert(1:nx, 0:nt-1) = r8_missval
  fit_rms(1:nx, 0:nt-1) = r8_missval

  ! -------------------------------------
  ! Re-initialize saved fitting variables
  ! -------------------------------------
  fitvar_rad_saved(1:n_max_fitpars ) = ctrvar%fitvar_rad_init(1:n_max_fitpars)
  
  ! ---------------
  ! Loop over lines
  ! ---------------
  ScanLines: DO iline = 0, nt-1

     ! -----------------------------------------
     ! Skip if we don't have anything to process
     ! -----------------------------------------
     IF ( .NOT. ( yn_process(iline) ) ) CYCLE

     ! ----------------------------
     ! Not too fancy security check
     ! ----------------------------
     IF (iline > nt-1 ) EXIT ScanLines

     ! ----------------------------------------------------------
     ! Skip this line if it isn't in the list of those to process
     ! ----------------------------------------------------------
     IF ( .NOT. yn_process(iline) ) CYCLE

     ! ------------------
     ! Report on progress
     ! ------------------
     addmsg = ''
     WRITE (addmsg,'(A,I5)') 'Working on scan line', iline
     estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), " ", vb_lev_omidebug )

     fpix = xtrange(iline,1)
     lpix = xtrange(iline,2)
     
     CALL xtrack_radiance_fitting_loop ( yn_common_fit, &
          n_max_rspec, fpix, lpix, iline, &
          n_fitvar_rad, &
          all_fitted_columns (1:n_fitvar_rad,fpix:lpix,iline),   &
          all_fitted_errors  (1:n_fitvar_rad,fpix:lpix,iline),   &
          correlation_columns(1:n_fitvar_rad,fpix:lpix,iline),   &
          target_var(1:ctrvar%n_fincol_idx,fpix:lpix), locerrstat, &
          fitspc_tmp(1:n_comm_wvl,fpix:lpix,1:4), n_comm_wvl )
     
     ipix = (fpix+lpix)/2
     addmsg = ''
     WRITE (addmsg,'(I5, I3,1P,(3E15.5),I5)') iline, ipix, &
          column_amount(ipix, iline), column_uncert(ipix, iline), &
          fit_rms(ipix, iline), MAX(INT(-1,KIND=2),itnum_flag(ipix, iline))
     estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(addmsg), " ", vb_lev_omidebug )
     IF ( pcfvar%verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(addmsg)

     ! ---------------------------------------------------------------
     ! Write out diagnostic results of radiance fitting (if yn_commit)
     ! ---------------------------------------------------------------
     npix = lpix-fpix+1
     IF ( yn_commit ) THEN 
        CALL he5_write_radfit_output (iline, npix, fpix, lpix, n_comm_wvl, &
             all_fitted_columns (1:n_fitvar_rad,fpix:lpix,iline), &
             all_fitted_errors  (1:n_fitvar_rad,fpix:lpix,iline), &
             correlation_columns(1:n_fitvar_rad,fpix:lpix,iline), &
             fitspc_tmp(1:n_comm_wvl,fpix:lpix,1:4), locerrstat )
        errstat = MAX ( errstat, locerrstat )
     END IF

          
  END DO ScanLines

  RETURN
END SUBROUTINE omi_pge_swathline_loop
