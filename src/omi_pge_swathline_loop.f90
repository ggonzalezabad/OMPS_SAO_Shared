SUBROUTINE omi_pge_swathline_loop ( &
     pge_idx, nt, nx, n_max_rspec, yn_process,  &
     xtrange, yn_radiance_reference, yn_remove_target, &
     yn_commit, errstat)


  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, maxchlen
  USE OMSAO_indices_module,    ONLY: n_max_fitpars
  USE OMSAO_variables_module,  ONLY:  &
       n_fitvar_rad, fitvar_rad_init, fitvar_rad_saved, &
       n_fincol_idx, pcfvar
  USE OMSAO_omidata_module,    ONLY:  &
       nlines_max, nUTCdim, omi_scanline_no, &
       omi_itnum_flag, omi_fitconv_flag, omi_column_amount, &
       omi_column_uncert, omi_time_utc, omi_fit_rms,  &
       nwavel_max
  USE OMSAO_prefitcol_module
  USE OMSAO_errstat_module
  USE OMSAO_radiance_ref_module, ONLY: remove_target_from_radiance

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nt, n_max_rspec
  INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
  LOGICAL,           DIMENSION (0:nt-1),      INTENT (IN) :: yn_process
  LOGICAL,           INTENT (IN) :: yn_commit, yn_radiance_reference, yn_remove_target

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: iline, fpix, lpix, ipix, estat, locerrstat
  CHARACTER (LEN=maxchlen) :: addmsg

  ! ---------------------------------------------------------------
  ! Variables to remove target gas from radiance reference spectrum
  ! ---------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (n_fincol_idx,1:nx) :: target_var, targsum, targcnt
  REAL (KIND=r8), DIMENSION (1:nx)              :: target_fit, target_col

  ! ---------------------------------------------------------------------------------
  ! CCM Array to hold (1) Fitted Spec (2) Observed Spec (3) Spec Pos (4) Weight flags
  ! ---------------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,4)        :: fitspc_tmp
!!$  REAL (KIND=r8), DIMENSION (n_comm_wvl,nxtrack_max,4,0:nt-1) :: omi_fitspc

  ! -------------------------------------
  ! Correlations with main output product 
  ! -------------------------------------
  REAL (KIND=r8), DIMENSION (n_fitvar_rad,nx,0:nlines_max-1) :: &
       all_fitted_columns, all_fitted_errors, correlation_columns

  locerrstat = pge_errstat_ok

  ! --------------------------------
  ! Initialize fitting output arrays
  ! --------------------------------
  all_fitted_columns =  r8_missval
  all_fitted_errors  =  r8_missval
  correlation_columns = r8_missval


  IF ( yn_radiance_reference .AND. yn_remove_target ) THEN
     target_var = 0.0_r8
     targsum    = 0.0_r8
     targcnt    = 0.0_r8
     target_fit = 0.0_r8
     target_col = 0.0_r8
  END IF

  ! ------------------------------------------
  ! Initialize output fields with MissingValue
  ! ------------------------------------------
  omi_itnum_flag   (1:nx,     0:nt-1) = i2_missval
  omi_fitconv_flag (1:nx,     0:nt-1) = i2_missval
  omi_column_amount(1:nx,     0:nt-1) = r8_missval
  omi_column_uncert(1:nx,     0:nt-1) = r8_missval
  omi_fit_rms      (1:nx,     0:nt-1) = r8_missval
  omi_time_utc     (1:nUTCdim,0:nt-1) = i2_missval

  ! ---------------
  ! Loop over lines
  ! ---------------
  ScanLines: DO iline = 0, nt-1

     ! -----------------------------------------
     ! Skip if we don't have anything to process
     ! -----------------------------------------
     IF ( .NOT. ( yn_process(iline) ) ) CYCLE

     ! -------------------------------------
     ! Re-initialize saved fitting variables
     ! -------------------------------------
     fitvar_rad_saved(1:n_max_fitpars ) = fitvar_rad_init(1:n_max_fitpars)

     ! --------------------------------------------------------------------
     ! Further down, in deeper layers of the algorithm, we require both the
     ! current line in the data block and the absolute swath line number.
     ! Both values are initialized here.
     ! --------------------------------------------------------------------
     omi_scanline_no  = iline

     IF ( omi_scanline_no > nt-1 ) EXIT ScanLines

     ! ----------------------------------------------------------
     ! Skip this line if it isn't in the list of those to process
     ! ----------------------------------------------------------
     IF ( .NOT. yn_process(omi_scanline_no) ) CYCLE

     ! ------------------
     ! Report on progress
     ! ------------------
     addmsg = ''
     WRITE (addmsg,'(A,I5)') 'Working on scan line', omi_scanline_no
     estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), " ", vb_lev_omidebug )

     fpix = xtrange(omi_scanline_no,1)
     lpix = xtrange(omi_scanline_no,2)
     
     CALL xtrack_radiance_fitting_loop (                         &
          n_max_rspec, fpix, lpix, pge_idx, iline,               &
          n_fitvar_rad,                              &
          all_fitted_columns (1:n_fitvar_rad,fpix:lpix,iline),   &
          all_fitted_errors  (1:n_fitvar_rad,fpix:lpix,iline),   &
          correlation_columns(1:n_fitvar_rad,fpix:lpix,iline),   &
          target_var(1:n_fincol_idx,fpix:lpix), locerrstat, fitspc_tmp, nwavel_max )
     ipix = (fpix+lpix)/2
     addmsg = ''
     WRITE (addmsg,'(I5, I3,1P,(3E15.5),I5)') omi_scanline_no, ipix, &
          omi_column_amount(ipix, iline), omi_column_uncert(ipix, iline), &
          omi_fit_rms   (ipix, iline), MAX(INT(-1,KIND=2),omi_itnum_flag(ipix, iline))
     estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(addmsg), " ", vb_lev_omidebug )
     IF ( pcfvar%verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(addmsg)
     
     ! ----------------------------------------------------------------
     ! AMF calculation and update of fitting statistics only need to be
     ! done for the final round throught the common mode iteration loop
     ! ----------------------------------------------------------------
     IF ( yn_commit ) THEN

        CALL he5_write_radfit_output (                       &
             pge_idx, iline, nx, fpix, lpix,         &
             all_fitted_columns (1:n_fitvar_rad,1:nx,iline), &
             all_fitted_errors  (1:n_fitvar_rad,1:nx,iline), &
             correlation_columns(1:n_fitvar_rad,1:nx,iline), &
             fitspc_tmp,locerrstat )
        errstat = MAX ( errstat, locerrstat )

     END IF
     
     
  END DO ScanLines

  RETURN
END SUBROUTINE omi_pge_swathline_loop
