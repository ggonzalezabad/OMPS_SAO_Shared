MODULE OMSAO_solmonthave_module

  ! ====================================================================
  ! This module defines variables associated with the monthly average
  ! solar spectrum and contains necessary subroutines to read and use it
  ! ====================================================================
  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8, C_LONG
  USE OMSAO_omidata_module, ONLY: nxtrack_max
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, &
       pge_errstat_error, file_read_ok, omsao_e_open_solmonave_file, &
       omsao_e_read_solmonave_file, omsao_w_close_solmonave_file, &
       pgs_smf_mask_lev_s, pgsd_io_gen_rseqfrm, f_sep, &
       vb_lev_default, error_check

  IMPLICIT NONE

CONTAINS

  SUBROUTINE omi_read_monthly_average_irradiance (ntimes, nxtrack, errstat)

    ! ---------------------------------------------------------------------
    ! At this point the subroutine is ready to read the solar monthly mean
    ! irradiance from the ASCII files supplied by X. Liu, maybe in the futu
    ! re we will move to an hdf file. However, having to produce a new spec
    ! tra for each day makes that option no too charming
    ! ---------------------------------------------------------------------

    USE OMSAO_variables_module, ONLY: pcfvar, ctrvar
    USE OMSAO_omidata_module, ONLY: &
         nwavel_max, nxtrack_max, omi_irradiance_spec,        &
         omi_irradiance_qflg, omi_irradiance_prec, omi_irradiance_wavl, omi_nwav_irrad, &
         omi_ccdpix_selection, omi_ccdpix_exclusion,             &
         omi_sol_wav_avg, EarthSunDistance
    USE OMSAO_indices_module, ONLY: solmonthave_lun

    IMPLICIT NONE

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    INTEGER (KIND=i4), INTENT (OUT)   :: ntimes, nxtrack

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: funit, ios, dummy, ix, jw, nwavel, nwavelcoef, imin,    &
         imax, icnt, iw, j
    REAL    (KIND=r4), DIMENSION (nwavel_max,nxtrack_max) :: tmp_spc, tmp_wvl, tmp_prc
    INTEGER (KIND=i4), DIMENSION (nwavel_max,nxtrack_max) :: tmp_flg, tmp_n
    ! ------------------------------
    ! Astronomical unit AU in meters
    ! ------------------------------
    REAL (KIND=r8), PARAMETER :: AU_m = 1.495978707d11
    REAL (KIND=r8)            :: sun_earth_distance
    ! --------------------------------------------------------------------
    ! Variable that holds the solar monthly averaged spectra
    ! --------------------------------------------------------------------
    INTEGER   (KIND = i4)     :: nx
    INTEGER   (KIND = i4)     :: nw
    REAL      (KIND = r4), DIMENSION(:,:,:), ALLOCATABLE :: com
    INTEGER   (KIND = i4), DIMENSION(:,:),   ALLOCATABLE :: ncom

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: version, locerrstat

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER (KIND=i4), EXTERNAL :: &
         pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=35), PARAMETER :: modulename = 'omi_read_monthly_average_irradiance' 

    locerrstat = pge_errstat_ok

    ntimes = 1 ; nxtrack = 0 ; nwavel = 0 ; nwavelcoef = 0
    
    ! -------------------------
    ! Open monthly average file
    ! -------------------------
    version = 1
    locerrstat = PGS_IO_GEN_OPENF ( solmonthave_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    CALL error_check ( &
         locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_SOLMONAVE_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(pcfvar%solmonthave_filename)), vb_lev_default, errstat )
    IF (  errstat /= pge_errstat_ok ) RETURN

    ! --------------------------------
    ! Reading the monthly average file
    ! --------------------------------
    ! UV1 channel
    ! -----------
    READ(UNIT=funit, FMT=*, IOSTAT=ios) nx, nw
    IF ( ios /= 0 ) THEN
       CALL error_check ( &
            ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
            modulename//f_sep//TRIM(ADJUSTL(pcfvar%solmonthave_filename)), vb_lev_default, errstat )
       IF (  errstat /= pge_errstat_ok ) RETURN
    END IF

    ALLOCATE (com(nw,3,nx))
    ALLOCATE (ncom(nw,nx))

    DO ix = 1, nx
       READ(UNIT=funit, FMT=*, IOSTAT=ios) dummy
       IF ( ios /= 0 ) THEN
          CALL error_check ( &
               ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
               modulename//f_sep//TRIM(ADJUSTL(pcfvar%solmonthave_filename)), vb_lev_default, errstat )
          IF (  errstat /= pge_errstat_ok ) RETURN
       END IF
       DO jw = 1, nw
          READ(UNIT=funit, FMT=*, IOSTAT=ios) com(jw,1,ix), &
               com(jw,2,ix), com(jw,3,ix), dummy, ncom(jw,ix)
          IF ( ios /= 0 ) THEN
             CALL error_check ( &
                  ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
                  modulename//f_sep//TRIM(ADJUSTL(pcfvar%solmonthave_filename)),      &
                  vb_lev_default, errstat )
             IF (  errstat /= pge_errstat_ok ) RETURN
          END IF
       END DO
     END DO

    ! -----------------------------------------------
    ! Close monthly average file, report SUCCESS read
    ! -----------------------------------------------
    locerrstat = PGS_IO_GEN_CLOSEF ( funit )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    CALL error_check ( &
         locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_SOLMONAVE_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(pcfvar%solmonthave_filename)), vb_lev_default, errstat )

    IF ( errstat >= pge_errstat_error ) RETURN

    ! ---------------------------------------------
    ! Now is time to convert EarthSunDistance to AU
    ! The averaged irradiance is normalized to 1 AU
    ! so below I have to apply the scaling factor:
    ! 1.0 / sun_earth_distance / sun_earth_distance
    ! ---------------------------------------------
    sun_earth_distance = EarthSunDistance / AU_m

    nwavel = nw ; nxtrack = nx
    tmp_wvl(1:nwavel,1:nxtrack) = com(1:nwavel,1,1:nxtrack)
    tmp_spc(1:nwavel,1:nxtrack) = com(1:nwavel,2,1:nxtrack) / sun_earth_distance / sun_earth_distance
    DO iw = 1, nwavel
       DO ix = 1, nxtrack
          tmp_prc(iw,ix) = com(iw,3,ix) / SQRT(REAL(ncom(iw,ix), KIND=r8))
       END DO
    END DO
    tmp_flg(1:nwavel,1:nxtrack) = 128_i2 ! No problem with the pixel set for all of them
    tmp_n(1:nwavel, 1:nxtrack)  = ncom(1:nwavel,1:nxtrack) !Number of irradiance averaged

    DEALLOCATE(com)
    DEALLOCATE(ncom)

    ! -----------------------------------------------------------------------------------------
    ! Writing the readed spectra to the irradiance as if they where comming from an OMI L1 file
    ! -----------------------------------------------------------------------------------------
    ! ----------------------------
    ! Initialize irradiance arrays
    ! ----------------------------
    omi_irradiance_spec (1:nwavel,1:nxtrack) = 0.0_r8 ! Irradiance
    omi_irradiance_prec (1:nwavel,1:nxtrack) = 0.0_r8 ! Irradiance precision
    omi_irradiance_qflg (1:nwavel,1:nxtrack) = 0_i2   ! Irradiance quality flag
    omi_irradiance_wavl (1:nwavel,1:nxtrack) = 0.0_r8 ! Irradiance wavelengths

    ! ----------------------------------------------------
    ! Limit irradiance arrays to fitting window. Check for
    ! strictly ascending wavelengths in the process.
    ! ----------------------------------------------------
    DO ix = 1, nxtrack

       ! -------------------------------------------------------------------------------
       ! Determine the CCD pixel numbers based on the selected wavelength fitting window
       ! -------------------------------------------------------------------------------
       DO j = 1, 3, 2
          CALL array_locate_r4 ( &
               nwavel, tmp_wvl(1:nwavel,ix), REAL(ctrvar%fit_winwav_lim(j  ),KIND=r4), 'LE', &
               omi_ccdpix_selection(ix,j  ) )
          CALL array_locate_r4 ( &
               nwavel, tmp_wvl(1:nwavel,ix), REAL(ctrvar%fit_winwav_lim(j+1),KIND=r4), 'GE', &
               omi_ccdpix_selection(ix,j+1) )
       END DO


       imin = omi_ccdpix_selection(ix,1)
       imax = omi_ccdpix_selection(ix,4)

       icnt = imax - imin + 1
       omi_irradiance_wavl(1:icnt,ix) = REAL(tmp_wvl(imin:imax,ix), KIND=r8)
       omi_irradiance_spec(1:icnt,ix) = REAL(tmp_spc(imin:imax,ix), KIND=r8)
       omi_irradiance_prec(1:icnt,ix) = REAL(tmp_prc(imin:imax,ix), KIND=r8)
       omi_irradiance_qflg(1:icnt,ix) = tmp_flg(imin:imax,ix)
       omi_nwav_irrad (ix) = icnt
       omi_sol_wav_avg(ix) = SUM( tmp_wvl(imin:imax,ix) ) / REAL(icnt, KIND=r8)

       ! ------------------------------------------------------------------------------
       ! If any window is excluded, find the corresponding indices. This has to be done
       ! after the array assignements above because we need to know which indices to
       ! exclude from the final arrays, not the complete ones read from the HE4 file.
       ! ------------------------------------------------------------------------------
       omi_ccdpix_exclusion(ix,1:2) = -1
       IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
          CALL array_locate_r4 ( &
               nwavel, tmp_wvl(1:nwavel,ix), REAL(ctrvar%fit_winexc_lim(1),KIND=r4), 'GE', omi_ccdpix_exclusion(ix,1) )
          CALL array_locate_r4 ( &
               nwavel, tmp_wvl(1:nwavel,ix), REAL(ctrvar%fit_winexc_lim(2),KIND=r4), 'LE', omi_ccdpix_exclusion(ix,2) )
       END IF

    END DO

  END SUBROUTINE omi_read_monthly_average_irradiance

END MODULE OMSAO_solmonthave_module
