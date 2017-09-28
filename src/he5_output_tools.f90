FUNCTION he5_init_swath ( file_name, swath_name, nTimes, nXtrack, nSwLevels ) RESULT ( he5stat )

  !------------------------------------------------------------------------------
  ! This function initializes the HE5 output swath.
  !
  ! Input:
  !   file_name  - Name of HE5 output file
  !   swath_name - Name of swath to be created
  !   nXtrack    - Number of cross-track positions (mind spatial zoom mode!)
  !
  ! Return: he5stat
  !
  ! Variables passed through MODULE:
  !   pge_swath_file_id - id number for HE5 output file (required for closing it)
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
  USE OMSAO_data_module,   ONLY: nclenfit, nUTCdim, nwavel_max
  USE OMSAO_he5_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY: n_fitvar_rad, ctrvar

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'he5_init_swath'
 
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4),      INTENT(IN) :: nTimes, nXtrack, nSwLevels
  CHARACTER (LEN=maxchlen), INTENT(IN) :: file_name, swath_name

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=C_LONG), PARAMETER :: onecl = 1, twocl = 2, fourcl = 4
  INTEGER   (KIND=i4) :: errstat

  he5stat = pge_errstat_ok
  errstat = pge_errstat_ok

  ! ---------------------------------------------------------------
  ! Open HE5 output file and check PGE_SWATH_FILE_ID ( -1 if error)
  ! ---------------------------------------------------------------
  pge_swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_trunc )
  IF ( pge_swath_file_id == he5_stat_fail ) THEN
     CALL error_check ( &
          0, 1, pge_errstat_fatal, OMSAO_F_HE5SWOPEN, modulename, vb_lev_default, he5stat )
     RETURN
  END IF

  ! ------------------------------------------------------
  ! Create HE5 swath and check PGE_SWATH_ID ( -1 if error)
  ! ------------------------------------------------------
  pge_swath_id = HE5_SWcreate ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( pge_swath_id == he5_stat_fail ) THEN
     CALL error_check ( &
          0, 1, pge_errstat_fatal, OMSAO_F_HE5SWCREATE, modulename, vb_lev_default, he5stat )
     RETURN
  END IF
  
  ! ---------------------------------------------------------
  ! Define new dimensions in HE5 swath and check error status
  ! ---------------------------------------------------------
  errstat = HE5_SWdefdim  ( pge_swath_id, ntc,   INT(nTimes,     KIND=C_LONG) )
  errstat = HE5_SWdefdim  ( pge_swath_id, nxc,   INT(nXtrack,    KIND=C_LONG) )
  errstat = HE5_SWdefdim  ( pge_swath_id, nlc,   INT(nSwLevels,  KIND=C_LONG) )
  errstat = HE5_SWdefdim  ( pge_swath_id, ntcp1, INT(nTimes+1,   KIND=C_LONG) )
  errstat = HE5_SWdefdim  ( pge_swath_id, nxcp1, INT(nXtrack+1,  KIND=C_LONG) )
  errstat = HE5_SWDefdim  ( pge_swath_id, nutcd, INT(nUTCdim,    KIND=C_LONG) )
  errstat = HE5_SWDefdim  ( pge_swath_id, "1",   onecl                        )
  errstat = HE5_SWDefdim  ( pge_swath_id, "2",   twocl                        )
  errstat = HE5_SWDefdim  ( pge_swath_id, "4",   fourcl                       )
  errstat = HE5_SWdefdim  ( pge_swath_id, nwalm, INT(nwavel_max, KIND=C_LONG) )

  ! ---------------------------------------------------------
  ! Dimensions for Diagnostic Fields
  ! ---------------------------------------------------------
  IF ( ctrvar%yn_diagnostic_run ) THEN
     errstat = HE5_SWdefdim  ( pge_swath_id, nfv,   INT(n_fitvar_rad,   KIND=C_LONG) )
     errstat = HE5_SWDefdim  ( pge_swath_id, ncv,   INT(nclenfit,       KIND=C_LONG) )
     errstat = HE5_SWdefdim  ( pge_swath_id, nwcp,  INT(max_calfit_idx, KIND=C_LONG) )

     ! CCM for refspec database
     errstat = HE5_SWdefdim  ( pge_swath_id, nrspc, INT(max_rs_idx,     KIND=C_LONG) )
  END IF

  CALL error_check ( &
       errstat, he5_stat_ok, pge_errstat_fatal, OMSAO_F_HE5SWDEFDIM, &
       modulename, vb_lev_default, he5stat )
  IF ( he5stat >= pge_errstat_error ) RETURN


  RETURN
END FUNCTION he5_init_swath

FUNCTION he5_define_fields ( pge_idx, swath_name, nTimes, nXtrack, nSwLevels ) RESULT ( he5stat )

  !------------------------------------------------------------------------------
  ! This function defines fields in the HE5 output swath.
  !
  ! Input:
  !   pge_idx    - Index of current PGE
  !
  ! Return: he5stat
  !
  ! Variables passed through MODULE:
  !   pge_swath_file_id - id number for HE5 output file (required for closing it)
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_indices_module, ONLY: sao_molecule_names
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_data_module, ONLY: n_field_maxdim
  USE OMSAO_he5_module
  USE OMSAO_he5_datafields_module, ONLY: geo_he5fields, sol_calfit_he5fields, &
       rad_calfit_he5fields, rad_reffit_he5fields, comdata_he5fields, &
       diagnostic_he5fields, voc_he5fields, sw_he5fields
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, he5_stat_ok, &
       error_check, vb_lev_default, he5_stat_fail, omsao_e_he5swattach, &
       omsao_e_he5swdeffld, pge_errstat_fatal

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=17), PARAMETER :: modulename = 'he5_define_fields'
 
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4),      INTENT(IN) :: pge_idx, nTimes, nXtrack, nSwLevels
  CHARACTER (LEN=maxchlen), INTENT(IN) :: swath_name

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                             :: i, errstat
  CHARACTER (LEN=4)                               :: molstr

  ! ---------------------------------------------
  ! Variables related to compression and chunking
  ! ---------------------------------------------
  INTEGER   (KIND=i4)                              :: n_chunk_dim
  LOGICAL                                          :: yn_compress_field
  INTEGER   (KIND=C_LONG), DIMENSION (n_field_maxdim) :: chunk_dim

  he5stat = pge_errstat_ok
  errstat = pge_errstat_ok

  ! ----------------------------------------------------
  ! Some fields have the PGE molecule name added to them
  ! ----------------------------------------------------
  molstr = TRIM(ADJUSTL(sao_molecule_names(pge_idx)))

  ! -------------------------------------------------------------
  ! Define geolocation fields in HE5 swath and check error status
  ! -------------------------------------------------------------
  DO i = 1, n_gfields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(geo_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )
     END IF

     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( geo_he5fields(i), errstat )

     geo_he5fields(i)%Swath_ID = pge_swath_id
     errstat = HE5_SWdefgfld (                               &
          geo_he5fields(i)%Swath_ID,                         &
          TRIM(ADJUSTL(geo_he5fields(i)%Name)),              &
          TRIM(ADJUSTL(geo_he5fields(i)%Dimensions)),   " ", &
          geo_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO

  ! ----------------------------------------------------------------------
  ! Create fields for solar and radiance wavelength calibration parameters
  ! ----------------------------------------------------------------------
  DO i = 1, n_solcal_fields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(sol_calfit_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )
     END IF

     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( sol_calfit_he5fields(i), errstat )

     sol_calfit_he5fields(i)%Swath_ID = pge_swath_id
     errstat = HE5_SWdefdfld (                                      &
          sol_calfit_he5fields(i)%Swath_ID,                         &
          TRIM(ADJUSTL(sol_calfit_he5fields(i)%Name)),              &
          TRIM(ADJUSTL(sol_calfit_he5fields(i)%Dimensions)),   " ", &
          sol_calfit_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO

  DO i = 1, n_radcal_fields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(rad_calfit_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )
     END IF

     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( rad_calfit_he5fields(i), errstat )

     rad_calfit_he5fields(i)%Swath_ID = pge_swath_id
     errstat = HE5_SWdefdfld (                                      &
          rad_calfit_he5fields(i)%Swath_ID,                         &
          TRIM(ADJUSTL(rad_calfit_he5fields(i)%Name)),              &
          TRIM(ADJUSTL(rad_calfit_he5fields(i)%Dimensions)),   " ", &
          rad_calfit_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO

  DO i = 1, n_radref_fields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(rad_reffit_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )
     END IF

     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( rad_reffit_he5fields(i), errstat )

     rad_reffit_he5fields(i)%Swath_ID = pge_swath_id
     errstat = HE5_SWdefdfld (                                      &
          rad_reffit_he5fields(i)%Swath_ID,                         &
          TRIM(ADJUSTL(rad_reffit_he5fields(i)%Name)),              &
          TRIM(ADJUSTL(rad_reffit_he5fields(i)%Dimensions)),   " ", &
          rad_reffit_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO


  ! --------------------
  ! Radiance data fields
  ! --------------------
  ! ----------------------
  ! (1) Common Data Fields
  ! ----------------------
  DO i = 1, n_cdfields

     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(comdata_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )

     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )! comp_par )
     END IF

     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( comdata_he5fields(i), errstat )

     errstat = HE5_SWdefdfld (                                   &
          comdata_he5fields(i)%Swath_ID,                         &
          TRIM(ADJUSTL(comdata_he5fields(i)%Name)),              &
          TRIM(ADJUSTL(comdata_he5fields(i)%Dimensions)),        &
          " ",                                                   &
          comdata_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )

  END DO

  ! --------------------------
  ! (2) Diagnostic Data Fields
  ! --------------------------
  IF ( ctrvar%yn_diagnostic_run ) THEN

     DO i = 1, n_diag_fields

    	! Is this field outputted? CCM
        IF ( yn_output_diag(i) ) THEN
 
           CALL he5_check_for_compressibility ( &
                nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(diagnostic_he5fields(i)%Dimensions)), &
                yn_compress_field, n_chunk_dim, chunk_dim )
           IF ( yn_compress_field ) THEN
              errstat = HE5_SWdefcomch ( &
                   pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
           ELSE
              errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
              errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )! comp_par )
           END IF
           
           ! -------------
           ! Set FillValue
           ! -------------
           CALL he5_set_fill_value ( diagnostic_he5fields(i), errstat )
           
           errstat = HE5_SWdefdfld (                                      &
                diagnostic_he5fields(i)%Swath_ID,                         &
                TRIM(ADJUSTL(diagnostic_he5fields(i)%Name)),              &
                TRIM(ADJUSTL(diagnostic_he5fields(i)%Dimensions)),        &
                " ",                                                      &
                diagnostic_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
        END IF
     END DO
  END IF
  
  ! --------------
  ! AMF Cloud Data
  ! --------------
  DO i = 1, n_voc_fields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(voc_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, he5_nocomp_par )! comp_par )
     END IF
     
     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( voc_he5fields(i), errstat )
     
     errstat = HE5_SWdefdfld (                           &
          voc_he5fields(i)%Swath_ID,                     &
          TRIM(ADJUSTL(voc_he5fields(i)%Name)),          &
          TRIM(ADJUSTL(voc_he5fields(i)%Dimensions)),    &
          " ",                                           &
          voc_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO

  ! -------------------------------------------------
  ! (4) Scattering weights, gas profile averaging kernels
  ! and albedo gga
  ! -------------------------------------------------
  DO i = 1, n_sw_fields
     CALL he5_check_for_compressibility ( &
          nTimes, nXtrack, nSwLevels, TRIM(ADJUSTL(sw_he5fields(i)%Dimensions)), &
          yn_compress_field, n_chunk_dim, chunk_dim )
     IF ( yn_compress_field ) THEN
        errstat = HE5_SWdefcomch ( &
             pge_swath_id, he5_comp_type, he5_comp_par, n_chunk_dim, &
             chunk_dim(1:n_chunk_dim) )
     ELSE
        errstat = HE5_SWdefchunk( pge_swath_id, n_chunk_dim, &
             chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp ( pge_swath_id, he5_nocomp_type, &
             he5_nocomp_par )! comp_par )
     END IF
        
     ! -------------
     ! Set FillValue
     ! -------------
     CALL he5_set_fill_value ( sw_he5fields(i), errstat )
     
     errstat = HE5_SWdefdfld (                           &
          sw_he5fields(i)%Swath_ID,                      &
          TRIM(ADJUSTL(sw_he5fields(i)%Name)),           &
          TRIM(ADJUSTL(sw_he5fields(i)%Dimensions)),     &
          " ",                                           &
          sw_he5fields(i)%HE5_DataType, he5_hdfe_nomerge )
  END DO

  ! ------------------------------------------
  ! Check error status of swath initialization
  ! ------------------------------------------
  CALL error_check ( &
       errstat, he5_stat_ok, pge_errstat_error, OMSAO_E_HE5SWDEFFLD, &
       modulename, vb_lev_default, he5stat )
  
  ! -------------------------------------------------------------------------------
  ! Detach from and re-attach to created swath (recommended before adding to swath)
  ! -------------------------------------------------------------------------------
  errstat       = HE5_SWdetach ( pge_swath_id )
  pge_swath_id  = HE5_SWattach ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( pge_swath_id == he5_stat_fail ) CALL error_check ( &
       0, 1, pge_errstat_fatal, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )

  RETURN
END FUNCTION he5_define_fields


SUBROUTINE he5_write_wavcal_output ( nXtloc, fpix, lpix, errstat )

  USE OMSAO_indices_module, ONLY: max_calfit_idx
  USE OMSAO_he5_module
  USE OMSAO_he5_datafields_module, ONLY: sol_calfit_he5fields, &
       rad_calfit_he5fields, rad_reffit_he5fields
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, &
       he5_stat_ok, error_check, vb_lev_default, &
       omsao_e_he5swwrfld
  USE OMSAO_data_module,   ONLY: &
       n_roff_dig,                                            &
       solcal_xflag,  radcal_xflag, radref_xflag, &
       solcal_pars,   radcal_pars,  radref_pars,  &
       radref_col,    radref_dcol,  radref_rms,   &
       radref_xtrcol
  USE OMSAO_variables_module, ONLY: ctrvar

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_wavcal_output'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nXtloc, fpix, lpix

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)                             :: locerrstat, i, j, npix
  REAL    (KIND=r8), DIMENSION (max_calfit_idx) :: tmpr8
  locerrstat = pge_errstat_ok

  ! -----------------------------------------------
  ! Number of cross-track pixels actually processed
  ! -----------------------------------------------
  npix = lpix - fpix + 1

  ! -------------------------------------------------------------------------------------------
  ! Write results for solar and radiance wavelength calibration, and the radiance reference fit
  ! -------------------------------------------------------------------------------------------
  DO j = fpix, lpix
     tmpr8(1:max_calfit_idx) = solcal_pars(1:max_calfit_idx,j)
     CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
     solcal_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)

     tmpr8(1:max_calfit_idx) = radcal_pars(1:max_calfit_idx,j)
     CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
     radcal_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)

     tmpr8(1:max_calfit_idx) = radref_pars(1:max_calfit_idx,j)
     CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
     radref_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)
  END DO

  ! --------------------------------------------------
  ! There is only one variable in the data field array
  ! --------------------------------------------------
  i = 1

  ! ------------------------------------------
  ! The remaining parameters are 1-dimensional
  ! ------------------------------------------
  he5_start_2d  = (/                0, 0 /)
  he5_stride_2d = (/                1, 0 /)
  he5_edge_2d   = (/ INT(nXtloc,KIND=C_LONG), 0_i8 /)

  ! -----------------------------
  ! Exit value of Fitting Routine
  ! -----------------------------
  locerrstat = HE5_SWWRFLD (  &
       sol_calfit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(swccf_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, solcal_xflag(1:nXtloc) )
  locerrstat = HE5_SWWRFLD (  &
       rad_calfit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rwccf_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, radcal_xflag(1:nXtloc) )
  locerrstat = HE5_SWWRFLD (  &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrcf_field )), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, radref_xflag(1:nXtloc) )

  ! --------------------------------------------------------------------------
  ! Special output for radiance wavelength calibration, and the possibility of
  ! the usage of a radiance reference spectrum.
  ! --------------------------------------------------------------------------
  he5_start_2d  = (/ 0, 0 /) 
  he5_stride_2d = (/ 1, 0 /)
  he5_edge_2d   = (/ 2, 0 /)

  locerrstat = HE5_SWWRFLD (                                                  &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrlr_field )),          &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), ctrvar%radref_latrange(1:2) )

  ! ---------------------------------------------------------------------------
  ! Write results for column and column uncertainty from radiance reference fit
  ! ---------------------------------------------------------------------------
  CALL roundoff_1darr_r8 ( n_roff_dig, npix, radref_col   (fpix:lpix) )
  CALL roundoff_1darr_r8 ( n_roff_dig, npix, radref_dcol  (fpix:lpix) )
  CALL roundoff_1darr_r8 ( n_roff_dig, npix, radref_rms   (fpix:lpix) )
  CALL roundoff_1darr_r8 ( n_roff_dig, npix, radref_xtrcol(fpix:lpix) )

  he5_start_2d  = (/ 0, 0 /) 
  he5_stride_2d = (/ 1, 0 /)
  he5_edge_2d   = (/ INT(nXtloc,KIND=C_LONG), INT(0, KIND=C_LONG) /)
  locerrstat = HE5_SWWRFLD ( &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrcol_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), radref_col(1:nXtloc) )
  locerrstat = HE5_SWWRFLD ( &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrdcol_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), radref_dcol(1:nXtloc) )
  locerrstat = HE5_SWWRFLD ( &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrxcol_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), radref_xtrcol(1:nXtloc) )
  locerrstat = HE5_SWWRFLD ( &
       rad_reffit_he5fields(i)%Swath_ID, TRIM(ADJUSTL(rrrms_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), radref_rms(1:nXtloc) )

  ! ------------------
  ! Check error status
  ! ------------------

  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, errstat )

  RETURN
END SUBROUTINE he5_write_wavcal_output


SUBROUTINE he5_write_radfit_output ( &
     pge_idx, iline, nXtrack, fpix, lpix, &
     all_fitted_columns, all_fitted_errors, correlation_columns,&
     omi_fitspc, errstat )

  USE OMSAO_variables_module,  ONLY: n_fitvar_rad, n_rad_wvl
  USE OMSAO_indices_module,    ONLY: sao_molecule_names, &
       corr_didx,  corrcol_didx, correrr_didx, itnum_didx,  &
       fitwt_didx, posobs_didx,  spcobs_didx,  spcfit_didx, &
       spcres_didx
  USE OMSAO_variables_module,  ONLY: ctrvar
  USE OMSAO_data_module, ONLY: nxtrack_max, n_roff_dig, nwavel_max, &
       itnum_flag
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, he5_stat_ok, &
       error_check, vb_lev_default, omsao_e_he5swwrfld

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_radfit_output'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, iline, nXtrack, fpix, lpix

  ! CCM 
  REAL (KIND=r8), INTENT (IN), DIMENSION(nwavel_max,nxtrack_max,4) :: omi_fitspc
  
  ! Temporary array 
  REAL (KIND=r8),DIMENSION(nwavel_max,nxtrack_max) :: tmp_fitspc

  ! ----------------------
  ! Modified variables gga
  ! --------------------------------------------------------------
  ! The reason these appear as INTENT (INOUT) is that the rounding
  ! will modify the values. Hence INTENT (IN) produces a conflict.
  ! --------------------------------------------------------------
  REAL (KIND=r8), INTENT (INOUT), DIMENSION (1:n_fitvar_rad,1:nXtrack) :: &
       all_fitted_columns, all_fitted_errors, correlation_columns

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! -------------------
  ! Local variables gga
  ! -------------------
  INTEGER   (KIND=i4)                      :: locerrstat, i, npix
  CHARACTER (LEN=4)                        :: molstr

  ! CCM loop indices
  INTEGER (KIND=i4) :: ii,jj

  locerrstat = pge_errstat_ok

  ! -------------------------------------------------------
  ! Some fields have the PGE molecule name appended to them
  ! -------------------------------------------------------
  molstr = sao_molecule_names(pge_idx)

  ! ---------------------------------------------------
  ! Write current data block fitting output to HE5 file
  ! ---------------------------------------------------

  ! -----------------------------------------------
  ! Number of cross-track pixels actually processed
  ! -----------------------------------------------
  npix = lpix - fpix + 1

  ! --------------------------------------------------
  ! Correlation Information (requires additional rank)
  ! --------------------------------------------------
  IF ( ctrvar%yn_diagnostic_run ) THEN

     he5_start_2d = (/ 0, iline /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nxtrack, 1 /)
     
     IF( yn_output_diag( itnum_didx ) ) THEN 
        locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(itnum_field)), &
             he5_start_2d, he5_stride_2d, he5_edge_2d, itnum_flag(1:nXtrack,iline) )
     ENDIF
   
     DO i = 1, n_fitvar_rad
        CALL roundoff_1darr_r8 ( n_roff_dig, npix, correlation_columns(i,fpix:lpix) )
        CALL roundoff_1darr_r8 ( n_roff_dig, npix, all_fitted_columns (i,fpix:lpix) )
        CALL roundoff_1darr_r8 ( n_roff_dig, npix, all_fitted_errors  (i,fpix:lpix) )
     END DO
     
     he5_start_3d  = (/ 0,            0      , iline /)
     he5_stride_3d = (/ 1,            1      ,     1 /)
     he5_edge_3d   = (/ n_fitvar_rad, nXtrack,     1 /)
     
     IF( yn_output_diag( corr_didx ) ) THEN
        locerrstat = HE5_SWWRFLD ( pge_swath_id, corr_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             correlation_columns(1:n_fitvar_rad,1:nXtrack) )
     ENDIF
     
     IF( yn_output_diag( corrcol_didx ) ) THEN
        locerrstat = HE5_SWWRFLD ( pge_swath_id, corrcol_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             all_fitted_columns(1:n_fitvar_rad,1:nXtrack) )
     ENDIF
     
     IF( yn_output_diag( correrr_didx ) ) THEN
        locerrstat = HE5_SWWRFLD ( pge_swath_id, correrr_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             all_fitted_errors(1:n_fitvar_rad,1:nXtrack) )
     ENDIF
     
  END IF
  
  ! -------------------------------
  ! CCM Write Fit residuals to disk
  ! -------------------------------
  IF( ctrvar%yn_diagnostic_run ) THEN
          
     ! Write to he5
     he5_start_3d  = (/ 0,                0, iline /)
     he5_stride_3d = (/ 1,                1,     1 /)
     he5_edge_3d   = (/ n_rad_wvl,  nXtrack,     1 /)
     
     ! (1) Model Spectrum
     IF( yn_output_diag( spcfit_didx ) ) THEN
        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcfit_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             omi_fitspc(1:n_rad_wvl,1:nXtrack,1) )
     ENDIF
     
     ! (2) Measured Spectrum
     IF( yn_output_diag( spcobs_didx ) ) THEN
        DO ii=1,n_rad_wvl
           DO jj=1,nXtrack
                 tmp_fitspc(ii,jj) = omi_fitspc(ii,jj,2)
           ENDDO
        ENDDO
        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcobs_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             tmp_fitspc(1:n_rad_wvl,1:nXtrack) )
     ENDIF
     
     ! (3) Measured Position
     IF( yn_output_diag( posobs_didx ) ) THEN
        DO ii=1,n_rad_wvl
           DO jj=1,nXtrack
                 tmp_fitspc(ii,jj) = omi_fitspc(ii,jj,3)
           ENDDO
           
           CALL roundoff_1darr_r8 ( n_roff_dig, nXtrack, tmp_fitspc(ii,1:nXtrack) )
           
        ENDDO
        
        locerrstat = HE5_SWWRFLD ( pge_swath_id,posobs_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             tmp_fitspc(1:n_rad_wvl,1:nXtrack) )
     ENDIF
     
     ! (4) Fit Weights
     IF( yn_output_diag( fitwt_didx ) ) THEN
        DO ii=1,n_rad_wvl
           DO jj=1,nXtrack
                 tmp_fitspc(ii,jj) = omi_fitspc(ii,jj,4)
           ENDDO
        ENDDO
        
        locerrstat = HE5_SWWRFLD ( pge_swath_id,fitwt_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
          tmp_fitspc(1:n_rad_wvl,1:nXtrack) )
     ENDIF
  
     ! Residual Spectrum 
     IF( yn_output_diag( spcres_didx ) ) THEN 
        DO ii=1,n_rad_wvl
           DO jj=1,nXtrack
              tmp_fitspc(ii,jj) = omi_fitspc(ii,jj,2) - omi_fitspc(ii,jj,1)
           ENDDO
        ENDDO
        
        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcres_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
             tmp_fitspc(1:n_rad_wvl,1:nXtrack) )
     ENDIF
     
  END IF
  ! CCM
  
  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, errstat )
  
  RETURN
END SUBROUTINE he5_write_radfit_output


SUBROUTINE he5_write_common_mode ( nXtrack, npts, errstat )

  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, he5_stat_ok, &
       error_check, vb_lev_default, omsao_e_he5swwrfld
  USE OMSAO_he5_module
  USE OMSAO_indices_module,   ONLY: commcnt_didx, commspc_didx, &
       commwvl_didx, ccdpix_didx
  USE OMSAO_variables_module, ONLY: common_mode_spec, ctrvar
  USE OMSAO_data_module,   ONLY: n_roff_dig

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_write_common_mode'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nXtrack, npts

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                            :: locerrstat, j1, j2
  INTEGER   (KIND=i2), DIMENSION (nXtrack,2)     :: locccd
  INTEGER   (KIND=i4), DIMENSION (nXtrack)       :: loccnt
  REAL      (KIND=r4), DIMENSION (nXtrack, npts) :: locwvl
  REAL      (KIND=r8), DIMENSION (nXtrack, npts) :: locspc

  locerrstat = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  ! Geolocation Fields: Latitiude, Longitude, Solar Zenith, Viewing Zenith, Time
  ! ----------------------------------------------------------------------------
  ! NOTE: The Aura Altitude and TIME fields are one-dimensional, and so require
  !       a different stride for writing. We write those ones first, then set
  !       the strides for the rest of the fields.
  ! ----------------------------------------------------------------------------
  j1                       =       1
  j2                       = nXtrack
  locccd(1:nXtrack,1:2)    =        common_mode_spec%CCDPixel    (j1:j2,1:2)
  loccnt(1:nXtrack)        =        common_mode_spec%RefSpecCount(j1:j2)
  locspc(1:nXtrack,1:npts) =        common_mode_spec%RefSpecData (j1:j2,1:npts)
  locwvl(1:nXtrack,1:npts) = REAL ( common_mode_spec%RefSpecWavs (j1:j2,1:npts), KIND=r4 )

  !DO j1 = 1, npts
  !   WRITE (22,'(0PF10.4, 1PE15.5)') locwvl(10,j1), locspc(10,j1)
  !END DO


  CALL roundoff_2darr_r4 ( n_roff_dig, nXtrack, npts, locwvl (1:nXtrack,1:npts) )
  CALL roundoff_2darr_r8 ( n_roff_dig, nXtrack, npts, locspc (1:nXtrack,1:npts) )


  ! ------------------------------------------
  ! Common Mode Spectrum - Count for Averaging
  ! ------------------------------------------
  IF( yn_output_diag(commcnt_didx) ) THEN
     he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ nXtrack, 0 /)
     locerrstat = HE5_SWWRFLD ( &
          pge_swath_id, commcnt_field, he5_start_2d, he5_stride_2d, he5_edge_2d, loccnt(1:nXtrack) )
  ENDIF
  
  ! ----------------------------------------------
  ! Common Mode Spectrum - Wavelengths and Spectra
  ! ----------------------------------------------
  he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, npts /)
  
  IF( yn_output_diag(commwvl_didx) ) THEN
     locerrstat = HE5_SWWRFLD ( &
          pge_swath_id, commwvl_field, he5_start_2d, he5_stride_2d, he5_edge_2d, locwvl(1:nXtrack,1:npts) )
  ENDIF
  IF( yn_output_diag(commspc_didx) ) THEN
     locerrstat = HE5_SWWRFLD ( &
          pge_swath_id, commspc_field, he5_start_2d, he5_stride_2d, he5_edge_2d, locspc(1:nXtrack,1:npts) )
  ENDIF
  ! --------------------------
  ! CCD Pixel - First and Last
  ! --------------------------
  IF ( ctrvar%yn_diagnostic_run .AND. yn_output_diag(ccdpix_didx)) THEN
     he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, 2 /)
     locerrstat = HE5_SWWRFLD ( &
          pge_swath_id, ccdpix_field, he5_start_2d, he5_stride_2d, he5_edge_2d, locccd(1:nXtrack,1:2) )
  END IF
  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, errstat )

  RETURN
END SUBROUTINE he5_write_common_mode

SUBROUTINE he5_write_ins_database ( database_he5, database_he5_wvl, nRefSpec,  npts, nXtrack, errstat )

  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, he5_stat_ok, &
       error_check, vb_lev_default, omsao_e_he5swwrfld
  USE OMSAO_he5_module
  USE OMSAO_variables_module, ONLY: refspecs_original
  USE OMSAO_indices_module,   ONLY: spdata_didx, spnrmf_didx, spdatw_didx
  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_ins_database'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN)                                   :: nXtrack, npts, nRefSpec
  REAL    (KIND=r8), INTENT (IN), DIMENSION(nRefSpec,npts,nXtrack) :: database_he5
  REAL    (KIND=r8), INTENT (IN), DIMENSION(npts,nXtrack)          :: database_he5_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                 :: locerrstat,ii
  
  REAL (KIND=r8), DIMENSION(nRefSpec) :: tmp_normfactor
  
  ! =============================================================================
  ! he5_write_ins_database starts here
  ! =============================================================================
  
  ! Set local error status
  locerrstat = pge_errstat_ok
	
  ! Datablock
  he5_start_3d  = (/ 0,            0,       0 /)
  he5_stride_3d = (/ 1,            1,       1 /)
  he5_edge_3d   = (/ nRefSpec,  npts, nXtrack /)
	
  ! Write refspec database
  IF( yn_output_diag(spdata_didx) ) THEN
     locerrstat = HE5_SWWRFLD ( pge_swath_id,spdata_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                database_he5 )
  ENDIF

  ! Datablock
  he5_start_2d  = (/    0,       0 /)
  he5_stride_2d = (/    1,       1 /)
  he5_edge_2d   = (/ npts, nXtrack /)  
  ! Write refspec database wavelength
  IF( yn_output_diag(spdatw_didx) ) THEN
     locerrstat = HE5_SWWRFLD ( pge_swath_id, spdatw_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
                                database_he5_wvl )
  ENDIF


  DO ii=1,nRefSpec
     tmp_normfactor(ii) = refspecs_original(ii)%NormFactor
  END DO
  
  ! Datablock
  he5_start_2d  = (/ 0, 0 /)
  he5_stride_2d = (/ 1, 0 /)
  he5_edge_2d = (/ nRefSpec, 0 /)
  
  ! Write Normalisation factors
  IF( yn_output_diag( spnrmf_didx ) ) THEN
     locerrstat = HE5_SWWRFLD ( pge_swath_id, spnrmf_field, he5_start_2d, he5_stride_2d, he5_edge_2d,&
          tmp_normfactor)
  ENDIF
  
  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, errstat )
	
  RETURN
  
END SUBROUTINE he5_write_ins_database

SUBROUTINE he5_write_fitting_statistics ( &
     pge_idx, nx, nt, saomqf, avg_col, avg_dcol, avg_rms, errstat )

  USE OMSAO_indices_module, ONLY: sao_molecule_names, correlm_didx
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, he5_stat_ok, &
       omsao_e_he5swwrfld, error_check, vb_lev_default
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_data_module, ONLY: correlation_names_concat, nclenfit, nlines_max

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=28), PARAMETER :: modulename = 'he5_write_fitting_statistics'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nx, nt
  REAL    (KIND=r8),                          INTENT (IN) :: avg_col, avg_dcol, avg_rms
  INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saomqf

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=C_LONG), PARAMETER :: zerocl = 0, onecl = 1 
  CHARACTER (LEN=4)        :: molstr
  INTEGER   (KIND=i4)      :: locerrstat, iline, nt_loop

  ! ---------------------------------------------------
  ! Write current data block fitting output to HE5 file
  ! ---------------------------------------------------

  locerrstat = pge_errstat_ok

  ! ----------------------------------------------------------
  ! Average fields have the PGE molecule name appended to them
  ! ----------------------------------------------------------
  molstr = sao_molecule_names(pge_idx)


  ! ----------------------------------------
  ! Write the Main Data Quality Flag to file
  ! ----------------------------------------
  DO iline = 0, nt-1, nlines_max

     ! --------------------------------------------------------
     ! Check if loop ends before n_times_loop max is exhausted.
     ! --------------------------------------------------------
     nt_loop = MIN ( nlines_max, nt-iline )

     ! --------------------------------------------
     ! Set Start, Stride, and Edge of writing block
     ! --------------------------------------------
     he5_start_2d  = (/            zerocl,   INT(iline,KIND=C_LONG) /)
     he5_stride_2d = (/            onecl,                 onecl /)
     he5_edge_2d   = (/ INT(nx,KIND=C_LONG), INT(nt_loop,KIND=C_LONG) /)
    
     ! ---------------
     ! Write MDQA flag
     ! ---------------
     locerrstat = HE5_SWwrfld ( pge_swath_id, TRIM(ADJUSTL(mainqa_field)),           &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saomqf(1:nx,iline:iline+nt_loop-1) )

  END DO


  ! ----------------------------------------------
  ! All quantities written below are single numbers
  ! ----------------------------------------------
  he5_start_2d  = (/ zerocl, zerocl /)
  he5_stride_2d = (/  onecl, zerocl /)
  he5_edge_2d   = (/  onecl, zerocl /)

  ! ---------------------------------------------------
  ! Average column amount, uncertainty, and fitting RMS
  ! ---------------------------------------------------
  locerrstat = HE5_SWwrfld ( &
       pge_swath_id,  TRIM(ADJUSTL(avgcol_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_col  )
  locerrstat = HE5_SWwrfld ( &
       pge_swath_id, TRIM(ADJUSTL(avgdcol_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_dcol )
  locerrstat = HE5_SWwrfld ( &
       pge_swath_id,  TRIM(ADJUSTL(avgrms_field)), &
       he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_rms  )

  ! -------------------------------------------
  ! Elements included in the Correlation Output
  ! -------------------------------------------
  ! This was added later, after the routine had been named to include AVERAGE.
  ! It is a bit confusing, but no matter where we places it at this point, it
  ! will be. And this is one of the few places outside the scan line loop, so
  ! the field won't be written over and over again.
  ! ----------------------------------------------------------------------------
  ! To make matters worse, STRINGS cannot be written to the swath (at least on
  ! non-TLCF implementations of HDF-EOS5). We have thus adopted the work-around
  ! solution of converting everything to INTEGERs first and write those to file.
  ! ----------------------------------------------------------------------------
  IF ( ctrvar%yn_diagnostic_run .AND. yn_output_diag(correlm_didx) ) THEN
     he5_start_2d  = (/ zerocl, zerocl /)
     he5_stride_2d = (/  onecl, zerocl /)
     he5_edge_2d   = (/ INT(nclenfit,KIND=C_LONG), zerocl /)
     locerrstat = HE5_SWWRFLD ( &
          pge_swath_id, correlm_field,  he5_start_2d, he5_stride_2d, he5_edge_2d, &
          TRIM(ADJUSTL(correlation_names_concat)) )
  END IF
  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, errstat )

  RETURN
END SUBROUTINE he5_write_fitting_statistics

FUNCTION he5_set_field_attributes ( pge_idx ) RESULT ( he5stat )

  !----------------------------------------------------------------------
  ! This function sets to FillValue all those entries in the swath fields
  ! that haven't been initialized properly at this point.
  ! 
  ! It is called after all processing is complete.
  ! ---------------------------------------------------------------------
  !
  ! Input:
  !   pge_idx    - Index for current PGE
  !
  ! Return: he5stat = OMI_E_SUCCESS if it is
  !
  ! Variables passed through MODULE:
  !   pge_swath_file_id - id number for HE5 output file (required for closing it)
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_indices_module, ONLY: sao_molecule_names
  USE OMSAO_variables_module, ONLY: ctrvar
  USE OMSAO_he5_module
  USE OMSAO_he5_datafields_module, ONLY: sw_he5fields, voc_he5fields, &
       diagnostic_he5fields, comdata_he5fields, rad_reffit_he5fields, &
       rad_calfit_he5fields, sol_calfit_he5fields, geo_he5fields
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, &
       vb_lev_default, he5_stat_ok, omsao_w_he5swrlattr, error_check

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'he5_set_field_attributes'
 
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4),      INTENT(IN) :: pge_idx

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4) :: i, locerrstat
  CHARACTER (LEN=4)   :: molstr

  he5stat    = pge_errstat_ok
  locerrstat = pge_errstat_ok

  ! ---------------------------------------------------
  ! Some fields have PGE molecule name appended to them
  ! ---------------------------------------------------
  molstr = sao_molecule_names (pge_idx)

  ! --------------------------------------------------------------------------
  ! Set the Title, Units, and MissingValue for all Geolocation and Data Fields.
  ! This is a non-critical operation, so for now we don't care whether the 
  ! error status might be "not O.K."
  !
  ! This is a lot of repeated code, so eventually we want to come up with
  ! something more elegant. But for now it will have to do.
  ! --------------------------------------------------------------------------


  ! ------------------
  ! Geolocation Fields
  ! ------------------
  DO i = 1, n_gfields
     CALL he5_write_local_attributes ( "", geo_he5fields(i), locerrstat )
  END DO

  ! ---------------------------------------------------------------------
  ! Data fields associated with solar and radiance wavelength calibration
  ! ---------------------------------------------------------------------
  DO i = 1, n_solcal_fields
     CALL he5_write_local_attributes ( "", sol_calfit_he5fields(i), locerrstat )
  END DO
  DO i = 1, n_radcal_fields
     CALL he5_write_local_attributes ( "", rad_calfit_he5fields(i), locerrstat )
  END DO
  DO i = 1, n_radref_fields
     CALL he5_write_local_attributes ( "", rad_reffit_he5fields(i), locerrstat )
  END DO

  ! ----------------
  ! Main Data Fields
  ! ----------------
  DO i = 1, n_cdfields
     CALL he5_write_local_attributes ( "", comdata_he5fields(i), locerrstat )
  END DO

  ! ----------------------
  ! Diagnostic Data Fields
  ! ----------------------
  IF ( ctrvar%yn_diagnostic_run ) THEN
     DO i = 1, n_diag_fields      
        ! Only do if we are outputting the field CCM
        IF( yn_output_diag(i) ) THEN
           CALL he5_write_local_attributes ( "", diagnostic_he5fields(i), locerrstat )
        ENDIF
     END DO
  END IF
  
  ! ---------------------
  ! AMF cloud data fields
  ! ---------------------
  DO i = 1, n_voc_fields
     CALL he5_write_local_attributes ( "", voc_he5fields(i), locerrstat )
  END DO
  
  ! --------------------------------------------------
  ! Scattering weights, gas profile, averaging kernels
  ! and albedo attributes. gga
  ! --------------------------------------------------
  DO i = 1, n_sw_fields
     CALL he5_write_local_attributes ( "", sw_he5fields(i), locerrstat )
  END DO


  ! -----------------------------------------------------------------------
  ! Assign the error status returned from this routine. If we have
  ! reached here, then everything has gone well and we can exit with "O.K."
  ! -----------------------------------------------------------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5SWRLATTR, &
       modulename, vb_lev_default, he5stat )

  RETURN
END FUNCTION he5_set_field_attributes


FUNCTION he5_write_global_attributes ( ) RESULT ( he5stat )

  !----------------------------------------------------------------------
  ! This function writes HDF Global Attributes defined by the PGE.
  ! It uses the HE5_EHwrtglatt function, which writes the attribute,
  ! and defines it beforehand should it not exists.
  !
  ! It is called after all processing is complete.
  ! ---------------------------------------------------------------------
  !
  ! Input:
  !   none
  !
  ! Return: he5stat = OMI_S_SUCCESS if it is
  !
  ! Variables passed through MODULE:
  !   pge_swath_file_id - id number for HE5 output file (required for closing it)
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_he5_module
  USE OMSAO_metadata_module, ONLY: n_mdata_dbl, n_mdata_str, mdata_string_fields, &
       mdata_string_values, mdata_double_fields, mdata_double_values
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, error_check, pge_errstat_warning, &
       omsao_w_he5ehwrglatt, f_sep, he5_stat_ok, pgs_s_success, vb_lev_default, &
       omsao_w_mdl2arc
  USE OMSAO_parameters_module, ONLY: maxchlen, n_fit_winwav
  USE OMSAO_indices_module, ONLY: &
       n_config_luns, yn_config_lun_autocopy, config_lun_strings, config_lun_values
  USE OMSAO_variables_module, ONLY: ctrvar

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=27), PARAMETER :: modulename = 'he5_write_global_attributes'
 
  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=C_LONG), PARAMETER ::  onecl = 1
  INTEGER   (KIND=i4)      :: locerr, i, imd
  INTEGER   (KIND=C_LONG)     :: nlen48
  CHARACTER (LEN=maxchlen) :: parname
  REAL      (KIND=r4), DIMENSION (n_fit_winwav+2) :: fitwinlim

  locerr  = pge_errstat_ok
  he5stat = pge_errstat_ok

  ! -------------------------------------------------------------------
  ! First Global Attributes required by the Aura File Format Guidelines
  ! -------------------------------------------------------------------
  ! The first loop auto-copies some strings to the L2 output file
  ! -------------------------------------------------------------
  DO i = 1, n_config_luns
     IF ( yn_config_lun_autocopy(i) ) THEN
        parname = TRIM(ADJUSTL(config_lun_strings(i)))
        nlen48 = INT( LEN_TRIM(ADJUSTL(config_lun_values(i))), KIND=C_LONG )
        locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
             HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(config_lun_values(i))) )
        CALL error_check ( &
             locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
             modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
     END IF
  END DO

  parname = "InputVersions"
  nlen48 = INT ( LEN_TRIM(ADJUSTL(input_versions)), KIND=C_LONG )
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(input_versions)) )

  parname = "OrbitData"
  nlen48 = INT ( LEN_TRIM(ADJUSTL(l1b_orbitdata)), KIND=C_LONG )
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(l1b_orbitdata))           )

  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
  parname = "GranuleDay"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_day )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
  parname = "GranuleMonth"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_month )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
  parname = "GranuleYear"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_year )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  ! -------------------------------------
  ! Finally the TAI time of 0z of Granule
  ! -------------------------------------
  parname = tai_attr
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_DOUBLE, onecl, TAI93At0zOfGranule )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  ! ------------------------------------------------------------------------
  ! Now the Global Attributes that are set by the PGE: PGE statistics mostly
  ! ------------------------------------------------------------------------
  parname = "FittingWindowLimits"
  fitwinlim(1:n_fit_winwav) = REAL( ctrvar%fit_winwav_lim(1:n_fit_winwav), KIND=r4 )
  fitwinlim(5:6)            = REAL( ctrvar%fit_winexc_lim(1:2),            KIND=r4 )
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_FLOAT, INT( n_fit_winwav+2, KIND=r8), fitwinlim(1:n_fit_winwav+2) )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )


  parname = "NumberOfCrossTrackPixels"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofCrossTrackPixels )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfScanLines"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofScanLines )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfInputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofInputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfGoodInputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofGoodInputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfGoodOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofGoodOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfSuspectOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofSuspectOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfBadOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofBadOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfConvergedSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofConvergedSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfFailedConvergenceSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofFailedConvergenceSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfExceededIterationsSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofExceededIterationsSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "NumberOfOutOfBoundsSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_INT, onecl, NrofOutofBoundsSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "PercentGoodOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_FLOAT, onecl, PercentGoodOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "PercentBadOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_FLOAT, onecl, PercentBadOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  parname = "PercentSuspectOutputSamples"
  locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
       HE5T_NATIVE_FLOAT, onecl, PercentSuspectOutputSamples )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
       modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

  ! ------------------------------------------------------------------
  ! And here we have the ArchiveMetadata, which have to be written as
  ! Global Attributes. This section somewhat mirrors the one in 
  ! SET_L2_METADATA, picking up the missing "arc" pieces. The Metadata
  ! section was written before it was known that the ArchivedMetadata
  ! could not be written as such.
  ! ------------------------------------------------------------------
  ! -----------------------------
  ! Write STRING MetaData to file
  ! -----------------------------
  DO imd = 1, n_mdata_str
     IF ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) == "arc" ) THEN
        nlen48 = INT ( LEN_TRIM(ADJUSTL(mdata_string_values(imd))), KIND=C_LONG )
        locerr = HE5_EHwrglatt ( &
             pge_swath_file_id, TRIM(ADJUSTL(mdata_string_fields(1,imd))),  &
             HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(mdata_string_values(imd))) )
        CALL error_check ( &
             locerr, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2ARC, &
             modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,imd))), &
             vb_lev_default, he5stat )
     END IF
  END DO
  ! -------------------
  ! Set DOUBLE MetaData
  ! -------------------
  DO imd = 1, n_mdata_dbl
     IF ( TRIM(ADJUSTL(mdata_double_fields(2,imd))) == "arc" ) THEN
        locerr = HE5_EHwrglatt ( &
             pge_swath_file_id, TRIM(ADJUSTL(mdata_double_fields(1,imd))),  &
             HE5T_NATIVE_DOUBLE, onecl, mdata_double_values(imd) )
        CALL error_check ( &
             locerr, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2ARC, &
             modulename//f_sep//TRIM(ADJUSTL(mdata_double_fields(1,imd))), &
             vb_lev_default, he5stat )
     END IF
  END DO

  RETURN
END FUNCTION he5_write_global_attributes

FUNCTION he5_close_output_file ( ) RESULT ( he5stat )

  !------------------------------------------------------------------------------
  ! This function detatches from the HE5 swath and closes the HE5 output file.
  !
  ! Input: NONE
  !
  !------------------------------------------------------------------------------

  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, &
       error_check, omsao_w_he5swclose, he5_stat_ok, vb_lev_default

  IMPLICIT NONE

  ! ---------------------------------------
  ! Name of this module/subroutine/function
  ! ---------------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_close_output_file'

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! --------------
  ! Local variable
  ! --------------
  INTEGER (KIND=i4) :: locerr

  he5stat = pge_errstat_ok
  ! -----------------------------------------------
  ! Detach from HE5 swath and close HE5 output file
  ! -----------------------------------------------
  locerr = HE5_SWDETACH ( pge_swath_id )
  locerr = HE5_SWCLOSE  ( pge_swath_file_id )
  CALL error_check ( &
       locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5SWCLOSE, &
       modulename, vb_lev_default, he5stat )

  ! -----------
  ! Prefit file
  ! -----------
  IF ( TRIM(ADJUSTL(prefit_swath_name)) /= 'undefined' .AND. &
       prefit_swath_id                  /= -1          .AND. &
       prefit_swath_file_id             /= -1              ) THEN
     locerr = HE5_SWDETACH ( prefit_swath_id )
     locerr = HE5_SWCLOSE  ( prefit_swath_file_id )
  END IF
  
  RETURN
END FUNCTION he5_close_output_file


SUBROUTINE he5_set_fill_value ( data_field, errstat )

  USE OMSAO_errstat_module, ONLY: pge_errstat_ok
  USE OMSAO_parameters_module,     ONLY: str_missval
  USE OMSAO_he5_datafields_module, ONLY: DataField_HE5
  USE OMSAO_he5_module

  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  TYPE (DataField_HE5), INTENT (IN) :: data_field

  ! ----------------
  ! Output Variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat
  INTEGER (KIND=i2) :: locfil_i2
  INTEGER (KIND=i4) :: locfil_i4
  REAL    (KIND=r4) :: locfil_r4

  locerrstat = pge_errstat_ok

  ! MissingValue and ValidRange depend on data type
  SELECT CASE ( data_field%HE5_DataType )
  CASE ( HE5T_NATIVE_DOUBLE )
     errstat = HE5_SWsetfill (                   &
          data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          data_field%HE5_DataType, data_field%FillValue )
  CASE ( HE5T_NATIVE_FLOAT )
     locfil_r4 = REAL (data_field%FillValue, KIND=r4)
     errstat = HE5_SWsetfill (                   &
          data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          data_field%HE5_DataType, locfil_r4 )
  CASE ( HE5T_NATIVE_INT )
     locfil_i4 = INT (data_field%FillValue, KIND=i4)
     errstat = HE5_SWsetfill (                   &
          data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          data_field%HE5_DataType, locfil_i4 )
  CASE ( HE5T_NATIVE_INT16 )
     locfil_i2 = INT (data_field%FillValue, KIND=i2)
     errstat = HE5_SWsetfill (                   &
          data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          data_field%HE5_DataType, locfil_i2 )
  CASE ( HE5T_NATIVE_CHAR )
     errstat = HE5_SWsetfill (                   &
          data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          data_field%HE5_DataType, str_missval )
  CASE DEFAULT
     ! We should never reach here
  END SELECT


  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE he5_set_fill_value



SUBROUTINE he5_write_local_attributes ( addstr, data_field, errstat )

  USE OMSAO_parameters_module,     ONLY: str_missval
  USE OMSAO_he5_datafields_module, ONLY: DataField_HE5
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok

  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  CHARACTER (LEN=*),    INTENT (IN) :: addstr
  TYPE (DataField_HE5), INTENT (IN) :: data_field

  ! ----------------
  ! Output Variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER (KIND=C_LONG), PARAMETER :: onecl = 1 , twocl = 2
  INTEGER (KIND=i4) :: locerrstat
  INTEGER (KIND=i2)                :: locmis_i2
  INTEGER (KIND=i4)                :: locmis_i4
  REAL    (KIND=r4)                :: locmis_r4
  INTEGER (KIND=i2), DIMENSION (2) :: locrov_i2
  INTEGER (KIND=i4), DIMENSION (2) :: locrov_i4
  INTEGER (KIND=C_LONG)               :: nlen48
  REAL    (KIND=r4), DIMENSION (2) :: locrov_r4

  locerrstat = pge_errstat_ok

  ! Title
  nlen48 = data_field%LenTitle
  locerrstat = HE5_SWwrlattr ( &
       data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), title_attr, &
       HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%Title)) )
  ! Units
  nlen48 = data_field%LenUnits
  locerrstat = HE5_SWwrlattr ( &
       data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), units_attr, &
       HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%Units)) )
  ! UniqueFieldDefinition
  nlen48 = data_field%LenUniqueFD
  locerrstat = HE5_SWwrlattr ( &
       data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), ufd_attr, &
       HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%UniqueFD)) )
  ! Offset and ScaleFactor
  ! (suppose to be always DoublePrecsion to avoid rounding errors and improve
  !  compressibility)
  locerrstat = HE5_SWwrlattr ( &
       data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), offset_attr, &
       HE5T_NATIVE_DOUBLE, onecl, data_field%Offset )
  locerrstat = HE5_SWwrlattr ( &
       data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), scafac_attr, &
       HE5T_NATIVE_DOUBLE, onecl, data_field%ScaleFactor )

  ! MissingValue and ValidRange depend on data type
  SELECT CASE ( data_field%HE5_DataType )
  CASE ( HE5T_NATIVE_DOUBLE )
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID,  &
          TRIM(ADJUSTL(data_field%Name)), missval_attr, &
          data_field%HE5_DataType, onecl, data_field%MissingValue )
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), valids_attr, &
          data_field%HE5_DataType, twocl, data_field%ValidRange )
  CASE ( HE5T_NATIVE_FLOAT )
     locmis_r4 = REAL (data_field%MissingValue, KIND=r4)
     locrov_r4 = REAL (data_field%ValidRange,   KIND=r4)
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_r4 )
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_r4 )
  CASE ( HE5T_NATIVE_INT )
     locmis_i4 = INT (data_field%MissingValue, KIND=i4)
     locrov_i4 = INT (data_field%ValidRange,   KIND=i4)
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_i4 )
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_i4 )
  CASE ( HE5T_NATIVE_INT16 )
     locmis_i2 = INT (data_field%MissingValue, KIND=i2)
     locrov_i2 = INT (data_field%ValidRange,   KIND=i2)
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_i2 )
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_i2 )
  CASE ( HE5T_NATIVE_CHAR )
     nlen48 = LEN(str_missval)
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
          TRIM(ADJUSTL(data_field%Name)), missval_attr, &
          HE5T_NATIVE_CHAR, nlen48, str_missval )
  CASE DEFAULT
     ! We should never reach here
  END SELECT

  ! -------------
  ! Special cases
  ! -------------
  IF ( TRIM(ADJUSTL(addstr)) == rstemp_attr ) THEN
     locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
          rstemp_attr, data_field%HE5_DataType, onecl, data_field%SpecTemp )
  END IF

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE he5_write_local_attributes


FUNCTION he5_write_swath_attributes ( pge_idx ) RESULT ( he5stat )

  !----------------------------------------------------------------------
  ! This function writes HE5 Swath Attributes defined by the PGE.
  ! It uses the HE5_SWwrtatt function, which writes the attribute,
  ! and defines it in case it does not exists.
  !
  ! It is called after all processing is complete.
  ! ---------------------------------------------------------------------
  !
  ! Input:  pge_idx
  !
  ! Return: he5stat = OMI_E_SUCCESS if it is
  !
  ! Variables passed through MODULE:
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_error, omsao_e_he5swwrattr, &
       vb_lev_default, he5_stat_ok, pge_errstat_ok, error_check
  USE OMSAO_data_module, ONLY: EarthSunDistance
  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=26), PARAMETER :: modulename = 'he5_write_swath_attributes'
 
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat
  INTEGER (KIND=C_LONG) :: count48

  ! ----------------------------------
  ! Write swath attributes to HE5 file
  ! ------------------------------------------------------------------
  ! So far there are only two: VerticalCoordinate and EarthSunDistance
  ! ------------------------------------------------------------------

  he5stat    = pge_errstat_ok
  locerrstat = pge_errstat_ok

  count48 = INT(LEN_TRIM(ADJUSTL(vertical_coordinate(pge_idx))), KIND=C_LONG)
  locerrstat = he5_swwrattr ( &
       pge_swath_id, vcoordinate_field, HE5T_NATIVE_CHAR, count48, &
       TRIM(ADJUSTL(vertical_coordinate(pge_idx))) )
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRATTR, &
       modulename, vb_lev_default, he5stat )

  count48 = 1
  locerrstat = he5_swwrattr ( &
       pge_swath_id, 'EarthSunDistance', HE5T_NATIVE_FLOAT, count48, EarthSunDistance )
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRATTR, &
       modulename, vb_lev_default, he5stat )

  RETURN
END FUNCTION he5_write_swath_attributes

SUBROUTINE he5_check_for_compressibility ( &
     nTimes, nXtrack, nSwLevels, field_dim, yn_compress_field, n_chunk_dim, chunk_dim )

  USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
  USE OMSAO_he5_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY: n_fitvar_rad
  USE OMSAO_data_module,   ONLY: nclenfit, nUTCdim, n_field_maxdim, &
       nwavel_max


  ! -----------------------------------------------------------------
  ! The purpose of this subroutine is to check whether compression
  ! should be enabled for a data field. The catch with compression
  ! is that
  !
  ! (1) it only works for 2- and higher dimensional arrays
  ! (2) it requires "chunking", i.e., the definition of the chunk
  !     of data that are being written at one time.
  !
  ! The second point causes some headache: For any fields that
  ! include the "nTimes" dimension,  we are limited by the number
  ! of swath lines that are read from the L1b file (usually 100), 
  ! and we have either 2 or 3 dimensions; other fields may be
  ! written in one go. This defines the particular way of chunking.
  !
  ! The distinction can be made by checking the dimensions of the
  ! data fields and define chunking based on the strings that go
  ! into the definition of the data field.
  !
  ! The dimensions for which compression can be enabled are:
  !
  !    "nXtrack,2"
  !    "nXtrack,nwavel_max"
  !    "nXtrack,nTimes"
  !    "nXtrack+1,nTimes+1"
  !    "nWavCalPars,nXtrack"
  !    "nUTCdim,nTimes"
  !    "nFitElements,nXtrack,nTimes"
  !    "nXtrack,nTimes,nSwLevels"
  !
  ! NO compression is possible for fields of dimension:
  !
  !    "1"
  !    "nTimes"
  !    "nXtrack"
  !    "nCharLenFitElements"
  !
  ! -----------------------------------------------------------------

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4), INTENT (IN) :: nTimes, nXtrack, nSwLevels !gga
  CHARACTER (LEN=*),   INTENT (IN) :: field_dim

  ! ----------------
  ! Output variables
  ! ----------------
  LOGICAL,                                       INTENT (OUT) :: yn_compress_field
  INTEGER (KIND=i4),                             INTENT (OUT) :: n_chunk_dim
  INTEGER (KIND=C_LONG), DIMENSION (n_field_maxdim), INTENT (OUT) :: chunk_dim

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  yn_compress_field = .FALSE.
  n_chunk_dim       = 0
  chunk_dim(1:3)    = 0

  ! -----------------------------------------------
  ! Non-compression cases: Set only chunking
  ! -----------------------------------------------
  IF ( field_dim == '1' ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = 1
  END IF
  IF ( field_dim == '2' ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = 2
  END IF
  IF ( field_dim == '4' ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = 4
  END IF
  IF ( field_dim == "nTimes" ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = INT ( nTimes, KIND=C_LONG )
  END IF
     
  IF ( field_dim == "nXtrack" ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = INT ( nXtrack, KIND=C_LONG )
  END IF
  IF ( field_dim == "nCharLenFitElements" ) THEN
     n_chunk_dim  = 1
     chunk_dim(1) = INT ( nclenfit, KIND=C_LONG )
  END IF


  ! -------------------------------------------------------------------------
  ! Now the compressible cases, one by one. Since exact matches are required,
  ! there is no harm if we miss to RETURN on a non-compressible field.
  ! -------------------------------------------------------------------------

  IF ( field_dim == "nXtrack,2" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(2,KIND=C_LONG) /)
  END IF
  IF ( field_dim == "nXtrack,nwavel_max" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(nwavel_max,KIND=C_LONG) /)
  END IF
  IF ( field_dim == "nXtrack,nTimes" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(nTimes,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nXtrack+1,nTimes+1" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack+1,KIND=C_LONG), INT(nTimes+1,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nWavCalPars,nXtrack" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(max_calfit_idx,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nUTCdim,nTimes" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = (/ INT(nUTCdim,KIND=C_LONG), INT(nTimes,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nFitElements,nXtrack,nTimes" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 3
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(n_fitvar_rad,KIND=C_LONG), INT(nXtrack,KIND=C_LONG), INT(nTimes,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nXtrack,nTimes,nLevels" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 3
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(nXtrack,KIND=C_LONG), INT(nTimes,KIND=C_LONG), INT(nSwLevels,KIND=C_LONG) /)
  END IF

	! CCM New fields
  IF ( field_dim == "nwavel_max,nXtrack,nTimes" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 3
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(nwavel_max,KIND=C_LONG), INT(nXtrack,KIND=C_LONG), INT(nTimes,KIND=C_LONG) /)
  END IF
	
  IF ( field_dim == "nRfSpec,nwavel_max,nXtrack" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 3
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(max_rs_idx,KIND=C_LONG), INT(nwavel_max,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
  END IF

  IF ( field_dim == "nwavel_max,nXtrack" ) THEN
     yn_compress_field = .TRUE.
     n_chunk_dim = 2
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(nwavel_max,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
  END IF

  
  IF ( field_dim == "nRfSpec" ) THEN
     yn_compress_field = .FALSE.
     n_chunk_dim = 1
     chunk_dim(1:n_chunk_dim) = &
          (/ INT(max_rs_idx,KIND=C_LONG) /)
  END IF

  RETURN
END SUBROUTINE he5_check_for_compressibility


FUNCTION he5_open_readwrite ( file_name, swath_name ) RESULT ( he5stat )

  !------------------------------------------------------------------------------
  ! This function initializes the HE5 output swath.
  !
  ! Input:
  !   file_name  - Name of HE5 output file
  !   swath_name - Name of swath to be created
  !
  ! Return: he5stat
  !
  ! Variables passed through MODULE:
  !   pge_swath_file_id - id number for HE5 output file (required for closing it)
  !   pge_swath_id      - id number for swath (required for writing to swath)
  !
  !------------------------------------------------------------------------------

  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: vb_lev_default, pge_errstat_fatal, error_check, &
       pge_errstat_ok, omsao_e_he5swattach, he5_stat_fail, omsao_f_he5swopen

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=18), PARAMETER :: modulename = 'he5_open_readwrite'
 
  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=maxchlen), INTENT(IN) :: file_name, swath_name

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: he5stat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4) :: errstat

  he5stat = pge_errstat_ok
  errstat = pge_errstat_ok

  ! ---------------------------------------------------------------
  ! Open HE5 output file and check PGE_SWATH_FILE_ID ( -1 if error)
  ! ---------------------------------------------------------------
  pge_swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdwr )
  IF ( pge_swath_file_id == he5_stat_fail ) THEN
     CALL error_check ( &
          0, 1, pge_errstat_fatal, OMSAO_F_HE5SWOPEN, modulename, vb_lev_default, he5stat )
     RETURN
  END IF

  ! ---------------
  ! Attach to swath
  ! ---------------
  pge_swath_id  = HE5_SWattach ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( pge_swath_id == he5_stat_fail ) CALL error_check ( &
       0, 1, pge_errstat_fatal, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )

  RETURN
END FUNCTION he5_open_readwrite

SUBROUTINE saopge_geofield_read ( &
     ntimes, nxtrack, geodata_field, geodata, errstat )

  USE OMSAO_precision_module, ONLY: i2, i4, r4
  USE OMSAO_data_module,   ONLY: nlines_max
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok
  
  IMPLICIT NONE
  
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
  CHARACTER (LEN=*), INTENT (IN) :: geodata_field

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Output variable
  ! ----------------
  REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: geodata

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: locerrstat, iline, nblock
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1) :: geodata_i2

  ! -------------------------------------
  ! Initialize output and local variables
  ! -------------------------------------
  geodata_i2 = i2_missval; geodata = r4_missval ; locerrstat = pge_errstat_ok

  ! -------------------------------------------------------
  ! Loop over all lines in the file in blocks of NLINES_MAX
  ! -------------------------------------------------------
  ScanLines: DO iline = 0, ntimes-1, nlines_max

     ! --------------------------------------------------------
     ! Check if loop ends before n_times_loop max is exhausted.
     ! --------------------------------------------------------
     nblock = MIN ( nlines_max, ntimes-iline )

     ! ----------------------------------------------------
     ! Read current data block fitting output from HE5 file
     ! ----------------------------------------------------
     he5_start_2d  = (/       0,  iline /)
     he5_stride_2d = (/       1,      1 /)
     he5_edge_2d   = (/ nxtrack, nblock /)

     ! --------------------------------------
     ! Geolocation field
     ! --------------------------------------
     IF (geodata_field .EQ. thgt_field .OR. geodata_field .EQ. extr_field) THEN
        locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(geodata_field)),         &
             he5_start_2d, he5_stride_2d, he5_edge_2d, geodata_i2(1:nxtrack,iline:iline+nblock-1) )
        geodata(1:nxtrack,iline:iline+nblock-1) = INT(geodata_i2(1:nxtrack,iline:iline+nblock-1), KIND = r4)
     ELSE
        locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(geodata_field)),         &
             he5_start_2d, he5_stride_2d, he5_edge_2d, geodata(1:nxtrack,iline:iline+nblock-1) )
     END IF
  END DO ScanLines

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE saopge_geofield_read

SUBROUTINE saopge_columninfo_read (              &
     ntimes, nxtrack, saocol, saodco, saorms, saoamf, saoqf, amfdiag, errstat )

  USE OMSAO_precision_module,  ONLY: i2, i4
  USE OMSAO_parameters_module, ONLY: r8_missval, i2_missval
  USE OMSAO_data_module,    ONLY: nlines_max
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, error_check

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Output variables
  ! ----------------
  REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saocol, saodco, saorms
  REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saoamf
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saoqf, amfdiag


  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: locerrstat, iline, nblock

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  saocol = r8_missval ; saodco  = r8_missval
  saorms = r8_missval ; saoamf  = r8_missval
  saoqf  = i2_missval ; amfdiag = i2_missval

  locerrstat = pge_errstat_ok

  ! -------------------------------------------------------
  ! Loop over all lines in the file in blocks of NLINES_MAX
  ! -------------------------------------------------------
  ScanLines: DO iline = 0, ntimes-1, nlines_max

     ! --------------------------------------------------------
     ! Check if loop ends before n_times_loop max is exhausted.
     ! --------------------------------------------------------
     nblock = MIN ( nlines_max, ntimes-iline )

     ! ----------------------------------------------------
     ! Read current data block fitting output from HE5 file
     ! ----------------------------------------------------
     he5_start_2d  = (/       0,  iline /)
     he5_stride_2d = (/       1,      1 /)
     he5_edge_2d   = (/ nxtrack, nblock /)

     ! ------------------------------------------
     ! Column amount, column uncertainty, and RMS
     ! ------------------------------------------
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(col_field)),         &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saocol(1:nxtrack,iline:iline+nblock-1) )
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(dcol_field)),        &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saodco(1:nxtrack,iline:iline+nblock-1) )
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(fitrms_field)),      &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saorms(1:nxtrack,iline:iline+nblock-1) )

     ! ---------------
     ! Air Mass Factor
     ! ---------------
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)),      &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saoamf(1:nxtrack,iline:iline+nblock-1) )

     ! -----------------------------
     ! Main quality flag convergence
     ! -----------------------------
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(mainqa_field)),      &
          he5_start_2d, he5_stride_2d, he5_edge_2d, saoqf(1:nxtrack,iline:iline+nblock-1) )

     ! -------------------
     ! AMF diagnostic flag
     ! -------------------
     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(amfdiag_field)),      &
          he5_start_2d, he5_stride_2d, he5_edge_2d, amfdiag(1:nxtrack,iline:iline+nblock-1) )
     
  END DO ScanLines

  errstat = MAX ( errstat, locerrstat )
     
  RETURN
END SUBROUTINE saopge_columninfo_read

SUBROUTINE he5_write_geolocation ( nTimes, nXtrack, &
     fpix, lpix, locerrstat)

  USE OMSAO_data_module, ONLY: spacecraft_alt, instrument_flag, &
       latitude, longitude, sazimuth, szenith, &
       vazimuth, vzenith, xtrflg, height, &
       snowicefraction
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: he5_stat_ok, omsao_e_he5swwrfld, &
       pge_errstat_ok, pge_errstat_error, vb_lev_default, error_check

  IMPLICIT NONE


  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_write_geolocation'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nTimes, nXtrack, fpix, lpix

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat

  ! --------------
  ! Local variable
  ! --------------
  INTEGER (KIND=i4) :: npix

  locerrstat = pge_errstat_ok

  ! -----------------------------------------------
  ! Number of cross-track pixels actually processed
  ! -----------------------------------------------
  npix = lpix - fpix + 1

  ! ---------------------------------------------------------------------------
  ! Geolocation Fields: Latitude, Longitude, Solar Zenith, Viewing Zenith,
  !                     Solar Azimuth, Viewing Azimuth, Altitude, XtrackQuality
  ! ---------------------------------------------------------------------------
  ! NOTE: The OMPS Altitude fields are one-dimensional, and so require
  !       a different stride for writing. We write those ones first, then set
  !       the strides for the rest of the fields.
  ! ----------------------------------------------------------------------------
  he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ nTimes, 0 /)

  locerrstat = HE5_SWWRFLD ( pge_swath_id, auraalt_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
       spacecraft_alt(0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, extr_field,      he5_start_2d, he5_stride_2d, he5_edge_2d, &
       INT(instrument_flag(0:nTimes-1), KIND=2) )

  he5_start_2d = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, nTimes /)

  locerrstat = HE5_SWWRFLD ( pge_swath_id, lat_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       latitude(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, lon_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       longitude(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, saa_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       sazimuth(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, sza_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       szenith(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, vaa_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       vazimuth(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, vza_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       vzenith(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, xtr_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       xtrflg(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, thgt_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
       height(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, sif_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       snowicefraction(1:nXtrack,0:nTimes-1) )

  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, locerrstat )

  RETURN

END SUBROUTINE he5_write_geolocation

SUBROUTINE he5_write_results ( nTimes, nXtrack, locerrstat)

  USE OMSAO_data_module, ONLY: column_amount, column_uncert, &
       fitconv_flag, fit_rms     
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: he5_stat_ok, omsao_e_he5swwrfld, &
       pge_errstat_ok, pge_errstat_error, vb_lev_default, error_check

  IMPLICIT NONE


  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=17), PARAMETER :: modulename = 'he5_write_results'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: nTimes, nXtrack

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat

  locerrstat = pge_errstat_ok

  ! ---------------------------------------------------------------------------
  ! Result Fields: Column Amount, column uncertainty, fit convergence flag and
  !                fitting RMS & solar savelenght calibration convergence flag
  !                which has a different size so I have to use 1d start,
  !                stride and edge.
  ! ---------------------------------------------------------------------------
  he5_start_2d = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, nTimes /)

  locerrstat = HE5_SWWRFLD ( pge_swath_id, col_field,    &
       he5_start_2d, he5_stride_2d, he5_edge_2d, &
       column_amount(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, dcol_field,   &
       he5_start_2d, he5_stride_2d, he5_edge_2d, &
       column_uncert(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, fitcon_field, &
       he5_start_2d, he5_stride_2d, he5_edge_2d, &
       fitconv_flag(1:nXtrack,0:nTimes-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, fitrms_field, &
       he5_start_2d, he5_stride_2d, he5_edge_2d, &
       fit_rms(1:nXtrack,0:nTimes-1) )

  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, locerrstat )


  RETURN

END SUBROUTINE he5_write_results

SUBROUTINE he5_write_solarwavcal ( nw, ip, shift, residual, locerrstat)

  USE OMSAO_data_module, ONLY: irradiance_wavl, &
       irradiance_spec, irradiance_wght, &
       solcal_xflag
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: he5_stat_ok, omsao_e_he5swwrfld, &
       pge_errstat_ok, pge_errstat_error, vb_lev_default, error_check

  IMPLICIT NONE


  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_write_solarwavcal'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN)                  :: nw, ip
  REAL    (KIND=r8), INTENT (IN)                  :: shift
  REAL    (KIND=r8), DIMENSION(1:nw), INTENT (IN) :: residual

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat

  locerrstat = pge_errstat_ok

  he5_start_2d  = (/ ip-1, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ 1, 0 /)
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swshi_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
       shift ) 
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swccf_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
       solcal_xflag(ip) )


  he5_start_2d = (/ ip-1, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ 1, nw /)
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swres_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       residual(1:nw) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swwav_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       irradiance_wavl(1:nw,ip) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swirr_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       irradiance_spec(1:nw,ip) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, swwei_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       irradiance_wght(1:nw,ip) )

  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, locerrstat )

  RETURN

END SUBROUTINE he5_write_solarwavcal

SUBROUTINE he5_write_radiancewavcal ( nw, ip, shift, residual, locerrstat)

  USE OMSAO_data_module, ONLY: radcal_xflag
  USE OMSAO_he5_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, &
       omsao_e_he5swwrfld, he5_stat_ok, vb_lev_default, error_check
  USE OMSAO_indices_module, ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_variables_module, ONLY: curr_rad_spec

  IMPLICIT NONE


  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_write_solarwavcal'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN)                  :: nw, ip
  REAL    (KIND=r8), INTENT (IN)                  :: shift
  REAL    (KIND=r8), DIMENSION(1:nw), INTENT (IN) :: residual

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat

  locerrstat = pge_errstat_ok

  he5_start_2d  = (/ ip-1, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ 1, 0 /)
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwshi_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
       shift ) 
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwccf_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
       radcal_xflag(ip) )


  he5_start_2d = (/ ip-1, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ 1, nw /)
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwres_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       residual(1:nw) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwwav_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       curr_rad_spec(wvl_idx,1:nw) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwrad_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       curr_rad_spec(spc_idx,1:nw) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, rwwei_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
       curr_rad_spec(sig_idx,1:nw) )

  ! ------------------
  ! Check error status
  ! ------------------
  CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
       modulename, vb_lev_default, locerrstat )

  RETURN

END SUBROUTINE he5_write_radiancewavcal
