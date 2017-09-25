MODULE OMSAO_prefitcol_module

  ! =================================================================
  ! This module collects all variables and subroutines required for
  ! the ingestion of prefitted columns as may be required by the
  ! OMI OMHCHO PGE. They have been collected in this module to
  ! unclutter the main OMI PGE fitting routine.
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_omidata_module,    ONLY: nxtrack_max, nlines_max
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_indices_module,    ONLY: &
       o3_t1_idx, o3_t2_idx, o3_t3_idx, bro_idx, lqh2o_idx, &
       pge_hcho_idx, pge_gly_idx
  USE OMSAO_variables_module, ONLY: ctrvar
  
  IMPLICIT NONE

  ! ------------------------------------------
  ! Total number of prefitted column variables
  ! ------------------------------------------
  INTEGER (KIND=i4) :: n_prefit_vars
  
  ! --------------------------------------------
  ! BrO prefitted column variables
  ! --------------------------------------------
  CHARACTER (LEN=maxchlen)                                    :: prefit_fname
  INTEGER   (KIND=i4)                                         :: prefit_fitidx
  REAL      (KIND=r8)                                         :: prefit_var
  REAL      (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: &
       prefit_col, prefit_dcol

CONTAINS

  SUBROUTINE init_prefit_files ( pge_idx, ntimes, nxtrack, errstat )

    USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error, &
         pge_errstat_fatal, f_sep, omi_s_success, omsao_e_prefitdim, &
         vb_lev_default, error_check
    USE OMSAO_he5_module, ONLY: &
         prefit_swath_id, prefit_swath_file_id, prefit_swath_name, &
         he5_init_input_file

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntimes, nxtrack

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, pre_ntimes, pre_nxtrack
    CHARACTER (LEN=17), PARAMETER :: modulename = 'init_prefit_files'

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    prefit_swath_id  = -1 ; prefit_swath_file_id  = -1

    ! --------------------
    ! Return if no prefits
    ! --------------------
    IF ( .NOT. ctrvar%yn_prefit(1) ) RETURN
		
    ! -------
    ! Prefits
    ! -------
    locerrstat = pge_errstat_ok
    IF ( ctrvar%yn_prefit(1) ) THEN
       CALL he5_init_input_file ( &
            prefit_fname, prefit_swath_name, prefit_swath_id, prefit_swath_file_id, &
            pre_ntimes, pre_nxtrack, errstat )
       IF ( ntimes /= pre_ntimes .OR. nxtrack /= pre_nxtrack ) THEN
          locerrstat = pge_errstat_error
          CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_PREFITDIM, &
               modulename//f_sep//"Prefit access failed.", vb_lev_default, errstat )
          ctrvar%yn_prefit = .FALSE.
       END IF
    END IF
    
    RETURN
  END SUBROUTINE init_prefit_files

 SUBROUTINE read_prefit_columns ( pge_idx, nxtrack, nloop, iline, errstat )

    USE OMSAO_parameters_module, ONLY: r8_missval
    USE OMSAO_variables_module,  ONLY: refspecs_original, ctrvar
    USE OMSAO_he5_module, ONLY: prefit_swath_id
    USE OMSAO_errstat_module, ONLY: pge_errstat_ok

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nxtrack, nloop, iline

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, i, iloop, it, nxtloc
    LOGICAL           :: yn_read_amf
    CHARACTER (LEN=12), PARAMETER :: col_str  = "ColumnAmount"
    CHARACTER (LEN=17), PARAMETER :: dcol_str = "ColumnUncertainty"
    INTEGER  (KIND=i4), PARAMETER :: lcolstr = LEN(col_str), ldcolstr = LEN(dcol_str)
     
    ! ------------------------------------------
    ! Prefitted columns and column uncertainties
    ! ------------------------------------------
    yn_read_amf = .TRUE. ; locerrstat = pge_errstat_ok
    IF ( ctrvar%yn_prefit(1) ) THEN
       CALL he5_read_prefit_columns (                                 &
            prefit_swath_id, nloop, nxtrack, iline,                   &
            lcolstr,   col_str, prefit_col (1:nxtrack,0:nloop-1), &
            ldcolstr, dcol_str, prefit_dcol(1:nxtrack,0:nloop-1), &
            yn_read_amf, locerrstat )
       errstat = MAX( errstat, locerrstat )
          
       ! -------------------------------------------------------------------
       ! Multiply columns with normalization factor to return to true values
       ! -------------------------------------------------------------------
       WHERE ( prefit_col (1:nxtrack,0:nloop-1) > r8_missval )
          prefit_col (1:nxtrack,0:nloop-1) = &
               prefit_col (1:nxtrack,0:nloop-1) * refspecs_original(ctrvar%prefit_idx)%NormFactor
          prefit_dcol(1:nxtrack,0:nloop-1) = &
               prefit_dcol(1:nxtrack,0:nloop-1) * refspecs_original(ctrvar%prefit_idx)%NormFactor
       END WHERE
    END IF
    
    ! --------------------------------------------------------------------------
    ! Shift the prefit-values to the proper index positions (e.g., spatial zoom)
    ! --------------------------------------------------------------------------
    DO iloop = 0, nloop-1

       it = iline + iloop

       ! ---------------------------------------------------------
       ! Set the number of total (available) cross-track positions
       ! ---------------------------------------------------------
       nxtloc = nxtrack

       ! ----------------
       ! Shift BrO arrays
       ! ----------------
       IF ( ctrvar%yn_prefit(1) ) THEN
          prefit_col (1:nxtloc,iloop) = prefit_col (1:nxtloc,iloop)
          prefit_dcol(1:nxtloc,iloop) = prefit_dcol(1:nxtloc,iloop)
       END IF

    END DO

    RETURN
  END SUBROUTINE read_prefit_columns


  SUBROUTINE he5_read_prefit_columns (                       &
       swath_id, ntimes_mol, nxtrack_mol, iline,             &
       mcol_len,  molcol_field, col_mol, mdcol_len,          &
       moldcol_field, dcol_mol, yn_read_amf, errstat )

    USE OMSAO_omidata_module
    USE OMSAO_he5_module
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_read_prefit_columns'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4),       INTENT (IN) :: iline, swath_id, ntimes_mol, nxtrack_mol
    INTEGER   (KIND=i4),       INTENT (IN) :: mcol_len, mdcol_len
    CHARACTER (LEN=mcol_len) , INTENT (IN) :: molcol_field
    CHARACTER (LEN=mdcol_len), INTENT (IN) :: moldcol_field
    LOGICAL,                   INTENT (IN) :: yn_read_amf

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4),                                          INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1), INTENT (OUT)   :: col_mol, dcol_mol

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    REAL    (KIND=r8), DIMENSION (nxtrack_mol, 0:ntimes_mol-1) :: amf


    locerrstat = pge_errstat_ok

    ! ----------------------------------------------------
    ! Read current data block fitting output from HE5 file
    ! ----------------------------------------------------
    he5_start_2d = (/ 0, iline /) ; he5_stride_2d = (/ 1, 1 /) 
    he5_edge_2d  = (/ nxtrack_mol, ntimes_mol /)

    ! -----------------------------
    ! Column amount and uncertainty
    ! -----------------------------
    locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(molcol_field)),     &
         he5_start_2d, he5_stride_2d, he5_edge_2d, col_mol(1:nxtrack_mol,0:ntimes_mol-1) )
    locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(moldcol_field)),    &
         he5_start_2d, he5_stride_2d, he5_edge_2d, dcol_mol(1:nxtrack_mol,0:ntimes_mol-1) )


    ! --------------------------------------
    ! Air Mass Factor, but only if requested
    ! (e.g., for BrO, but not O3)
    ! --------------------------------------
    IF ( yn_read_amf ) THEN
       locerrstat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(amfmol_field)),     &
            he5_start_2d, he5_stride_2d, he5_edge_2d, amf(1:nxtrack_mol,0:ntimes_mol-1) )
       WHERE ( amf > 0.0_r8 )
          col_mol  =  col_mol * amf
          dcol_mol = dcol_mol * amf
       END WHERE
    END IF

    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, &
         OMSAO_E_HE5SWRDFLD, modulename, vb_lev_default, errstat )

    RETURN
  END SUBROUTINE he5_read_prefit_columns

  FUNCTION he5_close_prefit_file ( swath_id, swath_file_id ) RESULT ( he5stat )

    !------------------------------------------------------------------------------
    ! This function detatches from the HE5 swath and closes the HE5 input file.
    !
    ! Input: 
    !
    !   swath_id       - ID of pre-fitted input swath      (for BrO or O3)
    !   swath_file_id  - ID of pre-fitted input swath file (for BrO or O3)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_he5_module, ONLY: HE5_SWdetach, HE5_SWclose
    USE OMSAO_errstat_module, ONLY: pge_errstat_ok, he5_stat_ok, &
         pge_errstat_warning, omsao_w_he5swclose, vb_lev_default, &
         error_check

    IMPLICIT NONE

    ! ---------------------------------------
    ! Name of this module/subroutine/function
    ! ---------------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_close_prefit_file'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id, swath_file_id

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4) :: locerr


    he5stat = pge_errstat_ok
    locerr  = pge_errstat_ok

    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    locerr = HE5_SWdetach ( swath_id )
    locerr = HE5_SWclose  ( swath_file_id )
    CALL error_check ( locerr, HE5_STAT_OK, pge_errstat_warning, &
         OMSAO_W_HE5SWCLOSE, modulename, vb_lev_default, he5stat )

    RETURN
  END FUNCTION he5_close_prefit_file


END MODULE OMSAO_prefitcol_module
