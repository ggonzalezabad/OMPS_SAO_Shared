SUBROUTINE OMSAO_main ( exit_value )

  ! ********************************************************************
  !
  ! This is the main program of the SAO Ozone Mapping and Profiler Suite
  ! Nadir Mapper (OMPS NM) Trace Gas retrievals.
  !
  ! Authors: Gonzalo Gonzalez Abad
  !          Smithsonian Astrophysical Observatory
  !          60 Garden Street (MS 50)
  !          Cambridge, MA 02138 (USA)
  !
  !          EMail: ggonzalezabad@cfa.harvard.edu
  !
  ! ********************************************************************

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_variables_module, ONLY: pcfvar
  USE OMSAO_pcf_file_module, ONLY: read_pcf_file
  USE OMSAO_reference_spectra, ONLY: read_reference_spectra
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, &
       pge_errstat_error, pge_errstat_fatal, omsao_a_subroutine, &
       omsao_f_subroutine, f_sep, error_check, &
       pge_error_status_exit, vb_lev_screen

  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: exit_value


  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=10), PARAMETER :: modulename = 'OMSAO_main'

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late 
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER   (KIND=i4)      :: errstat, pge_error_status

  ! --------------------------------------------------------------------
  ! Maximum number of points in any reference spectrum. This is used for
  ! automatic (i.e., SUBROUTINE argument) memory allocation and hence is
  ! defined here instead of in a MODULE.
  ! --------------------------------------------------------------------
  INTEGER (KIND=i4) :: n_max_rspec

  ! ----------------------------
  ! Set PGE_ERROR_STATUS to O.K.
  ! ----------------------------
  pge_error_status = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  CALL unbufferSTDout()                       ! Make PGE write STD/IO unbuffered
  ! ----------------------------------------------------------------------------
  errstat = pge_errstat_ok

  ! ---------------------------------------------------------------------------
  CALL read_pcf_file ( errstat )   ! Read PCF file
  ! ---------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"READ_PCF_FILE.", vb_lev_screen, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666
  
  ! ---------------------------------------------------
!  CALL init_metadata ( errstat )  ! Initialize MetaData
  ! ---------------------------------------------------
!  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_F_SUBROUTINE, &
!       modulename//f_sep//"INIT_METADATA.", vb_lev_screen, pge_error_status )
!  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! -------------------------------------------------------------------------------
  CALL read_reference_spectra ( n_max_rspec, errstat )     ! Read reference spectra
  ! -------------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"READ_REFERENCE_SPECTRA.", vb_lev_screen, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! ---------------------------------------------
  ! Set number of InputPointers and InputVersions
  ! ---------------------------------------------
  CALL set_input_pointer_and_versions ( errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_SUBROUTINE, &
       modulename//f_sep//"SET_INPUT_POINTER_AND_VERSIONS.", vb_lev_screen, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! --------------------------------------------------------------------------------------------
  CALL omi_pge_fitting_process  ( pcfvar%pge_idx, n_max_rspec, errstat )   ! Where all the work is done
  ! --------------------------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_warning, errstat, OMSAO_A_SUBROUTINE, &
       modulename//f_sep//"OMI_PGE_FITTING_PROCESS.", vb_lev_screen, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
666 CALL pge_error_status_exit ( pge_error_status, exit_value )

  RETURN
END SUBROUTINE OMSAO_main
