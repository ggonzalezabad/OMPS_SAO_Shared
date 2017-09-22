MODULE OMSAO_Reference_sector_module

  ! ------------------------------------------------------------------
  ! This module defines variables associated with the Reference Sector
  ! Correction and contains the subroutines needed to apply it to the
  ! HCHO retrieval.
  ! ------------------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE OMSAO_parameters_module,   ONLY: maxchlen, r8_missval, r4_missval, &
       i4_missval, i2_missval, i1_missval
  USE OMSAO_variables_module, ONLY: pcfvar
  USE OMSAO_omidata_module,      ONLY: omi_column_amount, &
       omi_column_uncert, omi_latitude, omi_longitude
  USE OMSAO_he5_module,          ONLY: granule_month, lat_field

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! maxngrid: Parameter to define maximum size of arrays; set to 1000
  ! grid_lat: Latitudes of GEOS-Chem background levels from the Pacif
  !            ic, Hawaii is not included
  ! background_correction: Background level-reference_sector_concentra
  !                        tion
  ! Reference_sector_concentration: Total column (molecules/cm2) obtai
  !                                 ned using the GEOS-Chem climatolog
  !                                 y by D. Millet
  ! ngridpoints: Number of points in grid_lat and Reference_Sector_con
  !              centration
  ! ------------------------------------------------------------------

  INTEGER (KIND=i2), PARAMETER               :: maxngrid = 1000
  REAL    (KIND=r8), DIMENSION (maxngrid)    :: grid_lat, background_correction
  REAL    (KIND=r8), DIMENSION (maxngrid,12) :: Reference_sector_concentration
  REAL    (KIND=r8), DIMENSION (maxngrid)    :: Ref_column_month
  INTEGER (KIND=i2)                          :: ngridpoints

  CONTAINS
    
    SUBROUTINE Reference_Sector_correction (ntimes, nxtrack, ntimesrr, nxtrackrr, refmqf, refamf, &
         amfflg, errstat)

      ! --------------------------------------------------------------------
      ! This subroutine is a wrapper for the Reference Background correction
      ! --------------------------------------------------------------------
      ! VARIABLE DECLARATION
      ! --------------------
      IMPLICIT NONE
      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack, ntimesrr, nxtrackrr
      INTEGER (KIND=i2), DIMENSION (1:nxtrackrr,0:ntimesrr-1), INTENT (IN) :: refmqf
      REAL    (KIND=r8), DIMENSION (1:nxtrackrr,0:ntimesrr-1), INTENT (IN) :: refamf
      INTEGER (KIND=i2), DIMENSION (1:nxtrackrr,0:ntimesrr-1), INTENT (IN) :: amfflg

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4), INTENT (INOUT) :: errstat

      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4) :: itrack, iline
      REAL    (KIND=r8), DIMENSION (1) :: one_background_correction, one_lat
      REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1) :: saocol, saodco, saorms
      REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1) :: saolat, saoamf
      INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1) :: saoqf, amfdiag, rscdiag


      ! ------------------------
      ! Error handling variables
      ! ------------------------
      INTEGER (KIND=i4) :: locerrstat
      
      locerrstat = pge_errstat_ok

      rscdiag = i2_missval
      ! -------------------------------------------------
      ! Read the concentrations from the Reference Sector
      ! Total column molecules/cm2
      ! -------------------------------------------------
      CALL Read_reference_sector_concentration(errstat)
      
      ! ------------------------------------------------
      ! Now we work out the median of the mem_correction
      ! for each grid bin and store it in:
      ! background_correction
      ! ------------------------------------------------
      CALL compute_background_median_and_correction(ntimesrr, nxtrackrr,               &
           omi_column_amount(1:nxtrackrr,0:ntimesrr-1),                                &
           refamf(1:nxtrackrr,0:ntimesrr-1), omi_latitude(1:nxtrackrr,0:ntimesrr-1),   &
           omi_longitude(1:nxtrackrr,0:ntimesrr-1), refmqf(1:nxtrackrr,0:ntimesrr-1),  &
           amfflg(1:nxtrackrr,0:ntimesrr-1), locerrstat)
      ! --------------------------------------------------------------------
      ! Now I have to recover the vcd columns and the amf for the granule of
      ! interest from the L2 file, also the latitudes
      ! --------------------------------------------------------------------
      CALL saopge_columninfo_read (ntimes, nxtrack, saocol, saodco, saorms, &
           saoamf, saoqf, amfdiag, locerrstat)
      CALL saopge_geofield_read(ntimes, nxtrack, lat_field, saolat, locerrstat)

      ! -----------------------------------------------------------------------------
      ! Apply the correction to the slant columns and covert back to vertical columns
      ! -----------------------------------------------------------------------------

      DO itrack = 1, nxtrack
         DO iline = 0, ntimes-1

            IF (saoqf(itrack,iline) .GE. 0 .AND. saoqf(itrack,iline) .LE. 1 .AND. &
                amfdiag(itrack,iline) .EQ. 0) THEN

               one_lat(1) = REAL(saolat(itrack,iline), KIND=8)

               CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),         &
                    grid_lat(1:ngridpoints), background_correction(1:ngridpoints), &
                    1, one_lat(1), one_background_correction(1), locerrstat )

               IF ( ABS(one_background_correction(1)) .LT. 1e18 .AND. &
                    ABS(one_background_correction(1)) .GT. 1e13 ) THEN

                  saocol(itrack,iline) = saocol(itrack,iline) + one_background_correction(1) / saoamf(itrack,iline)
                  rscdiag(itrack,iline) = 0
               ENDIF

            END IF

         END DO
      END DO

      ! ------------------------------------------------------
      ! Now I can write the reference sector correction fields
      ! saocol and rscdiag to the L2 file
      ! ------------------------------------------------------
      CALL he5_write_reference_sector_corrected_column( &
           ntimes, nxtrack, saocol, rscdiag, errstat)

      errstat = MAX ( errstat, locerrstat )

    END SUBROUTINE Reference_Sector_correction

    SUBROUTINE Read_reference_sector_concentration(errstat)

      USE OMSAO_indices_module,   ONLY: OMSAO_refseccor_lun
      
      IMPLICIT NONE
      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4),                                  INTENT (INOUT) :: errstat
      ! ---------------------------------
      ! External OMI and Toolkit routines
      ! ---------------------------------
      INTEGER (KIND=i4), EXTERNAL :: &
           pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef
      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)           :: funit, igrid
      CHARACTER(LEN=1), PARAMETER :: hstr='#'
      LOGICAL                     :: file_header
      CHARACTER (LEN=maxchlen)    :: header_line
      ! ------------------------
      ! Error handling variables
      ! ------------------------
      INTEGER (KIND=i4) :: version, locerrstat, ios
      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=35), PARAMETER :: modulename = 'Read_reference_sector_concentration'
      
      locerrstat = pge_errstat_ok

      ! --------------------
      ! Initialize variables
      ! --------------------
      ngridpoints = 0.0_r8
      grid_lat(1:maxngrid) =  0.0_r8
      Reference_sector_concentration(1:maxngrid,1:12) = 0.0_r8

      ! -----------------------------------------
      ! Open Reference Sector concentrations file
      ! -----------------------------------------
      version = 1
      locerrstat = PGS_IO_GEN_OPENF ( OMSAO_refseccor_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
      locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
      CALL error_check ( &
           locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_REFSECCOR_FILE, &
           modulename//f_sep//TRIM(ADJUSTL(pcfvar%refsec_filename)), vb_lev_default, errstat )
      IF (  errstat /= pge_errstat_ok ) RETURN

      ! ------------
      ! Reading file
      ! -----------------------------------------------------------------------------------
      ! Skip header lines. The header is done when the first character of the line is not #
      ! -----------------------------------------------------------------------------------
      file_header = .TRUE.
      skip_header: DO WHILE (file_header .EQV. .TRUE.)
         READ (UNIT=funit, FMT='(A)', IOSTAT=ios) header_line
         IF ( ios /= 0 ) THEN
            CALL error_check ( &
                 ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
                 modulename//f_sep//TRIM(ADJUSTL(pcfvar%refsec_filename)), vb_lev_default, errstat )
            IF (  errstat /= pge_errstat_ok ) RETURN
         END IF
         IF (header_line(1:1) /= hstr) THEN
            file_header = .FALSE.
         ENDIF
      END DO skip_header

      ! -----------------------------------------------
      ! Read number of grid points. Variable via module
      ! -----------------------------------------------
      READ (UNIT=funit, FMT='(I5)', IOSTAT=ios) ngridpoints
      IF ( ios /= 0 ) THEN
         CALL error_check ( &
              ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
              modulename//f_sep//TRIM(ADJUSTL(pcfvar%refsec_filename)), vb_lev_default, errstat )
         IF (  errstat /= pge_errstat_ok ) RETURN
      END IF

      ! ---------------------------------------------------------
      ! Read reference sector concentrations. Variable via module
      ! ---------------------------------------------------------
      DO igrid = 1, ngridpoints
         READ (UNIT=funit, FMT='(F6.2,2x,12(1x,E14.7))', IOSTAT=ios) grid_lat(igrid), &
              Reference_sector_concentration(igrid, 1:12)
         IF ( ios /= 0 ) THEN
            CALL error_check ( &
                 ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
                 modulename//f_sep//TRIM(ADJUSTL(pcfvar%refsec_filename)), vb_lev_default, errstat )
            IF (  errstat /= pge_errstat_ok ) RETURN
         END IF
      END DO

      ! -----------------------------------------------
      ! Close monthly average file, report SUCCESS read
      ! -----------------------------------------------
      locerrstat = PGS_IO_GEN_CLOSEF ( funit )
      locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
      CALL error_check ( &
           locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_REFSECCOR_FILE, &
           modulename//f_sep//TRIM(ADJUSTL(pcfvar%refsec_filename)), vb_lev_default, errstat )
      IF ( errstat >= pge_errstat_error ) RETURN

      Ref_column_month(1:maxngrid) = Reference_sector_concentration(1:maxngrid,granule_month)
     
    END SUBROUTINE Read_reference_sector_concentration

    SUBROUTINE compute_background_median_and_correction(nTimesRadRR, nXtrackRadRR, &
           mem_column_amount, mem_amf, mem_latitude, mem_longitude, mem_mqf, amfflg, &
           locerrstat)

      USE OMSAO_median_module, ONLY: median

      IMPLICIT NONE
      
      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4), INTENT(IN) :: nTimesRadRR, nXtrackRadRR
      REAL (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: mem_column_amount, mem_amf
      REAL (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: mem_latitude, mem_longitude
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: mem_mqf, amfflg

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4), INTENT(INOUT) :: locerrstat
      
      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i2)                           :: igrid, itrack, iline
      REAL    (KIND=r8), DIMENSION(ngridpoints+1) :: grid
      REAL    (KIND=r8), DIMENSION(1000)          :: into_median
      REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
                         :: column_amounts, amfs
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
                         :: yn_median
      INTEGER (KIND=i4)                           :: npixels, ipixel

      ! -------------------
      ! Routine starts here
      ! -------------------
      locerrstat = pge_errstat_ok
      background_correction = 0.0_r8

      ! ----------------------------------------------------------------
      ! Computing the median of all the pixels of the radiance reference
      ! granule within two consequtive latitudes of the background conce
      ! trations grid readed from the background_concentrations file.
      ! Only pixels between 180W and 140W are considered. Indeed because
      ! of the conditions impossed on the radiance reference granule, as
      ! close as possible to the 165W, they should be no pixels above
      ! 180W and below 140W.
      ! ----------------------------------------------------------------
      ! Generating grid
      ! ---------------
      DO igrid = 1, ngridpoints+1
         grid(igrid) = -90.0_r8 + (180.0_r8 / (REAL(ngridpoints, KIND=r8)) * &
                       (REAL(igrid-1, KIND=r8)))
      END DO

      DO igrid = 1, ngridpoints

         column_amounts        = 0.0_r8
         yn_median             = 0_i2
         into_median           = 0.0_r8

         ! -------------------------------------------------------------------
         ! Finging pixels in the granule between grid(igrid) and grid(igrid+1)
         ! Mid point is grid_lat(igrid).
         ! -------------------------------------------------------------------
         WHERE (mem_latitude .gt. grid(igrid) .AND. mem_latitude .lt. grid(igrid+1) &
                .AND. mem_mqf .EQ. 0 .AND. mem_longitude .GT. -180 .AND.            &
                mem_longitude .LT. -140 .AND. amfflg .EQ. 0 )
            column_amounts = mem_column_amount
            amfs           = mem_amf
            yn_median      = 1
         END WHERE

         npixels = 0_i4
         npixels = SUM(yn_median)
         ipixel  = 1_i4
         DO itrack = 1, nXtrackRadRR
            DO iline = 0, nTimesRadRR-1
               IF (yn_median(itrack,iline) .EQ. 1) THEN
                  into_median(ipixel) = (ref_column_month(igrid) - column_amounts(itrack,iline)) * amfs(itrack,iline)
                  ipixel = ipixel + 1
               END IF
            END DO
         END DO
         IF (npixels .GT. 0) THEN
            background_correction(igrid) = median(npixels, into_median(1:npixels))
         END IF

      END DO

    END SUBROUTINE COMPUTE_BACKGROUND_median_and_correction

    SUBROUTINE he5_write_reference_sector_corrected_column( &
         nt, nx, column, rscdiag, errstat)

      USE OMSAO_he5_module
      USE OMSAO_omidata_module, ONLY: n_roff_dig

      IMPLICIT NONE

      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
      REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: column
      INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: rscdiag

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4), INTENT (INOUT) :: errstat

      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=43), PARAMETER :: modulename = &
           'he5_write_reference_sector_corrected_column'

      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)                          :: locerrstat
      REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc
  
      ! -----------------------------------------------
      ! Total column corrected by the reference sector.
      ! -----------------------------------------------
      ! Only implemented for HCHO.
      ! -----------------------------------------------
      
      locerrstat = pge_errstat_ok

      he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)      

      colloc = column

      CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
      locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rscol_field)), &
           he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
      errstat = MAX ( errstat, locerrstat )
  
      locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rsfla_field)), &
           he5_start_2d, he5_stride_2d, he5_edge_2d, rscdiag(1:nx,0:nt-1) )
      errstat = MAX ( errstat, locerrstat )

      ! ------------------
      ! Check error status
      ! ------------------
      CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
           modulename, vb_lev_default, errstat )

      RETURN
  
    END SUBROUTINE he5_write_reference_sector_corrected_column

END MODULE OMSAO_Reference_sector_module
