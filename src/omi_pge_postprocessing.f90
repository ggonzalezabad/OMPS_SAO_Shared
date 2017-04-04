SUBROUTINE omi_pge_postprocessing ( &
     pge_idx, ntimes, nxtrack, xtrange, errstat )

  ! ---------------------------------------------------------
  ! In this subroutine we collect all those computations that
  ! are done "post fitting". These include:
  !
  ! (1) AMF calculation
  ! (2) Fitting statistics
  ! (3) Cross-track destriping
  ! (4) Ground-pixel corner computation
  ! (5) Reference Sector Background Correction for HCHO
  ! ---------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8
  USE OMSAO_he5_module, ONLY: lat_field, lon_field, sza_field, &
       thgt_field, vza_field, extr_field
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok
  USE OMSAO_wfamf_module, ONLY: amf_calculation_bis, climatology_allocate, &
       Cmlat, Cmlon, CmETA, CmEp1

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                              INTENT (IN) :: ntimes, nxtrack, pge_idx
  INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),  INTENT (IN) :: xtrange

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Local variables
  ! ----------------
  ! (1) OMI data
  ! ----------------
  REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1) :: lat, lon, sza, vza, thg, extr
  REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: saocol, saodco, saorms, saoamf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: saofcf, saomqf

  ! --------------
  ! Error handling
  ! --------------
  INTEGER (KIND=i4) :: locerrstat

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = pge_errstat_ok

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  CALL  saopge_geofield_read ( ntimes, nxtrack, lat_field,  lat,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, lon_field,  lon,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, sza_field,  sza,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, vza_field,  vza,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, thgt_field, thg,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, extr_field, extr, locerrstat )
 
  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  CALL saopge_columninfo_read (                 &
       ntimes, nxtrack, saocol, saodco, saorms, &
       saoamf, saofcf, locerrstat                 )

  ! ----------------------------------
  ! Compute average fitting statistics
  ! ----------------------------------
  CALL compute_fitting_statistics (                      &
       pge_idx, ntimes, nxtrack, xtrange,                &
       saocol, saodco, saorms, saofcf, saomqf, locerrstat )

  ! --------------------------------
  ! Deallocate Climatology variables
  ! --------------------------------
  CALL climatology_allocate ( "d", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )


  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE omi_pge_postprocessing
