SUBROUTINE omps_data_to_omi_variables (OMPS_data,nt,nx,nw)

  ! -----------------
  ! gga Febraury 2014
  ! -----------------
  ! ------------------------------------------------
  ! Modules containing the variables to be filled up
  ! ------------------------------------------------
  USE OMSAO_omidata_module
  USE OMSAO_OMPS_READER
  USE OMSAO_he5_module, ONLY: granule_month
  USE OMSAO_variables_module, ONLY: ctrvar

  IMPLICIT NONE
  
  ! ---------------
  ! Input variables
  ! ---------------
  TYPE (OMPS_NMEV_type), INTENT(IN) :: OMPS_data
  INTEGER (KIND=i4), INTENT(IN) :: nt, nx, nw

  ! ---------------
  ! Local varaibles
  ! ---------------
  INTEGER (KIND=i4)                      :: ix, it, j, imin, imax, icnt
  REAL    (KIND=r8), DIMENSION(nw,nx)    :: tmp_wvl_irra
  CHARACTER(LEN=50)                      :: dummy

  ! ------------------------------------------------
  ! I'm going to try to fill up as many variables as
  ! possible but the OMPS files we have are not as
  ! complete as the OMI l1b files and I don't have
  ! good documentation on the OMPS files
  ! ------------------------------------------------
  omi_auraalt(0:nt-1) = REAL(OMPS_data%SpacecraftAltitude(1:nt),KIND=4)
  omi_xtrflg(1:nx,0:nt-1) = OMPS_data%GroundPixelQualityFlags(1:nx,1:nt)
  omi_instrument_flag(0:nt-1) = OMPS_data%InstrumentQualityFlags(1:nt)

  omi_latitude(1:nx,0:nt-1)        = OMPS_data%Latitude(1:nx,1:nt)
  omi_longitude(1:nx,0:nt-1)       = OMPS_data%Longitude(1:nx,1:nt)
  omi_szenith(1:nx,0:nt-1)         = OMPS_data%SolarZenithAngle(1:nx,1:nt)
  omi_sazimuth(1:nx,0:nt-1)        = OMPS_data%SolarAzimuth(1:nx,1:nt)
  omi_vzenith(1:nx,0:nt-1)         = OMPS_data%SatelliteZenithAngle(1:nx,1:nt)
  omi_vazimuth(1:nx,0:nt-1)        = OMPS_data%SatelliteAzimuth(1:nx,1:nt)
  
  dummy = TRIM(OMPS_DATA%UTC_CCSDS_A(1))
  read(dummy(6:7),'(i2)') granule_month

  ! --------------------
  ! Initialize variables
  ! --------------------
  nwavel         = OMPS_data%nWavel
  omi_nwav_irrad = OMPS_data%nWavel
  omi_nwav_rad   = OMPS_data%nWavel

  ! ------------------------------------------------------------
  ! Fill up pixel selection variables using the set up of the
  ! control file (omi_ccdpix_selection and omi_ccdpix_exclusion)
  ! Limit the pixels used to the fitting window.
  ! Radiance and irradiance
  ! ------------------------------------------------------------
  ! First the irradiance, only loop on xtrack
  ! -----------------------------------------
  tmp_wvl_irra(1:nw,1:nx) = OMPS_data%SolarFluxWavelengths(1:nw,1:nx)
  DO ix = 1, nx
     DO j = 1, 3, 2
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winwav_lim(j  ),KIND=r8), 'LE', &
             omi_ccdpix_selection(ix,j  ) )
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winwav_lim(j+1),KIND=r8), 'GE', &
             omi_ccdpix_selection(ix,j+1) )
     END DO

     imin = omi_ccdpix_selection(ix,1)
     imax = omi_ccdpix_selection(ix,4)

     icnt = imax - imin + 1
     omi_irradiance_wavl(1:icnt,ix) = REAL(OMPS_data%SolarFluxWavelengths(imin:imax,ix), KIND=r8)
     omi_irradiance_spec(1:icnt,ix) = REAL(OMPS_data%SolarFlux(imin:imax,ix), KIND=r8)
     omi_nwav_irrad (ix) = icnt
     omi_sol_wav_avg(ix) = SUM( OMPS_data%SolarFluxWavelengths(imin:imax,ix) ) / REAL(icnt, KIND=r8)
     

     omi_ccdpix_exclusion(ix,1:2) = -1
     IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
             nw, omi_irradiance_wavl(1:nw,ix), REAL(ctrvar%fit_winexc_lim(1),KIND=r8), 'GE', &
             omi_ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nw, omi_irradiance_wavl(1:nw,ix), REAL(ctrvar%fit_winexc_lim(2),KIND=r8), 'LE', &
             omi_ccdpix_exclusion(ix,2) )
     END IF
  END DO

  ! ------------------------------------------
  ! Now the radiances, xtrack and ntimes loops
  ! We don't need to work out again the 
  ! omi_ccdpix_selection arrays
  ! ------------------------------------------
  DO it = 0, nt-1
     DO ix = 1, nx

        imin = omi_ccdpix_selection(ix,1)
        imax = omi_ccdpix_selection(ix,4)
        icnt = imax - imin + 1
        omi_radiance_wavl(1:icnt,ix,it) = REAL ( OMPS_data%BandCenterWavelengths(imin:imax,ix,it+1), KIND=r8 )
        omi_radiance_spec(1:icnt,ix,it) = REAL ( OMPS_data%Radiance(imin:imax,ix,it+1), KIND=r8 )
        omi_radiance_prec(1:icnt,ix,it) = REAL ( OMPS_data%RadianceError(imin:imax,ix,it+1), KIND=r8 )
        omi_radiance_qflg(1:icnt,ix,it) = OMPS_data%PixelQualityFlags(imin:imax,ix,it+1)
        omi_nwav_rad     (       ix,it) = icnt

     END DO
  END DO

END SUBROUTINE omps_data_to_omi_variables
