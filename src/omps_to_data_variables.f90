SUBROUTINE omps_to_data_variables (omps_data,nt,nx,nw)

  ! -----------------
  ! gga Febraury 2014
  ! -----------------
  ! ------------------------------------------------
  ! Modules containing the variables to be filled up
  ! ------------------------------------------------
  USE OMSAO_precision_module, ONLY: i4, r8
  USE OMSAO_data_module, ONLY: spacecraft_alt, time, xtrflg, instrument_flag, &
       latitude, longitude, szenith, sazimuth, vzenith, vazimuth, vzenith, &
       vazimuth, razimuth, nwav_irrad, nwav_rad, ccdpix_selection, &
       radiance_wavl, radiance_spec, radiance_prec, radiance_qflg, nwav_rad, &
       omi_sol_wav_avg, irradiance_wavl, irradiance_spec, ccdpix_exclusion, &
       EarthSunDistance, yn_process_pixel, ntime_rad, nxtrack_rad
  USE OMSAO_OMPS_READER, ONLY: omps_nmev_type
  USE OMSAO_he5_module, ONLY: granule_day, granule_month, granule_year
  USE OMSAO_variables_module, ONLY: ctrvar

  IMPLICIT NONE
  
  ! ---------------
  ! Input variables
  ! ---------------
  TYPE (OMPS_NMEV_type), INTENT(IN) :: omps_data
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
  ntime_rad = nt
  nxtrack_rad = nx
  spacecraft_alt(0:nt-1) = REAL(omps_data%SpacecraftAltitude(1:nt),KIND=4)
  time(0:nt-1) = REAL(omps_data%ImageMidPoint_TAI93(1:nt),KIND=4)
  xtrflg(1:nx,0:nt-1) = OMPS_data%GroundPixelQualityFlags(1:nx,1:nt)
  instrument_flag(0:nt-1) = OMPS_data%InstrumentQualityFlags(1:nt)

  yn_process_pixel = .FALSE.
  DO it = 0, nt-1
     DO ix = 1, nx
        IF (OMPS_data%ReportIntQualFlags(it+1) == 0 .AND. &
            xtrflg(ix,it) == 0) yn_process_pixel = .TRUE.
     END DO
  END DO

  latitude(1:nx,0:nt-1) = OMPS_data%Latitude(1:nx,1:nt)
  longitude(1:nx,0:nt-1) = OMPS_data%Longitude(1:nx,1:nt)
  szenith(1:nx,0:nt-1) = OMPS_data%SolarZenithAngle(1:nx,1:nt)
  sazimuth(1:nx,0:nt-1) = OMPS_data%SolarAzimuth(1:nx,1:nt)
  vzenith(1:nx,0:nt-1) = OMPS_data%SatelliteZenithAngle(1:nx,1:nt)
  vazimuth(1:nx,0:nt-1) = OMPS_data%SatelliteAzimuth(1:nx,1:nt)
  razimuth(1:nx,0:nt-1) = (sazimuth(1:nx,0:nt-1) + 180.0 - vazimuth(1:nx,0:nt-1))

  earthsundistance = REAL(OMPS_data%SunEarthDistance,KIND=4)

  dummy = TRIM(OMPS_DATA%UTC_CCSDS_A(FLOOR(REAL(nt)/2.0)))
  read(dummy(9:10),'(i2)') granule_day
  read(dummy(6:7),'(i2)') granule_month
  read(dummy(1:4),'(i4)') granule_year

  ! --------------------
  ! Initialize variables
  ! --------------------
  nwav_irrad = OMPS_data%nWavel
  nwav_rad = OMPS_data%nWavel

  ! ------------------------------------------------------------
  ! Fill up pixel selection variables using the set up of the
  ! control file (ccdpix_selection and ccdpix_exclusion)
  ! Limit the pixels used to the fitting window.
  ! Radiance and irradiance
  ! ------------------------------------------------------------
  ! First the irradiance, only loop on xtrack
  ! -----------------------------------------
  tmp_wvl_irra(1:nw,1:nx) = REAL(OMPS_data%SolarFluxWavelengths(1:nw,1:nx),KIND=r8)
  DO ix = 1, nx
     DO j = 1, 3, 2
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winwav_lim(j  ),KIND=r8), 'LE', &
             ccdpix_selection(ix,j  ) )
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winwav_lim(j+1),KIND=r8), 'GE', &
             ccdpix_selection(ix,j+1) )
     END DO

     imin = ccdpix_selection(ix,1)
     imax = ccdpix_selection(ix,4)

     icnt = imax - imin + 1
     irradiance_wavl(1:icnt,ix) = REAL(OMPS_data%SolarFluxWavelengths(imin:imax,ix), KIND=r8)
     irradiance_spec(1:icnt,ix) = REAL(OMPS_data%SolarFlux(imin:imax,ix), KIND=r8)
     nwav_irrad (ix) = icnt
     omi_sol_wav_avg(ix) = SUM(irradiance_wavl(1:icnt,ix) ) / REAL(icnt, KIND=r8)
     
     ccdpix_exclusion(ix,1:2) = -1
     IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winexc_lim(1),KIND=r8), 'GE', &
             ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nw, tmp_wvl_irra(1:nw,ix), REAL(ctrvar%fit_winexc_lim(2),KIND=r8), 'LE', &
             ccdpix_exclusion(ix,2) )
     END IF
  END DO

  ! ------------------------------------------
  ! Now the radiances, xtrack and ntimes loops
  ! We don't need to work out again the 
  ! ccdpix_selection arrays
  ! ------------------------------------------
  DO it = 0, nt-1
     DO ix = 1, nx
        imin = ccdpix_selection(ix,1)
        imax = ccdpix_selection(ix,4)
        icnt = imax - imin + 1
        radiance_wavl(1:icnt,ix,it) = REAL ( OMPS_data%BandCenterWavelengths(imin:imax,ix,it+1), KIND=r8 )
        radiance_spec(1:icnt,ix,it) = REAL ( OMPS_data%Radiance(imin:imax,ix,it+1), KIND=r8 )
        radiance_prec(1:icnt,ix,it) = REAL ( OMPS_data%RadianceError(imin:imax,ix,it+1), KIND=r8 )
        radiance_qflg(1:icnt,ix,it) = OMPS_data%PixelQualityFlags(imin:imax,ix,it+1)
        nwav_rad (ix,it) = icnt
     END DO
  END DO

END SUBROUTINE omps_to_data_variables
