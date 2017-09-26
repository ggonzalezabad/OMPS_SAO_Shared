!! MODULE to READ OMPS data
!!
! Created by gga January 2014
!!
! All he5 internals are done by the belgium (BIRA) module ReadH5dataset
! with small modifications
!!

MODULE OMSAO_omps_reader

  USE ReadH5dataset, ONLY: H5ReadDataset, H5ReadAttribute, ErrorFlag, &
       ErrorMessage ! Contains the generic HDF5 reading routines

  PRIVATE
  PUBLIC :: omps_nmev_type, omps_nmto3_type, omps_nmev_reader, omps_nmto3_reader
  INTEGER (KIND = 4), PARAMETER :: MAX_NAME_LENGTH = 256

  TYPE omps_nmev_type
     CHARACTER (LEN=MAX_NAME_LENGTH) :: filename
     INTEGER(KIND=4)                 :: nLines
     INTEGER(KIND=4)                 :: nXtrack
     INTEGER(KIND=4)                 :: nWavel
     ! ======================
     ! Data storage variables
     ! ======================
     ! Calibration data
     REAL(KIND=4), DIMENSION(:,:,:),   POINTER :: BandCenterWavelengths => NULL()
     INTEGER(KIND=2), DIMENSION(:,:,:),POINTER :: CCDRowColIndicies => NULL()
     REAL(KIND=4), DIMENSION(:,:),     POINTER :: DarkCurrentCorrection => NULL()
     REAL(KIND=4), DIMENSION(:,:),     POINTER :: RadianceCalCoeff => NULL()
     REAL(KIND=4), DIMENSION(:,:,:),   POINTER :: SmearCorrection => NULL()
     REAL(KIND=4), DIMENSION(:,:),     POINTER :: SolarFlux => NULL()
     REAL(KIND=4), DIMENSION(:,:),     POINTER :: SolarFluxWavelengths => NULL()
     REAL(KIND=4), DIMENSION(:,:,:),   POINTER :: StrayLightCorrection => NULL()
     ! Geolocation data
     REAL(KIND=8), DIMENSION(:), POINTER :: ECISolarDeclination => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: ECISolarRightAscension => NULL()
     REAL(KIND=8), DIMENSION(:,:), POINTER :: ECISolarUnitVector => NULL()
     REAL(KIND=8), DIMENSION(:,:), POINTER :: ECISpacecraftPosition => NULL()
     REAL(KIND=8), DIMENSION(:,:), POINTER :: ECISpacecraftVelocity => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: GoniometricSolarAzimuth => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: GoniometricSolarElevation => NULL()
     INTEGER(KIND=2), DIMENSION(:,:), POINTER :: GroundPixelQualityFlags => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: ImageMidpoint_TAI93 => NULL()
     INTEGER(KIND=2), DIMENSION(:), POINTER :: InstrumentQualityFlags => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Latitude => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: LatitudeCorner => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Longitude => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: LongitudeCorner => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SatelliteAzimuth => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SatelliteZenithAngle => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SolarAzimuth => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: SolarBeta => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SolarZenithAngle => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: SpacecraftAltitude => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: SpacecraftLatitude => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: SpacecraftLongitude => NULL()
     REAL(KIND=8), DIMENSION(:), POINTER :: SubSatelliteSolarZenithAngle => NULL()
     REAL(KIND=8), POINTER :: SunEarthDistance => NULL()
     CHARACTER(LEN=MAX_NAME_LENGTH), DIMENSION(:), POINTER :: UTC_CCSDS_A => NULL()
     ! Input pointers
     CHARACTER(LEN=7000), POINTER :: ControlFileContents => NULL() !7000 expecting it to be long enough
     ! Science data
     REAL(KIND=4), DIMENSION(:), POINTER :: ExposureTime => NULL()
     INTEGER(KIND=1), DIMENSION(:), POINTER :: NumberCoadds => NULL()
     INTEGER(KIND=2), DIMENSION(:,:,:), POINTER :: PixelQualityFlags => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: Radiance => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: RadianceError => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: RawCounts => NULL()
     INTEGER(KIND=4), DIMENSION(:), POINTER :: ReportIntQualFlags => NULL()
     INTEGER(KIND=4), DIMENSION(:), POINTER :: SensorStatusBits => NULL()
  END TYPE omps_nmev_type

  TYPE omps_nmto3_type
     CHARACTER (LEN=MAX_NAME_LENGTH) :: filename
     INTEGER(KIND=4) :: nLines
     INTEGER(KIND=4) :: nXtrack
     INTEGER(KIND=4) :: nWavel
     ! ======================
     ! Data storage variables
     ! ======================
     ! Ancillary data
     REAL(KIND=4), DIMENSION(:,:), POINTER :: CloudPressure => NULL()
     INTEGER(KIND=4), DIMENSION(:,:), POINTER :: SurfaceCategory => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: TerrainPressure => NULL()
     ! Geolocation data
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Latitude => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: LatitudeCorner => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Longitude => NULL()
     REAL(KIND=4), DIMENSION(:,:,:), POINTER :: LongitudeCorner => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: RelativeAzimuthAngle => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SolarAzimuthAngle => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: SolarZenithAngle => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: ViewingAzimuthAngle => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: ViewingZenithAngle => NULL()
     ! Science data
     REAL(KIND=4), DIMENSION(:,:), POINTER :: ColumnAmountO3 => NULL()
     INTEGER(KIND=1), DIMENSION(:), POINTER :: MeasurementQualityFlags => NULL()
     INTEGER(KIND=2), DIMENSION(:,:), POINTER :: QualityFlags => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: RadiativeCloudFraction => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Reflectivity331 => NULL()
     REAL(KIND=4), DIMENSION(:,:), POINTER :: Reflectivity360 => NULL()
  END TYPE omps_nmto3_type

CONTAINS

  FUNCTION omps_nmev_reader(this, filename) RESULT(status)

    ! I'm going to use the Belgiums hdf5 reader
    ! The most recent version of the package is available at:
    ! http://uv-vis.aeronomie.be/software/tools/hdf5read.php
    ! All information available on the OMPS-NPP-NPP-TC_SDR_EV_NASA
    ! type files is going to be readed for possible future use.
    ! Stand alone version.
    
    IMPLICIT NONE

    CHARACTER (LEN=*) :: filename
    TYPE (OMPS_NMEV_type), INTENT(INOUT) :: this
    INTEGER(KIND=4) :: status, ierr
    
    ! -------------------------------------
    ! Make sure pointers are not associated
    ! -------------------------------------
    IF (ASSOCIATED(this%BandCenterWavelengths)) DEALLOCATE(this%BandCenterWavelengths)
    IF (ASSOCIATED(this%CCDRowColIndicies)) DEALLOCATE(this%CCDRowColIndicies)
    IF (ASSOCIATED(this%DarkCurrentCorrection)) DEALLOCATE(this%DarkCurrentCorrection)
    IF (ASSOCIATED(this%RadianceCalCoeff)) DEALLOCATE(this%RadianceCalCoeff)
    IF (ASSOCIATED(this%SmearCorrection)) DEALLOCATE(this%SmearCorrection)
    IF (ASSOCIATED(this%SolarFlux)) DEALLOCATE(this%SolarFlux)
    IF (ASSOCIATED(this%SolarFluxWavelengths)) DEALLOCATE(this%SolarFluxWavelengths)
    IF (ASSOCIATED(this%StrayLightCorrection)) DEALLOCATE(this%StrayLightCorrection)
    
    IF (ASSOCIATED(this%ECISolarDeclination)) DEALLOCATE(this%ECISolarDeclination)
    IF (ASSOCIATED(this%ECISolarRightAscension)) DEALLOCATE(this%ECISolarRightAscension)
    IF (ASSOCIATED(this%ECISolarUnitVector)) DEALLOCATE(this%ECISolarUnitVector)
    IF (ASSOCIATED(this%ECISpacecraftPosition)) DEALLOCATE(this%ECISpacecraftPosition)
    IF (ASSOCIATED(this%ECISpacecraftVelocity)) DEALLOCATE(this%ECISpacecraftVelocity)
    IF (ASSOCIATED(this%GoniometricSolarAzimuth)) DEALLOCATE(this%GoniometricSolarAzimuth)
    IF (ASSOCIATED(this%GoniometricSolarElevation)) DEALLOCATE(this%GoniometricSolarElevation)
    IF (ASSOCIATED(this%GroundPixelQualityFlags)) DEALLOCATE(this%GroundPixelQualityFlags)
    IF (ASSOCIATED(this%ImageMidpoint_TAI93)) DEALLOCATE(this%ImageMidpoint_TAI93)
    IF (ASSOCIATED(this%InstrumentQualityFlags)) DEALLOCATE(this%InstrumentQualityFlags)
    IF (ASSOCIATED(this%Latitude)) DEALLOCATE(this%Latitude)
    IF (ASSOCIATED(this%LatitudeCorner)) DEALLOCATE(this%LatitudeCorner)
    IF (ASSOCIATED(this%Longitude)) DEALLOCATE(this%Longitude)
    IF (ASSOCIATED(this%LongitudeCorner)) DEALLOCATE(this%LongitudeCorner)
    IF (ASSOCIATED(this%SatelliteAzimuth)) DEALLOCATE(this%SatelliteAzimuth)
    IF (ASSOCIATED(this%SatelliteZenithAngle)) DEALLOCATE(this%SatelliteZenithAngle)
    IF (ASSOCIATED(this%SolarAzimuth)) DEALLOCATE(this%SolarAzimuth)
    IF (ASSOCIATED(this%SolarBeta)) DEALLOCATE(this%SolarBeta)
    IF (ASSOCIATED(this%SolarZenithAngle)) DEALLOCATE(this%SolarZenithAngle)
    IF (ASSOCIATED(this%SpacecraftAltitude)) DEALLOCATE(this%SpacecraftAltitude)
    IF (ASSOCIATED(this%SpacecraftLatitude)) DEALLOCATE(this%SpacecraftLatitude)
    IF (ASSOCIATED(this%SpacecraftLongitude)) DEALLOCATE(this%SpacecraftLongitude)
    IF (ASSOCIATED(this%SubSatelliteSolarZenithAngle)) DEALLOCATE(this%SubSatelliteSolarZenithAngle)
    IF (ASSOCIATED(this%SunEarthDistance)) DEALLOCATE(this%SunEarthDistance)
    IF (ASSOCIATED(this%UTC_CCSDS_A)) DEALLOCATE(this%UTC_CCSDS_A)

    IF (ASSOCIATED(this%ControlFileContents)) DEALLOCATE(this%ControlFileContents)
    
    IF (ASSOCIATED(this%ExposureTime)) DEALLOCATE(this%ExposureTime)
    IF (ASSOCIATED(this%NumberCoadds)) DEALLOCATE(this%NumberCoadds)
    IF (ASSOCIATED(this%PixelQualityFlags)) DEALLOCATE(this%PixelQualityFlags)
    IF (ASSOCIATED(this%Radiance)) DEALLOCATE(this%Radiance)
    IF (ASSOCIATED(this%RadianceError)) DEALLOCATE(this%RadianceError)
    IF (ASSOCIATED(this%RawCounts)) DEALLOCATE(this%RawCounts)
    IF (ASSOCIATED(this%ReportIntQualFlags)) DEALLOCATE(this%ReportIntQualFlags)
    IF (ASSOCIATED(this%SensorStatusBits)) DEALLOCATE(this%SensorStatusBits)
    
    ! ----------------------------------------------------------
    ! Filling up the dimension variables nLines, nXtrack, nWavel
    ! ----------------------------------------------------------
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/Radiance", this%nLines, this%nXtrack, this%nWavel)
    
    ALLOCATE(this%BandCenterWavelengths(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Band Center Wavelengths'
       GOTO 990
    END IF
    ALLOCATE(this%CCDRowColIndicies(this%nWavel,this%nXtrack,2), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating CCD Row Col Indices'
       GOTO 990
    END IF
    ALLOCATE(this%DarkCurrentCorrection(this%nWavel,this%nXtrack), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Dark Current Correction'
       GOTO 990
    END IF
    ALLOCATE(this%RadianceCalCoeff(this%nWavel,this%nXtrack), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Radiance Cal Coeff'
       GOTO 990
    END IF
    ALLOCATE(this%SmearCorrection(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Smear Correction'
       GOTO 990
    END IF
    ALLOCATE(this%SolarFlux(this%nWavel,this%nXtrack), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Flux'
       GOTO 990
    END IF
    ALLOCATE(this%SolarFluxWavelengths(this%nWavel,this%nXtrack), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Flux Wavelengths'
       GOTO 990
    END IF
    ALLOCATE(this%StrayLightCorrection(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Stray Light Correction'
       GOTO 990
    END IF

    ALLOCATE(this%ECISolarDeclination(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating ECI Solar Declination'
       GOTO 990
    END IF
    ALLOCATE(this%ECISolarRightAscension(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating ECI Solar Right Ascension'
       GOTO 990
    END IF
    ALLOCATE(this%ECISolarUnitVector(3,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating ECI Unit Vector'
       GOTO 990
    END IF
    ALLOCATE(this%ECISpacecraftPosition(3,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating ECI Spacecraft Position'
       GOTO 990
    END IF
    ALLOCATE(this%ECISpacecraftVelocity(3,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating ECI Spacecraft Velocity'
       GOTO 990
    END IF
    ALLOCATE(this%GoniometricSolarAzimuth(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Goniometric Solar Azimuth'
       GOTO 990
    END IF
    ALLOCATE(this%GoniometricSolarElevation(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Goniometric Solar Elevation'
       GOTO 990
    END IF
    ALLOCATE(this%GroundPixelQualityFlags(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Ground Pixel Quality Flags'
       GOTO 990
    END IF
    ALLOCATE(this%ImageMidpoint_TAI93(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Image Mid point'
       GOTO 990
    END IF
    ALLOCATE(this%InstrumentQualityFlags(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Instrument Quality Flags'
       GOTO 990
    END IF
    ALLOCATE(this%Latitude(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Latitude'
       GOTO 990
    END IF
    ALLOCATE(this%LatitudeCorner(4,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Latitude Corner'
       GOTO 990
    END IF
    ALLOCATE(this%Longitude(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Longitude'
       GOTO 990
    END IF
    ALLOCATE(this%LongitudeCorner(4,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Longitude Corner'
       GOTO 990
    END IF
    ALLOCATE(this%SatelliteAzimuth(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Satellite Azimuth'
       GOTO 990
    END IF
    ALLOCATE(this%SatelliteZenithAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Satellite Zenith Angle'
       GOTO 990
    END IF
    ALLOCATE(this%SolarAzimuth(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Azimuth'
       GOTO 990
    END IF
    ALLOCATE(this%SolarBeta(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Beta'
       GOTO 990
    END IF
    ALLOCATE(this%SolarZenithAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Zenith Angle'
       GOTO 990
    END IF
    ALLOCATE(this%SpacecraftAltitude(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Spacecraft Altitude'
       GOTO 990
    END IF
    ALLOCATE(this%SpacecraftLatitude(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Spacecraft Latitude'
       GOTO 990
    END IF
    ALLOCATE(this%SpacecraftLongitude(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Spacecraft Longitude'
       GOTO 990
    END IF
    ALLOCATE(this%SubSatelliteSolarZenithAngle(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Sub Satellite Solar Zenith Angle'
       GOTO 990
    END IF
    ALLOCATE(this%SunEarthDistance, STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Sun Earth Distance'
       GOTO 990
    END IF
    ALLOCATE(this%UTC_CCSDS_A(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating UTC_CCSDS_A'
       GOTO 990
    END IF

    ALLOCATE(this%ControlFileContents, STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Control File Contents'
       GOTO 990
    END IF


    ALLOCATE(this%ExposureTime(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Exposure Time'
       GOTO 990
    END IF
    ALLOCATE(this%NumberCoadds(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Number Coadds'
       GOTO 990
    END IF
    ALLOCATE(this%PixelQualityFlags(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Pixel Quality Flags'
       GOTO 990
    END IF
    ALLOCATE(this%Radiance(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Radiance'
       GOTO 990
    END IF
    ALLOCATE(this%RadianceError(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Radiance Error'
       GOTO 990
    END IF
    ALLOCATE(this%RawCounts(this%nWavel,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Raw Counts'
       GOTO 990
    END IF
    ALLOCATE(this%ReportIntQualFlags(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Report Int Qual Flags'
       GOTO 990
    END IF
    ALLOCATE(this%SensorStatusBits(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Sensor Status Bits'
       GOTO 990
    END IF
    

    ! ============
    ! Reading data
    ! ============
    ! ----------------
    ! Calibration data
    ! ----------------
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/BandCenterWavelengths", &
         this%BandCenterWavelengths)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/CCDRowColIndicies", this%CCDRowColIndicies)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/DarkCurrentCorrection", &
         this%DarkCurrentCorrection)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/RadianceCalCoeff", this%RadianceCalCoeff)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/SmearCorrection", this%SmearCorrection)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/SolarFlux", this%SolarFlux)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/SolarFluxWavelengths", &
         this%SolarFluxWavelengths)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/CalibrationData/StrayLightCorrection", &
         this%StrayLightCorrection)  
    IF (ErrorFlag.NE.0) GOTO 990

    ! ----------------
    ! Geolocation data
    ! ----------------
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ECISolarDeclination", &
         this%ECISolarDeclination)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ECISolarRightAscension", &
         this%ECISolarRightAscension)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ECISolarUnitVector", &
         this%ECISolarUnitVector)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ECISpacecraftPosition", &
         this%ECISpacecraftPosition)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ECISpacecraftVelocity", &
         this%ECISpacecraftVelocity)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/GoniometricSolarAzimuth", &
         this%GoniometricSolarAzimuth)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/GoniometricSolarElevation", &
         this%GoniometricSolarElevation)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/GroundPixelQualityFlags", &
         this%GroundPixelQualityFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/ImageMidpoint_TAI93", this%ImageMidpoint_TAI93)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/InstrumentQualityFlags", &
         this%InstrumentQualityFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/Latitude", &
         this%Latitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/LatitudeCorner", &
         this%LatitudeCorner)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/Longitude", &
         this%Longitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/LongitudeCorner", &
         this%LongitudeCorner)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SatelliteAzimuthAngle", &
         this%SatelliteAzimuth)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SatelliteZenithAngle", &
         this%SatelliteZenithAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SolarAzimuthAngle", &
         this%SolarAzimuth)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SolarBetaAngle", this%SolarBeta)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SolarZenithAngle", &
         this%SolarZenithAngle) 
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SpacecraftAltitude", this%SpacecraftAltitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SpacecraftLatitude", this%SpacecraftLatitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SpacecraftLongitude", this%SpacecraftLongitude)    
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/SubSatelliteSolarZenithAngle", &
         this%SubSatelliteSolarZenithAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadAttribute(filename, &
         "/BinScheme1/GeolocationData/SunEarthDistance", this%SunEarthDistance)    
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/GeolocationData/UTC_CCSDS_A", this%UTC_CCSDS_A)
    IF (ErrorFlag.NE.0) GOTO 990
    
    ! --------------
    ! Input pointers
    ! --------------
    CALL H5ReadDataset(filename, &
         "/InputPointers/ControlFileContents", this%ControlFileContents)
    IF (ErrorFlag.NE.0) GOTO 990
    
    ! ------------
    ! Science data
    ! ------------
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/ExposureTime", this%ExposureTime)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/NumberCoadds", this%NumberCoadds)    
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/PixelQualityFlags", &
         this%PixelQualityFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/Radiance", this%Radiance)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/RadianceError", &
         this%RadianceError)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/RawCounts", &
         this%RawCounts)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/ReportIntervalQualityFlags", this%ReportIntQualFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/BinScheme1/ScienceData/SensorStatusBits", this%SensorStatusBits)
    IF (ErrorFlag.NE.0) GOTO 990
    
    status = ErrorFlag
    RETURN
990 WRITE(*, '(A)') ErrorMessage
    status = 2
    RETURN
  END FUNCTION omps_nmev_reader

  FUNCTION omps_nmto3_reader(this, filename) RESULT(status)

    IMPLICIT NONE

    CHARACTER (LEN=*) :: filename
    TYPE(omps_nmto3_type), INTENT(INOUT) :: this
    INTEGER(KIND=4) :: status, ierr

    ! -------------------------------------
    ! Make sure pointers are not associated
    ! -------------------------------------
     ! Ancillary data
    IF (ASSOCIATED(this%CloudPressure)) DEALLOCATE(this%CloudPressure)
    IF (ASSOCIATED(this%SurfaceCategory)) DEALLOCATE(this%SurfaceCategory)
    IF (ASSOCIATED(this%TerrainPressure)) DEALLOCATE(this%TerrainPressure)
    ! Geolocation data
    IF (ASSOCIATED(this%Latitude)) DEALLOCATE(this%Latitude)
    IF (ASSOCIATED(this%LatitudeCorner)) DEALLOCATE(this%LatitudeCorner)
    IF (ASSOCIATED(this%Longitude)) DEALLOCATE(this%Longitude)
    IF (ASSOCIATED(this%LongitudeCorner)) DEALLOCATE(this%LongitudeCorner)
    IF (ASSOCIATED(this%RelativeAzimuthAngle)) DEALLOCATE(this%RelativeAzimuthAngle)
    IF (ASSOCIATED(this%SolarAzimuthAngle)) DEALLOCATE(this%SolarAzimuthAngle)
    IF (ASSOCIATED(this%SolarZenithAngle)) DEALLOCATE(this%SolarZenithAngle)
    IF (ASSOCIATED(this%ViewingAzimuthAngle)) DEALLOCATE(this%ViewingAzimuthAngle)
    IF (ASSOCIATED(this%ViewingZenithAngle)) DEALLOCATE(this%ViewingZenithAngle)
     ! Science data
    IF (ASSOCIATED(this%ColumnAmountO3)) DEALLOCATE(this%ColumnAmountO3)
    IF (ASSOCIATED(this%MeasurementQualityFlags)) DEALLOCATE(this%MeasurementQualityFlags)
    IF (ASSOCIATED(this%QualityFlags)) DEALLOCATE(this%QualityFlags)
    IF (ASSOCIATED(this%RadiativeCloudFraction)) DEALLOCATE(this%RadiativeCloudFraction)
    IF (ASSOCIATED(this%Reflectivity331)) DEALLOCATE(this%Reflectivity331)
    IF (ASSOCIATED(this%Reflectivity360)) DEALLOCATE(this%Reflectivity360)

    ! ----------------------------------------------------------
    ! Filling up the dimension variables nLines, nXtrack, nWavel
    ! ----------------------------------------------------------
    CALL H5ReadDataset(filename, &
         "/AncillaryData/CloudPressure", this%nLines, this%nXtrack)

    ALLOCATE(this%CloudPressure(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating cloud pressure'
       GOTO 990
    END IF
    ALLOCATE(this%SurfaceCategory(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating surface category'
       GOTO 990
    END IF
    ALLOCATE(this%TerrainPressure(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating terrain pressure'
       GOTO 990
    END IF
    ALLOCATE(this%Latitude(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Latitude'
       GOTO 990
    END IF
    ALLOCATE(this%LatitudeCorner(4,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Latitude Corner'
       GOTO 990
    END IF
    ALLOCATE(this%Longitude(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Longitude'
       GOTO 990
    END IF
    ALLOCATE(this%LongitudeCorner(4,this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Longitude Corner'
       GOTO 990
    END IF
    ALLOCATE(this%RelativeAzimuthAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Relative Azimuth Angle'
       GOTO 990
    END IF
    ALLOCATE(this%SolarAzimuthAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Azimuth Angle'
       GOTO 990
    END IF
    ALLOCATE(this%SolarZenithAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Solar Zenith Angle'
       GOTO 990
    END IF
    ALLOCATE(this%ViewingAzimuthAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Viewing Azimuth Angle'
       GOTO 990
    END IF
    ALLOCATE(this%ViewingZenithAngle(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Viewing Zenith Angle'
       GOTO 990
    END IF
    ALLOCATE(this%ColumnAmountO3(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Column Amount O3'
       GOTO 990
    END IF
    ALLOCATE(this%MeasurementQualityFlags(this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Measurement Quality Flags'
       GOTO 990
    END IF
    ALLOCATE(this%QualityFlags(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Quality Flags'
       GOTO 990
    END IF
    ALLOCATE(this%RadiativeCloudFraction(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Radiatice Cloud Fraction'
       GOTO 990
    END IF
    ALLOCATE(this%Reflectivity331(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Reflectivity 331'
       GOTO 990
    END IF
    ALLOCATE(this%Reflectivity360(this%nXtrack,this%nLines), STAT = ierr)
    IF ( ierr .NE. 0 ) THEN
       ErrorMessage = 'Error allocating Reflectivity 360'
       GOTO 990
    END IF

    ! ============
    ! Reading data
    ! ============
    ! ----------------
    ! Ancillary data
    ! ----------------
    CALL H5ReadDataset(filename, &
         "/AncillaryData/CloudPressure", &
         this%CloudPressure)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/AncillaryData/SurfaceCategory", &
         this%SurfaceCategory)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/AncillaryData/TerrainPressure", &
         this%TerrainPressure)
    IF (ErrorFlag.NE.0) GOTO 990

    ! ----------------
    ! Geolocation data
    ! ----------------
    CALL H5ReadDataset(filename, &
         "/GeolocationData/Latitude", &
         this%Latitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/LatitudeCorner", &
         this%LatitudeCorner)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/Longitude", &
         this%Longitude)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/LongitudeCorner", &
         this%LongitudeCorner)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/RelativeAzimuthAngle", &
         this%RelativeAzimuthAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/SolarAzimuthAngle", &
         this%SolarAzimuthAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/SolarZenithAngle", &
         this%SolarZenithAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/ViewingAzimuthAngle", &
         this%ViewingAzimuthAngle)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/GeolocationData/ViewingZenithAngle", &
         this%ViewingZenithAngle)
    IF (ErrorFlag.NE.0) GOTO 990

    ! ------------
    ! Science data
    ! ------------
    CALL H5ReadDataset(filename, &
         "/ScienceData/ColumnAmountO3", &
         this%ColumnAmountO3)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/ScienceData/MeasurementQualityFlags", &
         this%MeasurementQualityFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/ScienceData/QualityFlags", &
         this%QualityFlags)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/ScienceData/RadiativeCloudFraction", &
         this%RadiativeCloudFraction)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/ScienceData/Reflectivity331", &
         this%Reflectivity331)
    IF (ErrorFlag.NE.0) GOTO 990
    CALL H5ReadDataset(filename, &
         "/ScienceData/Reflectivity360", &
         this%Reflectivity360)
    IF (ErrorFlag.NE.0) GOTO 990

    status = ErrorFlag
    RETURN
990 WRITE(*, '(A)') ErrorMessage
    status = 2
    RETURN
  END FUNCTION omps_nmto3_reader


END MODULE OMSAO_omps_reader
