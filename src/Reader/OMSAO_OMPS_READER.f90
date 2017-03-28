!! MODULE to READ OMPS data
!!
! Created by gga January 2014
!!
! All he5 internals are done by the belgioum (BIRA) module ReadH5dataset
! with small modifications
!!

MODULE OMSAO_OMPS_READER

  USE ReadH5dataset ! Contains the generic HDF5 reading routines

  IMPLICIT NONE
  
  PUBLIC :: TC_SDR_OMPS_READER


  INTEGER (KIND = 4), PARAMETER, PUBLIC :: MAX_NAME_LENGTH = 256

  TYPE, PUBLIC :: TC_SDR_OMPS_type
     CHARACTER (LEN=MAX_NAME_LENGTH) :: filename
     INTEGER*4                       :: nLines
     INTEGER*4                       :: nXtrack
     INTEGER*4                       :: nWavel
     
     ! Data storage variables
     ! =================================================================
     ! Ancillary data
     REAL*4, DIMENSION(:,:),     POINTER :: CloudPressure => NULL()
     REAL*4, DIMENSION(:,:,:,:), POINTER :: O3Profile => NULL()
     REAL*4, DIMENSION(:),       POINTER :: PressureProfile => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SnowIceFraction => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: StandardDeviation => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SurfaceReflectivity => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: TemperatureProfile => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: TerrainPressure => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: TropopausePressure => NULL()
     REAL*4, DIMENSION(:,:,:,:), POINTER :: TroposphericO3 => NULL()

     ! Calibration data
     REAL*4, DIMENSION(:,:,:),   POINTER :: BandCenterWavelengths => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: CCDRowColIndicies => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: DarkCurrentCorrection => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: RadianceCalCoeff => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: SmearCorrection => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SolarFlux => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SolarFluxWavelengths => NULL()
     
     ! Geolocation data
     REAL*8, DIMENSION(:),       POINTER :: GoniometricSolarAzimuth => NULL()
     REAL*8, DIMENSION(:),       POINTER :: GoniometricSolarElevation => NULL()
     INTEGER*4, DIMENSION(:,:),  POINTER :: GroundPixelQualityFlags => NULL()
     REAL*8, DIMENSION(:),       POINTER :: ImageMidpoint_TAI93 => NULL()
     INTEGER*4, DIMENSION(:),    POINTER :: InstrumentQualityFlags => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: Latitude => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: Longitude => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SatelliteAzimuth => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SatelliteZenithAngle => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SolarAzimuth => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SolarBeta => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SolarDeclination_ECI => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SolarRightAscension_ECI => NULL()
     REAL*8, DIMENSION(:,:),     POINTER :: SolarUnitVectorECI => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: SolarZenithAngle => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SpacecraftAltitude => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SpacecraftLatitude => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SpacecraftLongitude => NULL()
     REAL*8, DIMENSION(:,:),     POINTER :: SpacecraftPositionECI => NULL()
     REAL*8, DIMENSION(:,:),     POINTER :: SpacecraftVelocityECI => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SubSatelliteSolarZenithAngle => NULL()
     REAL*8, DIMENSION(:),       POINTER :: SunEarthDistance => NULL()
     CHARACTER(LEN=MAX_NAME_LENGTH), &
          DIMENSION(:), POINTER :: UTC_CCSDS_A => NULL()
     REAL*4, DIMENSION(:),       POINTER :: expose => NULL()

     ! Input pointers
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_L1A_EV
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_STB
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_mCBC
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_mDRK
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_mIRF
     CHARACTER(LEN=MAX_NAME_LENGTH) :: TC_mRAD
     
     ! Science data
     REAL*4, DIMENSION(:,:),     POINTER :: ExposureTime => NULL()
     REAL*4, DIMENSION(:,:),     POINTER :: HousekeepingData => NULL()
     INTEGER*4, DIMENSION(:,:,:),POINTER :: IntegrationTimeUsed => NULL()
     INTEGER*4, DIMENSION(:,:),  POINTER :: NumberCoadds => NULL()
     INTEGER*4, DIMENSION(:,:,:),POINTER :: PixelQualityFlags => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: Radiance => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: RadianceError => NULL()
     REAL*4, DIMENSION(:,:,:),   POINTER :: RawCounts => NULL()
     INTEGER*4, DIMENSION(:),    POINTER :: ReportIntQualFlags => NULL()
     INTEGER*4, DIMENSION(:),    POINTER :: SensorStatusBits => NULL()

     END TYPE TC_SDR_OMPS_type

     CONTAINS
       FUNCTION TC_SDR_OMPS_READER(this, filename) RESULT(status)

         ! I'm going to use the Belgiums hdf5 reader
         ! The most recent version of the package is available at:
         ! http://uv-vis.aeronomie.be/software/tools/hdf5read.php
         ! All information available on the OMPS-NPP-NPP-TC_SDR_EV_NASA
         ! type files is going to be readed for possible future use.
         ! Stand alone version.

         CHARACTER (LEN=MAX_NAME_LENGTH)        :: filename
         TYPE (TC_SDR_OMPS_type), INTENT(INOUT) :: this
         INTEGER*2                              :: status

         ! ----------------
         ! Code starts here
         ! ----------------
         ! Save filename in the TC_SDR_OMPS_type
         this%filename = filename

         ! ----------------------------
         ! Read elements of the h5 file
         ! ----------------------------
         ! Ancillary data
         ! --------------
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/CloudPressure", this%CloudPressure)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/O3Profile", this%O3Profile)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/PressureProfile", this%PressureProfile)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/SnowIceFraction", this%SnowIceFraction)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/StandardDeviation", this%StandardDeviation)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/SurfaceReflectivity", this%SurfaceReflectivity)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/TemperatureProfile", this%TemperatureProfile)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/TerrainPressure", this%TerrainPressure)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/TropopausePressure", this%TropopausePressure)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/ANC_DATA/TroposphericO3", this%TroposphericO3)
         IF (ErrorFlag.lt.0) goto 990
         ! ----------------
         ! Calibration data
         ! ----------------
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/BandCenterWavelengths", this%BandCenterWavelengths)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/CCDRowColIndicies", this%CCDRowColIndicies)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/DarkCurrentCorrection", this%DarkCurrentCorrection)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/RadianceCalCoeff", this%RadianceCalCoeff)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/SmearCorrection", this%SmearCorrection)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/SolarFlux", this%SolarFlux)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/CALIBRATION_DATA/SolarFluxWavelengths", this%SolarFluxWavelengths)
         IF (ErrorFlag.lt.0) goto 990
         ! ----------------
         ! Geolocation data
         ! ----------------
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/GoniometricSolarAzimuth", this%GoniometricSolarAzimuth)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/GoniometricSolarElevation", this%GoniometricSolarElevation)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/GroundPixelQualityFlags", this%GroundPixelQualityFlags)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/ImageMidpoint_TAI93", this%ImageMidpoint_TAI93)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/InstrumentQualityFlags", this%InstrumentQualityFlags)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/Latitude", this%Latitude)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/Longitude", this%Longitude)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SatelliteAzimuth", this%SatelliteAzimuth)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SatelliteZenithAngle", this%SatelliteZenithAngle)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarAzimuth", this%SolarAzimuth)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarBeta", this%SolarBeta)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarDeclination_ECI", this%SolarDeclination_ECI)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarRightAscension_ECI", this%SolarRightAscension_ECI)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarUnitVectorECI", this%SolarUnitVectorECI)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SolarZenithAngle", this%SolarZenithAngle)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SpacecraftAltitude", this%SpacecraftAltitude)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SpacecraftLatitude", this%SpacecraftLatitude)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SpacecraftLongitude", this%SpacecraftLongitude)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SpacecraftPositionECI", this%SpacecraftPositionECI)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SpacecraftVelocityECI", this%SpacecraftVelocityECI)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SubSatelliteSolarZenithAngle", this%SubSatelliteSolarZenithAngle)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/SunEarthDistance", this%SunEarthDistance)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/UTC_CCSDS_A", this%UTC_CCSDS_A)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/GEOLOCATION_DATA/expose", this%expose)
         IF (ErrorFlag.lt.0) goto 990
         ! --------------
         ! Input pointers
         ! --------------
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_L1A_EV", this%TC_L1A_EV)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_STB", this%TC_STB)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_mCBC", this%TC_mCBC)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_mDRK", this%TC_mDRK)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_mIRF", this%TC_mIRF)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/InputPointers/TC_mRAD", this%TC_mRAD)
         IF (ErrorFlag.lt.0) goto 990       
         ! ------------
         ! Science data
         ! ------------
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/ExposureTime", this%ExposureTime)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/HousekeepingData", this%HousekeepingData)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/IntegrationTimeUsed", this%IntegrationTimeUsed)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/NumberCoadds", this%NumberCoadds)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/PixelQualityFlags", this%PixelQualityFlags)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/Radiance", this%Radiance)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/RadianceError", this%RadianceError)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/RawCounts", this%RawCounts)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/ReportIntQualFlags", this%ReportIntQualFlags)
         IF (ErrorFlag.lt.0) goto 990
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/SensorStatusBits", this%SensorStatusBits)
         IF (ErrorFlag.lt.0) goto 990

         ! ----------------------------------------------------------
         ! Filling up the dimension variables nLines, nXtrack, nWavel
         ! ----------------------------------------------------------
         CALL H5ReadDataset(filename, &
              "/SCIENCE_DATA/Radiance", this%nLines, this%nXtrack, this%nWavel)
         IF (ErrorFlag.lt.0) goto 990 

 
990      WRITE(*, '(A)') ErrorMessage
         status = ErrorFlag
         RETURN

       END FUNCTION TC_SDR_OMPS_READER
     END MODULE OMSAO_OMPS_READER
