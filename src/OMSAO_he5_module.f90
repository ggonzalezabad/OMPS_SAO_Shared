MODULE OMSAO_he5_module
  
  ! ==============================
  ! Module for HDF-EOS5 parameters
  ! ==============================

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: &
       maxchlen, i2_missval, i4_missval, r4_missval, r8_missval, blank13,       &
       valid_max_i2, valid_max_i4, valid_min_r8, valid_max_r8, zero_r8, one_r8, &
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8,                          &
       main_qa_min_flag_r8, main_qa_max_flag_r8
  USE OMSAO_indices_module,    ONLY: &
       n_sao_pge, max_calfit_idx, sao_pge_min_idx, sao_pge_max_idx, &
       o3_t1_idx, o3_t3_idx
  USE HDF5

  IMPLICIT NONE

  CHARACTER (LEN=12), PARAMETER :: omspec = "OMI-Specific"

  ! -----------------------------------------------------------------
  ! Blank strings of various lengths.
  !
  ! Some compilers don't allow to define CHARACTER PARAMETER arrays
  ! which are initialized with field of unequal length. The following
  ! are a few padding strings that we attach to shorter entries. Not
  ! very stylish, but effective. Note that what matters is the LEN
  ! declaration - we don't need to initialize the strings with the
  ! appropriate number of blanks. Using "" is perfectly fine.
  ! -----------------------------------------------------------------
  CHARACTER (LEN=12), PARAMETER :: blk12 = ""
  CHARACTER (LEN=13), PARAMETER :: blk13 = ""
  CHARACTER (LEN=14), PARAMETER :: blk14 = ""
  CHARACTER (LEN=16), PARAMETER :: blk16 = ""
  CHARACTER (LEN=17), PARAMETER :: blk17 = ""
  CHARACTER (LEN=19), PARAMETER :: blk19 = ""
  CHARACTER (LEN=20), PARAMETER :: blk20 = ""
  CHARACTER (LEN=21), PARAMETER :: blk21 = ""
  CHARACTER (LEN=23), PARAMETER :: blk23 = ""
  CHARACTER (LEN=25), PARAMETER :: blk25 = ""
  CHARACTER (LEN=26), PARAMETER :: blk26 = ""
  CHARACTER (LEN=27), PARAMETER :: blk27 = ""
  CHARACTER (LEN=28), PARAMETER :: blk28 = ""
  CHARACTER (LEN=29), PARAMETER :: blk29 = ""
  CHARACTER (LEN=30), PARAMETER :: blk30 = ""
  CHARACTER (LEN=31), PARAMETER :: blk31 = ""
  CHARACTER (LEN=32), PARAMETER :: blk32 = ""
  CHARACTER (LEN=33), PARAMETER :: blk33 = ""
  CHARACTER (LEN=35), PARAMETER :: blk35 = ""
  CHARACTER (LEN=38), PARAMETER :: blk38 = ""
  CHARACTER (LEN=39), PARAMETER :: blk39 = ""
  CHARACTER (LEN=40), PARAMETER :: blk40 = ""
  CHARACTER (LEN=41), PARAMETER :: blk41 = ""
  CHARACTER (LEN=43), PARAMETER :: blk43 = ""
  CHARACTER (LEN=45), PARAMETER :: blk45 = ""
  CHARACTER (LEN=48), PARAMETER :: blk48 = ""
  CHARACTER (LEN=49), PARAMETER :: blk49 = ""

  ! ----------------------------------------------
  ! Swath IDs for current HE5 swath and swath file
  ! ----------------------------------------------
  INTEGER (KIND=i4) :: pge_swath_file_id, pge_swath_id

  ! -----------------------------------------------
  ! Variables that will hold Global File Attributes
  ! -----------------------------------------------
  CHARACTER (LEN=maxchlen) :: &
       pge_swath_name, process_level, instrument_name, pge_version, l1b_orbitdata

  ! --------------------------------------------------------------
  ! Swath IDs for pre-fitted BrO and O3 HE5 swaths and swath files
  ! --------------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: &
       prefit_swath_name    = 'undefined'
  INTEGER   (KIND=i4) :: &
       prefit_swath_id    = -1, prefit_swath_file_id    = -1

  ! --------------------
  ! Fill Value attribute
  ! --------------------
  CHARACTER (LEN= 9), PARAMETER :: fillval_attr = "FillValue"
  CHARACTER (LEN=12), PARAMETER :: missval_attr = "MissingValue"
  CHARACTER (LEN= 6), PARAMETER :: offset_attr  = "Offset"
  CHARACTER (LEN=11), PARAMETER :: scafac_attr  = "ScaleFactor"
  CHARACTER (LEN= 5), PARAMETER :: title_attr   = "Title"
  CHARACTER (LEN=10), PARAMETER :: valids_attr  = "ValidRange"
  CHARACTER (LEN= 5), PARAMETER :: units_attr   = "Units"
  CHARACTER (LEN=28), PARAMETER :: rstemp_attr  = "ReferenceSpectrumTemperature"
  CHARACTER (LEN=21), PARAMETER :: ufd_attr     = "UniqueFieldDefinition"

  ! ----------------------------------
  ! Strings for Global File Attributes
  ! ----------------------------------
  INTEGER   (KIND=i4)      :: granule_day, granule_month, granule_year

  ! ----------------------
  ! Swath Level Attributes
  ! ----------------------
  CHARACTER (LEN=18), PARAMETER :: vcoordinate_field = "VerticalCoordinate"

  ! ----------------
  ! Swath dimensions
  ! ----------------
  CHARACTER (LEN= 6), PARAMETER :: ntc   = "nTimes"
  CHARACTER (LEN= 7), PARAMETER :: nxc   = "nXtrack"
  CHARACTER (LEN= 7), PARAMETER :: nlc   = "nLevels"
  CHARACTER (LEN=12), PARAMETER :: nfv   = "nFitElements"
  CHARACTER (LEN=19), PARAMETER :: ncv   = "nCharLenFitElements"
  CHARACTER (LEN= 7), PARAMETER :: nutcd = "nUTCdim"
  CHARACTER (LEN=11), PARAMETER :: nwcp  = "nWavCalPars"
  CHARACTER (LEN=14), PARAMETER :: nxtc  = nxc//","//ntc
  CHARACTER (LEN=27), PARAMETER :: nfxtc = nfv//","//nxc//","//ntc
  CHARACTER (LEN=19), PARAMETER :: nxwcp = nwcp//","//nxc
  CHARACTER (LEN=22), PARAMETER :: nxtlc = nxc//","//ntc//","//nlc !GGA

  ! CCM Add one for the refspec database
  CHARACTER (LEN=7), PARAMETER :: nrspc = "nRfSpec"
  CHARACTER (LEN=6), PARAMETER :: nwalm = "nWaves"

  ! ------------------------
  ! Swath geolocation fields
  ! ------------------------
  CHARACTER (LEN= 8), PARAMETER ::  lat_field   = "Latitude"
  CHARACTER (LEN= 9), PARAMETER ::  lon_field   = "Longitude"
  CHARACTER (LEN=14), PARAMETER ::  latcor_field= "LatitudeCorner"
  CHARACTER (LEN=15), PARAMETER ::  loncor_field= "LongitudeCorner"
  CHARACTER (LEN=17), PARAMETER ::  saa_field   = "SolarAzimuthAngle"
  CHARACTER (LEN=16), PARAMETER ::  sza_field   = "SolarZenithAngle"
  CHARACTER (LEN=19), PARAMETER ::  vaa_field   = "ViewingAzimuthAngle"
  CHARACTER (LEN=18), PARAMETER ::  vza_field   = "ViewingZenithAngle"
  CHARACTER (LEN=20), PARAMETER ::  raa_field   = "RelativeAzimuthAngle"
  CHARACTER (LEN=23), PARAMETER ::  xtr_field   = "GroundPixelQualityFlags"
  CHARACTER (LEN=22), PARAMETER ::  extr_field  = "InstrumentQualityFlags"
  CHARACTER (LEN=18), PARAMETER ::  alt_field   = "SpacecraftAltitude"
  CHARACTER (LEN= 5), PARAMETER ::  time_field  = "TAI93"
  CHARACTER (LEN= 7), PARAMETER ::  utc_field   = "TimeUTC"

  ! -----------------
  ! Swath data fields
  ! -----------------
  CHARACTER (LEN=27), PARAMETER ::  amfdiag_field  = "AirMassFactorDiagnosticFlag"
  CHARACTER (LEN=22), PARAMETER ::  amfgeo_field   = "AirMassFactorGeometric"
  CHARACTER (LEN=13), PARAMETER ::  amfmol_field   = "AirMassFactor"
  CHARACTER (LEN=19), PARAMETER ::  avgcol_field   = "AverageColumnAmount"
  CHARACTER (LEN=24), PARAMETER ::  avgdcol_field  = "AverageColumnUncertainty"
  CHARACTER (LEN=17), PARAMETER ::  avgrms_field   = "AverageFittingRMS"
  CHARACTER (LEN=12), PARAMETER ::  col_field      = "ColumnAmount"
  CHARACTER (LEN=18), PARAMETER ::  commspc_field  = "CommonModeSpectrum"
  CHARACTER (LEN=21), PARAMETER ::  commwvl_field  = "CommonModeWavelengths"
  CHARACTER (LEN=17), PARAMETER ::  dcol_field     = "ColumnUncertainty"
  CHARACTER (LEN=21), PARAMETER ::  dstrcol_field  = "ColumnAmountDestriped"
  CHARACTER (LEN=18), PARAMETER ::  fitcon_field   = "FitConvergenceFlag"
  CHARACTER (LEN=06), PARAMETER ::  fitrms_field   = "FitRMS"
  CHARACTER (LEN=19), PARAMETER ::  mainqa_field   = "MainDataQualityFlag"
  CHARACTER (LEN=20), PARAMETER ::  pxclat_field   = "PixelCornerLatitudes"
  CHARACTER (LEN=21), PARAMETER ::  pxclon_field   = "PixelCornerLongitudes"

  ! CCM Add field for fit residuals
  CHARACTER (LEN=11), PARAMETER ::  spcfit_field   = "FitSpectrum"
  CHARACTER (LEN=19), PARAMETER ::  spcobs_field   = "FitMeasuredSpectrum"
  CHARACTER (LEN=21), PARAMETER ::  posobs_field   = "FitMeasuredWavelength"
  CHARACTER (LEN=20), PARAMETER ::  fitwt_field    = "FitSpectralFitWeight"
  CHARACTER (LEN=19), PARAMETER ::  spcres_field   = "FitResidualSpectrum"
  ! Spectral Database
  CHARACTER (LEN=12), PARAMETER ::  spdata_field   = "DatabaseSpec"
  CHARACTER (LEN=12), PARAMETER ::  spdatw_field   = "DatabaseWavl"
  CHARACTER (LEN=17), PARAMETER ::  spname_field   = "DatabaseSpecNames"
  CHARACTER (LEN=18), PARAMETER ::  spnrmf_field   = "DatabaseNormFactor"

  ! -------------------------------------------------------
  ! Special data fields for wavelength-modified AMF fitting
  ! -------------------------------------------------------
  CHARACTER (LEN=20), PARAMETER ::  adalb_field     = "AdjustedSceneAlbedo"
  CHARACTER (LEN=18), PARAMETER ::  scol_field      = "SlantColumnAmount"
  CHARACTER (LEN=22), PARAMETER ::  sdcol_field     = "SlantColumnUncertainty"
  CHARACTER (LEN=26), PARAMETER ::  sdstrcol_field  = "SlantColumnAmountDestriped"
  CHARACTER (LEN=23), PARAMETER ::  sfitcon_field   = "SlantFitConvergenceFlag"
  CHARACTER (LEN=15), PARAMETER ::  sfitrms_field   = "SlantFittingRMS"


  CHARACTER (LEN=29), PARAMETER ::  rrcol_field    = "RadianceReferenceColumnAmount"
  CHARACTER (LEN=34), PARAMETER ::  rrdcol_field   = "RadianceReferenceColumnUncertainty"
  CHARACTER (LEN=29), PARAMETER ::  rrxcol_field   = "RadianceReferenceColumnXTRFit"
  CHARACTER (LEN=27), PARAMETER ::  rrrms_field    = "RadianceReferenceFittingRMS"
  CHARACTER (LEN=20), PARAMETER ::  rrfi_field     = "RadianceReferenceFit"
  CHARACTER (LEN=32), PARAMETER ::  rrcf_field     = "RadianceReferenceConvergenceFlag"
  CHARACTER (LEN=31), PARAMETER ::  rric_field     = "RadianceReferenceIterationCount"
  CHARACTER (LEN=30), PARAMETER ::  rrlr_field     = "RadianceReferenceLatitudeRange"

  CHARACTER (LEN=29), PARAMETER ::  rwccf_field    = "RadianceWavCalConvergenceFlag"
  CHARACTER (LEN=18), PARAMETER ::  rwwav_field    = "RadianceWavCalWavl"
  CHARACTER (LEN=18), PARAMETER ::  rwrad_field    = "RadianceWavCalRadi"
  CHARACTER (LEN=18), PARAMETER ::  rwwei_field    = "RadianceWavCalWght"
  CHARACTER (LEN=18), PARAMETER ::  rwres_field    = "RadianceWavCalResi"
  CHARACTER (LEN=18), PARAMETER ::  rwshi_field    = "RadianceWavCalShif"
  CHARACTER (LEN=18), PARAMETER ::  rwsqu_field    = "RadianceWavCalSque"

  CHARACTER (LEN=26), PARAMETER ::  swccf_field    = "SolarWavCalConvergenceFlag"
  CHARACTER (LEN=15), PARAMETER ::  swwav_field    = "SolarWavCalWavl"
  CHARACTER (LEN=15), PARAMETER ::  swirr_field    = "SolarWavCalIrra"
  CHARACTER (LEN=15), PARAMETER ::  swwei_field    = "SolarWavCalWght"
  CHARACTER (LEN=15), PARAMETER ::  swres_field    = "SolarWavCalResi"
  CHARACTER (LEN=15), PARAMETER ::  swshi_field    = "SolarWavCalShif"
  CHARACTER (LEN=15), PARAMETER ::  swsqu_field    = "SolarWavCalSque"
  CHARACTER (LEN=15), PARAMETER ::  swhw1_field    = "SolarWavCalHw1e"
  CHARACTER (LEN=15), PARAMETER ::  swasy_field    = "SolarWavCalAsym"
  CHARACTER (LEN=15), PARAMETER ::  swsha_field    = "SolarWavCalShap"

  ! --------------------------------------------------
  ! Swath data fields, additional in "diagnostic" runs
  ! --------------------------------------------------
  CHARACTER (LEN=13), PARAMETER ::  ccdpix_field   = "CCDPixelRange"
  CHARACTER (LEN=15), PARAMETER ::  commcnt_field  = "CommonModeCount"
  CHARACTER (LEN=24), PARAMETER ::  corr_field     = "FitParameterCorrelations"
  CHARACTER (LEN=19), PARAMETER ::  corrcol_field  = "FitParameterColumns"
  CHARACTER (LEN=17), PARAMETER ::  correlm_field  = "FitParameterNames"
  CHARACTER (LEN=23), PARAMETER ::  correrr_field  = "FitParameterUncertainty"
  CHARACTER (LEN=17), PARAMETER ::  itnum_field    = "FitIterationCount"
  CHARACTER (LEN=26), PARAMETER ::  xtrcor_field   = "CrossTrackStripeCorrection"


  ! -----------------------------------------------
  ! Swath data fields unique to OMHCHO and OMCHOCHO
  ! -----------------------------------------------
  CHARACTER (LEN=16), PARAMETER ::  amfcfr_field   = "AMFCloudFraction"
  CHARACTER (LEN=16), PARAMETER ::  amfctp_field   = "AMFCloudPressure"

  ! -------------------------------------- GGA
  ! Swath data fields for Reference Sector GGA
  ! -------------------------------------- GGA
  CHARACTER (LEN=21), PARAMETER :: rscol_field     = "ColumnAmountCorrected"
  CHARACTER (LEN=25), PARAMETER :: rsfla_field     = "ColumnAmountCorrectedFlag"

  ! ---------------------------------------------------------------------- GGA
  ! Swath data field for Scattering Weights, Gas Profile Averaging Kernels GGA
  ! and albedo                                                             GGA
  ! ---------------------------------------------------------------------- GGA
  CHARACTER (LEN=17), PARAMETER :: scaweights_field = "ScatteringWeights"
  CHARACTER (LEN=10), PARAMETER :: gasprofile_field = "GasProfile"
  CHARACTER (LEN=19), PARAMETER :: clialtgrid_field = "ClimatologyLevels"
  CHARACTER (LEN= 6), PARAMETER :: albedo_field     = "Albedo"

  ! ------------------
  ! Geolocation fields
  ! ------------------
  INTEGER (KIND=i4), PARAMETER :: n_gfields = 14
  CHARACTER (LEN=41), DIMENSION (2,n_gfields), PARAMETER ::  &
       geo_field_names = RESHAPE ( (/ &
       "Latitude                                  ","Geodetic Latitude                         ",    &
       "Longitude                                 ","Geodetic Longitude                        ",    &
       "LatitudeCorner                            ","Geodetic Latitude of Corner Points        ",    &
       "LongitudeCorner                           ","Geodetic Longitude of Corner Points       ",    &
       "SolarAzimuthAngle                         ","Solar Azimuth Angle                       ",    &
       "SolarZenithAngle                          ","Solar Zenith Angle                        ",    &
       "ViewingAzimuthAngle                       ","Viewing Azimuth Angle                     ",    &
       "ViewingZenithAngle                        ","Viewing Zenith Angle                      ",    &
       "RelativeAzimuthAngle                      ","Relative Azimuth Angle                    ",    &
       "GroundPixelQualityFlags                   ","Pixel Level Geolocation Data Quality Flags",    &
       "InstrumentQualityFlags                    ","Instrument Quality Flags                  ",    &
       "SpacecraftAltitude                        ","Altitude of Spacecraft                    ",    &
       "TAI93                                     ","TAI93 Image Midpoint Time                 ",    &
       "TimeUTC                                   ","UTC Image Midpoint Time                   "  /),&
       (/ 2, n_gfields /) )
  CHARACTER (LEN=16), DIMENSION (3,n_gfields), PARAMETER ::  &
       geo_field_specs = RESHAPE ( (/ &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","4,nXtrack,nTimes","r4              ",    &
       "deg             ","4,nXtrack,nTimes","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "deg             ","nXtrack,nTimes  ","r4              ",    &
       "                ","nXtrack,nTimes  ","i2              ",    &
       "                ","nTimes          ","i4              ",    &
       "m               ","nTimes          ","r4              ",    &
       "Seconds         ","nTimes          ","r8              ",    &
       "                ","nTimes,nUTCdim  ","ch              "  /),&
       (/ 3, n_gfields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_gfields ), PARAMETER :: &
       geo_valids = RESHAPE (       (/  &
        -90.0_r8,     +90.0_r8,     & !Latitude
       -180.0_r8,   +180.0_r8,     & !Longitude
        -90.0_r8,     +90.0_r8,     & !Latitude corner
       -180.0_r8,   +180.0_r8,     & !Longitude corner
       -180.0_r8,   +180.0_r8,     & !SAA
         zero_r8,     +180.0_r8,     & !SZA
       -180.0_r8,   +180.0_r8,     & !VAA
         zero_r8,     +180.0_r8,     & !VZA
       -180.0_r8,   +180.0_r8,     & !RAA
          1.0_r8, 1.0_r8,       & !Ground Pixel Quality Falgs
          1.0_r8, 1.0_r8,       & !Instrument Quality Flags
       400000.0_r8, 900000.0_r8,       & !Spacecraft altitutude
       valid_min_r8, valid_max_r8, & !TAI93
       zero_r8,   zero_r8      /), & !UTCTime
       (/ 2, n_gfields /) )

  ! --------------------------------------------  
  ! Data fields for Solar Wavelength Calibration
  ! --------------------------------------------  
  INTEGER (KIND=i4), PARAMETER :: n_solcal_fields = 10
  CHARACTER (LEN=45), DIMENSION ( 2, n_solcal_fields ), PARAMETER :: &
       solcal_field_names = RESHAPE ( (/ &
       "SolarWavCalConvergenceFlag                   ","Solar Wavelength Calibration Convergence Flag", &
       "SolarWavCalWavl                              ","Solar Wavelenght                             ", &
       "SolarWavCalIrra                              ","Solar Irradiance                             ", &
       "SolarWavCalWght                              ","Solar Weight                                 ", &
       "SolarWavCalResi                              ","Solar Residual                               ", &
       "SolarWavCalShif                              ","Solar Shift                                  ", &
       "SolarWavCalHw1e                              ","Solar Half Width at 1/e                      ", &
       "SolarWavCalSque                              ","Solar Squeeze                                ", &
       "SolarWavCalAsym                              ","Solar Asymmetry Parameter                    ", &
       "SolarWavCalShap                              ","Solar Shape Factor                           "/), &
       (/ 2, n_solcal_fields /) )
  CHARACTER (LEN=18), DIMENSION ( 3, n_solcal_fields ), PARAMETER :: &
       solcal_field_specs = RESHAPE ( (/ &
       "NoUnits           ","nXtrack           ","i2                " ,&
       "NoUnits           ","nXtrack,nWaves    ","r8                " ,&
       "NoUnits           ","nXtrack,nWaves    ","r8                " ,&
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                "/), &
       (/ 3, n_solcal_fields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_solcal_fields ), PARAMETER :: &
       solcal_valids = RESHAPE ( (/ &
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8, &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8           /), &
       (/ 2, n_solcal_fields /) )

  ! -----------------------------------------------  
  ! Data fields for Radiance Wavelength Calibration
  ! ----------------------------------------------------------------------
  ! NOTE: First N_SOLCAL_FIELDS must be the same type as in SOLCAL_FIELDS,
  !       since this is what is assumed in the output routine.
  ! ----------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_radcal_fields = 7
  CHARACTER (LEN=48), DIMENSION ( 2, n_radcal_fields ), PARAMETER :: &
       radcal_field_names = RESHAPE ( (/ &
       "RadianceWavCalConvergenceFlag                   ", "Radiance Wavelength Calibration Convergence Flag", &
       "RadianceWavCalRadi                              ", "Radiance Wavelength Calibration Radiance        ", &
       "RadianceWavCalWght                              ", "Radiance Wavelength Calibration Weight          ", &
       "RadianceWavCalResi                              ", "Radiance Wavelength Calibration Residual        ", &
       "RadianceWavCalWavl                              ", "Radiance Wavelenght Calibration Wavelenght      ", &
       "RadianceWavCalShif                              ", "Radiance Wavelenght Calibration Shift           ", &
       "RadianceWavCalSque                              ", "Radiance Wavelenght Calibration Squeeze         "/), &
       (/ 2, n_radcal_fields /) )
  CHARACTER (LEN=18), DIMENSION ( 3, n_radcal_fields ), PARAMETER :: &
       radcal_field_specs = RESHAPE ( (/ &
       "NoUnits           ","nXtrack           ","i2                ", &
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack,nWaves    ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                ", &
       "NoUnits           ","nXtrack           ","r8                " /), &
       (/ 3, n_radcal_fields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_radcal_fields ), PARAMETER :: &
       radcal_valids = RESHAPE ( (/ &
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8, &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8          , &
       valid_min_r8          , valid_max_r8           /), &
       (/ 2, n_radcal_fields /) )

  ! -----------------------------------------------  
  ! Data fields for Radiance Reference Fit
  ! ----------------------------------------------------------------------
  ! NOTE: First N_SOLCAL_FIELDS must be the same type as in SOLCAL_FIELDS,
  !       since this is what is assumed in the output routine.
  ! ----------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_radref_fields = 6
  CHARACTER (LEN=42), DIMENSION ( 2, n_radref_fields ), PARAMETER :: &
       radref_field_names = RESHAPE ( (/ &
       "RadianceReferenceConvergenceFlag         ", "Radiance Reference Fit Convergence Flag  ",    &
       "RadianceReferenceLatitudeRange           ", "Radiance Reference Fit Latitude Range    ",    &
       "RadianceReferenceColumnAmount            ", "Radiance Reference Fit Colunm Amount     ",    &
       "RadianceReferenceColumnUncertainty       ", "Radiance Reference Fit Colunm Uncertainty",    &
       "RadianceReferenceColumnXTRFit            ", "Radiance Reference Fit Colunm XTR Fit    ",    &
       "RadianceReferenceFittingRMS              ", "Radiance Reference Fit RMS               "  /),&
       (/ 2, n_radref_fields /) )
  CHARACTER (LEN=12), DIMENSION ( 3, n_radref_fields ), PARAMETER :: &
       radref_field_specs = RESHAPE ( (/ &
       "NoUnits     ", "nXtrack     ", "i2          ",   &
       "NoUnits     ", "2           ", "r4          ",   &
       "molec/cm2   ", "nXtrack     ", "r8          ",   &
       "molec/cm2   ", "nXtrack     ", "r8          ",   &
       "molec/cm2   ", "nXtrack     ", "r8          ",   &
       "NoUnits     ", "nXtrack     ", "r8          " /),&
       (/ 3, n_radref_fields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_radref_fields ), PARAMETER :: &
       radref_valids = RESHAPE ( (/ &
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8,    &
       -90.0_r8,               90.0_r8,                   &
       valid_min_r8,           valid_max_r8,              &
       zero_r8,                valid_max_r8,              &
       valid_min_r8,           valid_max_r8,              &
       zero_r8,                valid_max_r8           /), &
       (/ 2, n_radref_fields /) )

  ! -----------
  ! Data Fields
  ! -----------
  ! -----------------------------
  ! Main common output quantities
  ! -----------------------------
  INTEGER (KIND=i4), PARAMETER :: n_cdfields = 17
  CHARACTER (LEN=41), DIMENSION (2,n_cdfields), PARAMETER ::  &
       comdata_field_names = RESHAPE ( (/ &
       "AirMassFactor                            ", "Molecule Specific Air Mass Factor (AMF)  ", &
       "AirMassFactorDiagnosticFlag              ", "Diagnostic Flag for Molecule Specific AMF", &
       "AirMassFactorGeometric                   ", "Geometric Air Mass Factor (AMF)          ", &
       "AirMassFactorCloudFraction               ", "Cloud Fraction for AMF Computation       ", &
       "AirMassFactorCloudPressure               ", "Cloud Pressure for AMF Computation       ", &
       "Albedo                                   ", "Scene Albedo for AMF Computation         ", &
       "AprioriGasProfile                        ", "A-priori Gas profile for AMF Computation ", &
       "AverageColumnAmount                      ", "Average Column Amount                    ", &
       "AverageColumnUncertainty                 ", "Average Column Uncertainty               ", &
       "AverageFittingRMS                        ", "Average Fitting RMS                      ", &
       "ColumnAmount                             ", "Column Amount                            ", &
       "ColumnUncertainty                        ", "Column Uncertainty                       ", &
       "FitRMS                                   ", "Fitting RMS                              ", &
       "FitConvergenceFlag                       ", "Fitting Convergence Flag                 ", &
       "MainDataQualityFlag                      ", "Main Data Quality Flag                   ", &
       "ScatteringWeights                        ", "Scattering Weights                       ", &
       "SurfacePressure                          ", "Surface Pressure for AMF Computation     "/),&
       (/ 2, n_cdfields /) )
  CHARACTER (LEN=22), DIMENSION (3,n_cdfields), PARAMETER ::  &
       comdata_field_specs = RESHAPE ( (/ &
       "NoUnits               ", "nXtrack,nTimes        ","r4                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","i2                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","r4                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","r4                    ", &
       "hPa                   ", "nXtrack,nTimes        ","r4                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","r4                    ", &
       "molecules cm-2        ", "nXtrack,nTimes,nLevels","r8                    ", &
       "molec/cm2             ", "1                     ","r8                    ", &
       "molec/cm2             ", "1                     ","r8                    ", &
       "NoUnits               ", "1                     ","r8                    ", &
       "molec/cm2             ", "nXtrack,nTimes        ","r8                    ", &
       "molec/cm2             ", "nXtrack,nTimes        ","r8                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","r8                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","i2                    ", &
       "NoUnits               ", "nXtrack,nTimes        ","i2                    ", &
       "NoUnits               ", "nXtrack,nTimes,nLevels","r8                    ", &
       "hPa                   ", "nXtrack,nTimes        ","r8                    "/), &
       (/3, n_cdfields /) )

  ! -----------------------------------------------------------------
  ! Note that the VALIDS for FittingCorrelationElements are not used;
  ! they are only here to make this array of the same dimension as
  ! the others, and so eliminates the need for additional CASE logic.
  ! -----------------------------------------------------------------
  REAL (KIND=r8), DIMENSION ( 2, n_cdfields ), PARAMETER ::    &
       comdata_valids = RESHAPE ( (/                           &
       zero_r8,                valid_max_r8,                   &  ! AirMassFactor
       -2.0_r8,                13127.0_r8,                     &  ! AirMassFactorDiagnosticFlag
       zero_r8 ,               valid_max_r8,                   &  ! AirMassFactorGeometric
       zero_r8,                +1.0_r8,                        &  ! AMFCloudFraction
       zero_r8,                valid_max_r8,                   &  ! AMFCloudPressure
       valid_min_r8,         valid_max_r8,                     &  ! Albedo
       valid_min_r8,         valid_max_r8,                     &  ! Gas profile
       valid_min_r8,           valid_max_r8,                   &  ! AverageColumnAmount
       zero_r8,                valid_max_r8,                   &  ! AverageColumnUncertainty
       zero_r8,                valid_max_r8,                   &  ! AverageFittingRMS
       valid_min_r8,           valid_max_r8,                   &  ! ColumnAmount
       zero_r8,                valid_max_r8,                   &  ! ColumnUncertainty
       zero_r8,                valid_max_r8,                   &  ! FittingRMS
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8,         &  ! FitConvergenceFlag
       main_qa_min_flag_r8,    main_qa_max_flag_r8,            &  ! Main data quality flag
       valid_min_r8,         valid_max_r8,                     &  ! Scattering Weights
       valid_min_r8,         valid_max_r8                  /), &  ! Surface Pressure
       (/ 2, n_cdfields /) )

  ! ------------------------------------------
  ! Output for the Reference Sector correction
  ! ------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_rsfields = 5
  CHARACTER (LEN=40), DIMENSION (2,n_rsfields), PARAMETER ::  &
       rs_field_names = RESHAPE ( (/ &
       "ReferenceSectorSlantColumn              ","Reference Sector Slant Column           ", &
       "ReferenceSectorAirMassFactor            ","Reference Sector AMF                    ", &
       "ReferenceSectorApriori                  ","Reference Sector Apriori Vertical Column", &
       "ReferenceSectorFlag                     ","Reference Sector Correction Flag        ", &
       "ColumnAmountReferenceSectorCorrected    ","Column Amount Reference Sector Corrected" /), &
       (/ 2, n_rsfields /) )
  CHARACTER (LEN=14), DIMENSION (3,n_rsfields), PARAMETER ::  &
       rs_field_specs = RESHAPE ( (/ &
       "molec/cm2     ","nXtrack,nTimes","r8            ", &
       "NoUnits       ","nXtrack,nTimes","r8            ", &
       "molec/cm2     ","nXtrack,nTimes","r8            ", &
       "NoUnits       ","nXtrack,nTimes","r8            ", &
       "molec/cm2     ","nXtrack,nTimes","r8            "  /), &
       (/3, n_rsfields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_rsfields ), PARAMETER :: &
       rs_valids = RESHAPE ( (/                              &
       valid_min_r8,         valid_max_r8,                   & ! Reference Sector Slant Columns
       valid_min_r8,         valid_max_r8,                   & ! Reference Sector AMF
       valid_min_r8,         valid_max_r8,                   & ! Reference Sector apriori Vertical Columns
       zero_r8,              1.0_r8,                         & ! Reference Sector Correction Flags
       valid_min_r8,         valid_max_r8 /),                & ! Reference Sector Corrected Vertical Column
       (/2, n_rsfields/) )

  ! --------------------------------------
  ! Output for the de-stripping correction
  ! --------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_defields = 1
  CHARACTER (LEN=24), DIMENSION (2,n_defields), PARAMETER ::  &
       de_field_names = RESHAPE ( (/ &
       "ColumnAmountDestriped   ","Column Amount de striped" /), &
       (/ 2, n_defields /) )
  CHARACTER (LEN=14), DIMENSION (3,n_defields), PARAMETER ::  &
       de_field_specs = RESHAPE ( (/ &
       "molec/cm2     ","nXtrack,nTimes","r8            "  /), &
       (/3, n_defields /) )
  REAL (KIND=r8), DIMENSION ( 2, n_defields ), PARAMETER :: &
       de_valids = RESHAPE ( (/                              &
       valid_min_r8,         valid_max_r8 /),                & ! Column Amount de striped
       (/2, n_defields/) )
    
  ! -------------------------------------------
  ! (2) Additional fields for "diagnostic" runs
  ! -------------------------------------------
  
  ! CCM n_diag_fields 10 -> 18
  INTEGER (KIND=i4), PARAMETER :: n_diag_fields = 19
  
  ! CCM Add All fit residuals
  CHARACTER (LEN=41), DIMENSION (2,n_diag_fields), PARAMETER ::  &
       diagnostic_field_names = RESHAPE ( (/ &
       "CCDPixelRange                            ", "First and Last CCD Pixel Number Fitted   ",    &
       "CommonModeCount                          ", "Common Mode Spectrum Averaging Count     ",    &
       "CommonModeSpectrum                       ", "Common Mode Spectrum                     ",    &
       "CommonModeWavelengths                    ", "Common Mode Spectrum Wavelengths         ",    &
       "CrossTrackStripeCorrection               ", "Correction Factor for Cross Track Stripes",    &
       "FitParameterColumns                      ", "Colum Values of all Fitting Parameters   ",    &
       "FitParameterCorrelations                 ", "Correlations with Main Fitting Parameter ",    &
       "FitSpectrum                              ", "Fitted Spectrum                          ",    &
       "FitMeasuredSpectrum                      ", "Observed L1B Radiance                    ",    &
       "FitMeasuredWavelength                    ", "Spectral Position of L1B Radiance        ",    &
       "FitSpectralFitWeight                     ", "Weighting of L1B Radiance Pixel          ",    &
       "FitParameterNames                        ", "Names of all Fitting Parameters          ",    &
       "FitParameterUncertainty                  ", "Uncertainties of all Fitting Parameters  ",    &
       "DatabaseSpec                             ", "Reference Spectra used in fitting process",    & 
       "DatabaseWavl                             ", "Reference Spectra Wavelengths            ",    &
       "DatabaseNormFactor                       ", "Normalisation factors for DatabaseSpec   ",    &
       "DatabaseSpecNames                        ", "Species Names for DatabaseSpec           ",    &
       "FitIterationCount                        ", "Radiance Fit Iteration Count             ",    &
       "FitResidualSpectrum                      ", "Residual Spectrum                        " /), &
       (/ 2, n_diag_fields /) )
  CHARACTER (LEN=27), DIMENSION (3,n_diag_fields), PARAMETER ::  &
       diagnostic_field_specs = RESHAPE ( (/ &
       "NoUnits                    ","nXtrack,2                  ","i2                         ",   &
       "NoUnits                    ","nXtrack                    ","i4                         ",   &
       "NoUnits                    ","nXtrack,nWaves             ","r8                         ",   &
       "nm                         ","nXtrack,nWaves             ","r4                         ",   &
       "NoUnit                     ","nXtrack,nTimes             ","r8                         ",   &
       "NoUnits                    ","nFitElements,nXtrack,nTimes","r8                         ",   &
       "NoUnits                    ","nFitElements,nXtrack,nTimes","r8                         ",   &
       "NoUnits                    ","nWaves,nXtrack,nTimes      ","r8                         ",   &
       "NoUnits                    ","nWaves,nXtrack,nTimes      ","r8                         ",   &
       "NoUnits                    ","nWaves,nXtrack,nTimes      ","r8                         ",   &
       "NoUnits                    ","nWaves,nXtrack,nTimes      ","r8                         ",   &
       "NoUnits                    ","nCharLenFitElements        ","ch                         ",   &
       "NoUnits                    ","nFitElements,nXtrack,nTimes","r8                         ",   &
       "NoUnits                    ","nRfSpec,nWaves,nXtrack     ","r8                         ",   &
       "nm                         ","nWaves,nXtrack             ","r8                         ",   &
       "NoUnits                    ","nRfSpec                    ","r8                         ",   &
       "NoUnits                    ","nRfSpec                    ","ch                         ",   &
       "NoUnits                    ","nXtrack,nTimes             ","i2                         ",   &
       "NoUnits                    ","nWaves,nXtrack,nTimes      ","r8                         " /),&
       (/3, n_diag_fields /) )
       
  ! Output control logicals CCM
  LOGICAL, DIMENSION ( n_diag_fields ) :: &
       yn_output_diag = (/ .TRUE. ,&  ! CCDPixelRange
       .TRUE. ,&  ! CommonModeCount
       .TRUE. ,&  ! CommonModeSpectrum
       .TRUE. ,&  ! CommonModeWavelengths
       .TRUE. ,&  ! CrossTrackStripeCorrection
       .TRUE. ,&  ! FittingParameterColumns
       .TRUE. ,&  ! FittingParameterCorrelations
       .TRUE. ,&  ! FittedSpectrum
       .TRUE. ,&  ! MeasuredSpectrum
       .TRUE. ,&  ! MeasuredWavelength
       .TRUE. ,&  ! SpectralFitWeight
       .TRUE. ,&  ! FittingParameterNames
       .TRUE. ,&  ! FittingParameterUncertainty
       .TRUE. ,&  ! DatabaseSpec
       .TRUE. ,&  ! DatabaseWavl
       .TRUE. ,&  ! DatabaseNormFactor
       .TRUE. ,&  ! DatabaseSpecNames
       .TRUE. ,&  ! IterationCount
       .TRUE. /)  ! ResidualSpectrum

  ! -----------------------------------------------------------------
  ! Note that the VALIDS for FittingCorrelationElements are not used;
  ! they are only here to make this array of the same dimension as
  ! the others, and so eliminates the need for additional CASE logic.
  ! -----------------------------------------------------------------
  REAL (KIND=r8), DIMENSION ( 2, n_diag_fields ), PARAMETER :: &
       diagnostic_valids = RESHAPE ( (/          &
       1.0_r8,                 780.0_r8,         &  ! CCDPixelRange
       zero_r8,                valid_max_i4,     &  ! CommonModeCount
       valid_min_r8,           valid_max_r8,     &  ! CommonModeSpectrum
       valid_min_r8,           valid_max_r8,     &  ! CommonModeWavelengths
       valid_min_r8,           valid_max_r8,     &  ! CrossTrackStripeCorrection
       valid_min_r8,           valid_max_r8,     &  ! FittingParameterColumns
       valid_min_r8,           valid_max_r8,     &  ! FittingParameterCorrelations
       valid_min_r8,           valid_max_r8,     &  ! FittedSpectrum
       valid_min_r8,           valid_max_r8,     &  ! MeasuredSpectrum
       valid_min_r8,           valid_max_r8,     &  ! MeasuredWavelength
       valid_min_r8,           valid_max_r8,     &  ! SpectralFitWeight
       0.0_r8,                 valid_max_i2,     &  ! FittingParameterNames
       zero_r8,                valid_max_r8,     &  ! FittingParameterUncertainty
       valid_min_r8,           valid_max_r8,     &  ! DatabaseSpec
       valid_min_r8,           valid_max_r8,     &  ! DatabaseWavl
       valid_min_r8,           valid_max_r8,     &  ! DatabaseNormFactor
       0.0_r8,                 valid_max_i2,     &  ! DatabaseSpecNames
       0.0_r8,                 valid_max_i2,     &  ! IterationCount
       valid_min_r8,           valid_max_r8  /), &  ! ResidualSpectrum
       (/ 2, n_diag_fields /) )

  ! ------------------------------------------------------------------
  ! Integer variables for writing to swath data and geolocation fields
  ! ------------------------------------------------------------------
  INTEGER (KIND=C_LONG)                :: he5_start_1d, he5_stride_1d, he5_edge_1d
  INTEGER (KIND=C_LONG), DIMENSION (2) :: he5_start_2d, he5_stride_2d, he5_edge_2d
  INTEGER (KIND=C_LONG), DIMENSION (3) :: he5_start_3d, he5_stride_3d, he5_edge_3d
  INTEGER (KIND=C_LONG), DIMENSION (4) :: he5_start_4d, he5_stride_4d, he5_edge_4d
  INTEGER (KIND=C_LONG), DIMENSION (5) :: he5_start_5d, he5_stride_5d, he5_edge_5d
  INTEGER (KIND=C_LONG), DIMENSION (6) :: he5_start_6d, he5_stride_6d, he5_edge_6d
  INTEGER (KIND=C_LONG), DIMENSION (7) :: he5_start_7d, he5_stride_7d, he5_edge_7d

  ! -----------------------------------
  ! The Vertical Coordinate of the PGEs
  ! -----------------------------------
  ! Since we are doing nadir, this is either "Slant Column" or "Total Column"
  ! -------------------------------------------------------------------------
  CHARACTER (LEN=12), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       vertical_coordinate = (/                                                 &
       "Slant Column", &  ! OClO
       "Total Column", &  ! BrO
       "Total Column", &  ! HCHO
       "Slant Column", &  ! O3
       "Slant Column", &  ! NO2
       "Slant Column", &  ! SO2
       "Total Column", &  ! CHO-CHO
       "Slant Column", &  ! IO
       "Slant Column", &  ! H2O
       "Slant Column", &  ! HONO
       "Slant Column", &  ! O2O2
       "Slant Column", &  ! LqH2O
       "Slant Column"  /) ! NO2(D)
  
  ! -------------------------------
  ! Base names for L2 output swaths
  ! ------------------------------------------------------------------
  ! These Base Name are appended with the PGE target molecules during
  ! the initialization of the HE5 output swath name. They will be
  ! accessed as SWATH_BASE_NAME(pge_idx), and so have to correspond to
  ! the order of PGE indices defined in OMSAO_indices_module:
  !
  !        [1] OClO  [2] BrO  [3] HCHO  [4] O3
  !
  ! The "unofficial" ones are
  !
  !        [5] NO2   [6] SO2  [7] C2H2O2 (Glyoxal)  [8] IO  [9] NO2
  !
  ! Note that these strings will be used IF AND ONLY IF we cannot find
  ! the Swath Name in the PCF file.
  ! ------------------------------------------------------------------
  CHARACTER (LEN=24), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       swath_base_name = (/ &
       "OMPS Slant Column Amount", &  ! OClO    
       "OMPS Slant Column Amount", &  ! BrO
       "OMPS Slant Column Amount", &  ! HCHO
       "OMPS Slant Column Amount", &  ! O3
       "OMPS Slant Column Amount", &  ! NO2
       "OMPS Slant Column Amount", &  ! SO2
       "OMPS Slant Column Amount", &  ! CHO-CHO
       "OMPS Slant Column Amount", &  ! IO
       "OMPS Slant Column Amount", &  ! H2O
       "OMPS Slant Column Amount", &  ! HONO
       "OMPS Slant Column Amount", &  ! O2O2
       "OMPS Slant Column Amount", &  ! LqH2O
       "OMPS Slant Column Amount"  /) ! NO2(D)
  

  ! -------------------------
  ! PGE HDF Global Attributes
  ! ---------------------------------------------------------------
  ! NOTE: These values MIGHT ventually become PSAs, because we want
  !       to search for them in the data set. For now we keep them
  !       as Global Attributes, since this is a much more painless
  !       state than any newly defined PSA.
  ! ---------------------------------------------------------------
  INTEGER (KIND=i4) ::                                &
       NrofScanLines                   = i4_missval,  &
       NrofCrossTrackPixels            = i4_missval,  &
       NrofInputSamples                = i4_missval,  &
       NrofGoodInputSamples            = i4_missval,  &
       NrofGoodOutputSamples           = i4_missval,  &
       NrofMissingSamples              = i4_missval,  &
       NrofSuspectOutputSamples        = i4_missval,  &
       NrofBadOutputSamples            = i4_missval,  &
       NrofConvergedSamples            = i4_missval,  &
       NrofFailedConvergenceSamples    = i4_missval,  &
       NrofExceededIterationsSamples   = i4_missval,  &
       NrofOutofBoundsSamples          = i4_missval
  REAL    (KIND=r4) ::                                &
       PercentGoodOutputSamples        = r4_missval,  &
       PercentSuspectOutputSamples     = r4_missval,  &
       PercentBadOutputSamples         = r4_missval,  &
       AbsolutePercentMissingSamples   = r4_missval

  ! --------------------------------------------
  ! Variables for InputPointer and InputVersions
  ! --------------------------------------------
  INTEGER   (KIND=i4), PARAMETER                 :: n_lun_inp_max = 20
  INTEGER   (KIND=i4), DIMENSION (n_lun_inp_max) :: lun_input
  INTEGER   (KIND=i4)                            :: n_lun_inp
  CHARACTER (LEN=maxchlen)                       :: input_versions

  ! ---------------------------------------------------------------
  ! Finally some system definitions that come with the HE5 Library.
  ! We include this here so that we don't have to worry about it
  ! inside the subroutines that use it.
  ! ---------------------------------------------------------------
  INCLUDE 'hdfeos5.inc'
  INTEGER (KIND = i4), EXTERNAL :: &
       he5_ehrdglatt,  he5_ehwrglatt,  he5_swattach,  he5_swclose,    he5_swcreate,  &
       he5_swdefdfld,  he5_swdefdim,   he5_swdefgfld, he5_swdetach,   he5_swopen,    &
       he5_swrdattr,   he5_swrdfld,    he5_swrdgattr, he5_swrdlattr,  he5_swwrattr,  &
       he5_swwrfld,    he5_swwrgattr,  he5_swwrlattr, HE5_SWfldinfo,                 &
       HE5_SWdefcomp,  HE5_SWdefchunk, HE5_SWsetfill, HE5_SWdefcomch,  HE5_GDOPEN,   &
       HE5_GDATTACH,   HE5_GDRDFLD,    HE5_GDCLOSE,   HE5_GDINQFLDS,   HE5_GDDETACH, &
       HE5_GDINQLATTRS,HE5_GDRDLATTR,  HE5_SWinqdflds,HE5_SWwrcharfld

  INTEGER (KIND=C_LONG), EXTERNAL :: HE5_SWinqswath, HE5_SWinqdims

  ! ---------------------------------------------------------------------------
  ! Parameters for HE5 compression. From the limited experimentation performed,
  ! this combination has shown to produce the smallest file sizes for SAO PGEs.
  ! ---------------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       he5_comp_type   = HE5_HDFE_COMP_SHUF_DEFLATE, &
       he5_comp_par    = 9,                          &
       he5_nocomp_type = HE5_HDFE_COMP_NONE,         &
       he5_nocomp_par  = 0

    ! ------
    ! Groups
    ! ------
    CHARACTER(LEN=12), PARAMETER :: geolocation = '/Geolocation'
    CHARACTER(LEN=5),  PARAMETER :: datafields  = '/Data'
    
    ! ----------
    ! Data field
    ! ----------
    CHARACTER(LEN=17), PARAMETER :: fitres = 'Fitting residuals'
    CHARACTER(LEN=27), PARAMETER :: solres = 'Solar calibration residuals'
    CHARACTER(LEN=29), PARAMETER :: solwav = 'Solar calibration wavelengths'
    CHARACTER(LEN=30), PARAMETER :: radres = 'Radiance calibration residuals'
    CHARACTER(LEN=40), PARAMETER :: rrdres = 'Radiance reference calibration residuals'
    CHARACTER(LEN=25), PARAMETER :: radshi = 'Radiance wavelenght shift'
    CHARACTER(LEN=35), PARAMETER :: rrdshi = 'Radiance reference wavelenght shift'
    CHARACTER(LEN=28), PARAMETER :: solshi = 'Solar wavelenght shift'
    CHARACTER(LEN=30), PARAMETER :: asygau = 'Gaussian asymetry factor'
    CHARACTER(LEN=14), PARAMETER :: widgau = 'Gaussian width'
    CHARACTER(LEN=32), PARAMETER :: databa = 'Database of cross sections'
    CHARACTER(LEN=18), PARAMETER :: solflg = 'Solar fitting flag'
    CHARACTER(LEN=21), PARAMETER :: radflg = 'Radiance fitting flag'
    CHARACTER(LEN=31), PARAMETER :: rrdflg = 'Radiance reference fitting flag'
    CHARACTER(LEN=12), PARAMETER :: fitflg = 'Fitting flag'
    CHARACTER(LEN=14), PARAMETER :: varfit = 'Fitting values'
    CHARACTER(LEN=16), PARAMETER :: solirr = 'Solar Irradiance'
    CHARACTER(LEN=18), PARAMETER :: radref = 'Radiance Reference'
    CHARACTER(LEN= 8), PARAMETER :: radian = 'Radiance'
    CHARACTER(LEN=19), PARAMETER :: piqflg = 'Pixel Quality Flags'
    

    ! ----------------
    ! Geolocation fields
    ! ------------------
    CHARACTER(LEN=8),  PARAMETER :: latfld  = 'Latitude'
    CHARACTER(LEN=9),  PARAMETER :: lonfld  = 'Longitude'
    CHARACTER(LEN=17), PARAMETER :: saafld  = 'SolarAzimuthAngle'
    CHARACTER(LEN=16), PARAMETER :: szafld  = 'SolarZenithAngle'
    CHARACTER(LEN=18), PARAMETER :: scafld  = 'SpacecraftAltitude'
    CHARACTER(LEN=15), PARAMETER :: trpfld  = 'TerrainPressure'
    CHARACTER(LEN=19), PARAMETER :: vaafld  = 'ViewingAzimuthAngle'
    CHARACTER(LEN=18), PARAMETER :: vzafld  = 'ViewingZenithAngle'
    CHARACTER(LEN=23), PARAMETER :: gpqffld = 'GroundPixelQualityFlags'
    CHARACTER(LEN=22), PARAMETER :: iqffld  = 'InstrumentQualityFlags'
    CHARACTER(LEN=16), PARAMETER :: sedfld  = 'SunEarthDistance'

    INTEGER (HID_T)             :: file_id, group_id, dspace_id, dset_id, list_id, type_id, &
         attr_id, atype_id, aspace_id

    ! -----------------
    ! Define dimensions
    ! -----------------
    INTEGER(HSIZE_T), DIMENSION(1) :: he5_1D, adims, start_1d, count_1d, stride_1d, dims_1d, block_1d
    INTEGER(HSIZE_T), DIMENSION(2) :: he5_2D, start_2d, count_2d, stride_2d, dims_2d, block_2d
    INTEGER(HSIZE_T), DIMENSION(3) :: he5_3D, start_3d, count_3d, stride_3d, dims_3d, block_3d

CONTAINS

  SUBROUTINE he5_init_input_file ( &
       file_name, swath_name, swath_id, swath_file_id, ntimes_aux, nxtrack_aux, he5stat )

    !------------------------------------------------------------------------------
    ! This subroutine initializes the HE5 input files for reading
    ! prefitted BrO and O3 that are used in OMHCHO
    !
    ! Input:
    !   file_name         - Name of HE5 input file
    !
    ! Output:
    !   swath_name.............Name of existing swath in file
    !   swath_file_id_inp......id number for HE5 input file (required for closing it)
    !   swath_id_inp...........id number for swath (required for reading from swath)
    !   ntimes_aux.............nTimes  as given in product file
    !   nxtrack_aux............nXtrack as given in product file
    !   he5stat................OMI_E_SUCCESS if everything went well
    !
    ! No Swath ID Variables passed through MODULE.
    !
    !------------------------------------------------------------------------------

    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=19), PARAMETER :: modulename = 'he5_init_input_file'

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: file_name

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER   (KIND=i4),      INTENT (INOUT) :: he5stat
    INTEGER   (KIND=i4),      INTENT (OUT)   :: swath_id, swath_file_id, nxtrack_aux, ntimes_aux
    CHARACTER (LEN=maxchlen), INTENT (OUT)   :: swath_name

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                       :: nsep, ndim
    INTEGER   (KIND=i4)                       :: i, errstat, swlen, iend, istart
    INTEGER   (KIND=C_LONG)                   :: ndimcl, errstatcl
    INTEGER   (KIND=i4),      DIMENSION(0:12) :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG),  DIMENSION(0:12) :: dim_arraycl
    CHARACTER (LEN=maxchlen)                  :: dim_chars

    errstat = pge_errstat_ok

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdonly )
    IF ( swath_file_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, he5stat )
       IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    errstatcl = HE5_SWinqswath  ( TRIM(ADJUSTL(file_name)), swath_name, swlen )
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
       IF ( he5stat >= pge_errstat_error ) RETURN
    END IF

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:12) )
    ndim   = INT ( ndimcl, KIND=i4 )
    IF ( ndim <= 0 ) THEN
       he5stat = MAX ( he5stat, pge_errstat_error )
       RETURN
    END IF
    dim_array(0:12) = INT ( dim_arraycl(0:12), KIND=i4 )

    ! -----------------------------------------
    ! Extract nTimes and nXtrack from the swath
    ! -----------------------------------------
    nxtrack_aux = 0  ;  ntimes_aux = 0
    dim_chars = TRIM(ADJUSTL(dim_chars))
    swlen = LEN_TRIM(ADJUSTL(dim_chars))
    istart = 1  ;  iend = 1

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 ; dim_seps = 0 ; dim_seps(0) = 0
    getseps: DO i = 1, swlen 
       IF ( dim_chars(i:i) == ',' ) THEN
          nsep = nsep + 1
          dim_seps(nsep) = i
       END IF
    END DO getseps
    nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO i = 0, nsep-1
       istart = dim_seps(i)+1 ; iend = dim_seps(i+1)-1
       IF  ( dim_chars(istart:iend) == "nTimes"  ) ntimes_aux  = dim_array(i)
       IF  ( dim_chars(istart:iend) == "nXtrack" ) nxtrack_aux = dim_array(i)
       IF ( ntimes_aux > 0 .AND. nxtrack_aux > 0 ) EXIT getdims
    END DO getdims

    RETURN
  END SUBROUTINE he5_init_input_file

  SUBROUTINE write_solar_calibration(ip, nw, fit_res_out, flag)

    USE OMSAO_variables_module,  ONLY: fitvar_cal
    USE OMSAO_data_module
    USE OMSAO_indices_module,    ONLY: shi_idx

  
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=4), INTENT(IN)      :: ip,nw
    REAL    (KIND=8), DIMENSION(1:nw), INTENT(IN) :: fit_res_out
    INTEGER (KIND=4), INTENT(IN)      :: flag
    ! --------------
    ! Error variable
    ! --------------
    INTEGER (KIND=4) :: errstat
    
    CHARACTER(LEN=8), PARAMETER :: file_test='test.he5'

    INTEGER(HID_T) :: fid, mid
    INTEGER        :: rank
    REAL(KIND=8), DIMENSION(nw,1) :: data_2d

    ! Initialize HDF5 interface
    CALL h5open_f(errstat)

    ! Open output file
    CALL h5fopen_f(file_test, H5F_ACC_RDWR_F, file_id, errstat)
    ! Open data group
    CALL h5gopen_f(file_id, datafields, group_id, errstat)

    ! Open Solar wavelenght shift
    CALL h5dopen_f(group_id, solshi, dset_id, errstat)
    ! Get dataset's dataspace identifier
    CALL h5dget_space_f(dset_id, fid, errstat)
    ! Select file hyper slab
    start_1d(1)  = ip-1
    count_1d(1)  = 1
    stride_1d(1) = 1
    dims_1d(1)   = 1
    rank         = 1
    CALL h5sselect_hyperslab_f(fid, H5S_SELECT_SET_F, start_1d, count_1d, errstat, stride_1d)
    CALL h5screate_simple_f(rank,dims_1d,mid,errstat)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REAL(fitvar_cal(shi_idx), KIND=8), &
         dims_1d, errstat, mid, fid)
    ! Close solar wavelenght shift
    CALL h5sclose_f(fid, errstat)
    CALL h5sclose_f(mid, errstat)
    CALL h5dclose_f(dset_id, errstat)

    ! Open Solar fitting quality flag
    CALL h5dopen_f(group_id, solflg, dset_id, errstat)
    ! Get dataset's dataspace identifier
    CALL h5dget_space_f(dset_id, fid, errstat)
    ! Select file hyper slab
    start_1d(1)  = ip-1
    count_1d(1)  = 1
    stride_1d(1) = 1
    dims_1d(1)   = 1
    rank         = 1
    CALL h5sselect_hyperslab_f(fid, H5S_SELECT_SET_F, start_1d, count_1d, errstat, stride_1d)
    CALL h5screate_simple_f(rank,dims_1d,mid,errstat)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, flag, dims_1d, errstat, mid, fid)
    ! Close solar flag
    CALL h5sclose_f(fid, errstat)
    CALL h5sclose_f(mid, errstat)
    CALL h5dclose_f(dset_id, errstat)

    ! Open Solar residuals
    CALL h5dopen_f(group_id, solres, dset_id, errstat)
    ! Get dataset's dataspace identifier
    CALL h5dget_space_f(dset_id, fid, errstat)
    ! Create memory dataspace
    rank         = 2
    dims_2d(1)   = nw
    dims_2d(2)   = 1
    CALL h5screate_simple_f(rank,dims_2d,mid,errstat)
    ! Select file hyperslab
    start_2d(1)  = 0
    start_2d(2)  = ip-1
    count_2d(1)  = nw
    count_2d(2)  = 1
    stride_2d(1) = 1
    stride_2d(2) = 1
    block_2d(1)  = 1
    block_2d(2)  = 1
    dims_1d(1)   = nw
    data_2d(1:nw,1) = fit_res_out(1:nw)
    CALL h5sselect_hyperslab_f(fid, H5S_SELECT_SET_F, start_2d, count_2d, errstat)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REAL(data_2d, KIND=8), &
         dims_2d, errstat, mid, fid)
    ! Close solar residual
    CALL h5sclose_f(fid, errstat)
    CALL h5sclose_f(mid, errstat)
    CALL h5dclose_f(dset_id, errstat)
    


    ! Close data group
    CALL h5gclose_f(group_id, errstat)
    ! Close output file
    CALL h5fclose_f(file_id, errstat)
    
    

  END SUBROUTINE write_solar_calibration
END MODULE OMSAO_he5_module
