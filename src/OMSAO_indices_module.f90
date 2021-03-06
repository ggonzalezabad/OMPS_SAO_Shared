MODULE OMSAO_indices_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  IMPLICIT NONE

  ! ===============================================================
  ! In this MODULE we collect all indices that we have defined to
  ! generalize the match between reference spectra and their
  ! associated fitting parameters. Some of the indices overlap with
  ! those required for unique idenfication of I/O files for the OMI
  ! PGE PCF file.
  ! ===============================================================

  ! ------------------------------------------------------------------
  ! The following list represents all the reference spectra we can
  ! encounter during OMPS fitting of OClO, BrO, and HCHO (the "PGEs").
  ! Not all entries are required for all three molecules but each of 
  ! them is present in at least one PGE. 
  !
  ! The entry for the fitting input control file (ICF) is not a 
  ! reference spectrum, and is therefore assigned entry "0".
  !
  ! Each index (except ICF) has an identification string  
  ! associated with it, so that it can be used for the identification
  ! of reference spectra and fitting parameters in the ICF.
  !
  !     icf:    Input control file (not a reference spectrum)
  !
  !     solar:  Solar reference (usually Kitt Peak)
  !     ring:   Ring
  !     o3_t1:  O3, first  temperature
  !     o3_t2:  O3, second temperature
  !     o3_t3:  O3, third  temperature
  !     no2_t1: NO2, first  temperature
  !     no2_t2: NO2, second temperature
  !     o2o2:   O2-O2 collision complex (usually Greenblatt)
  !     so2:    SO2
  !     bro:    BrO
  !     oclo:   OClO
  !     hcho:   HCHO
  !     h2o:    H2O
  !     lqh2o:  Liquid Water Absorption
  !     glyox:  Glyoxal (CHOCHO)
  !     io:     Iodine Monoxide (IO)
  !     hono:   Nitrous Acid (HONO)
  !     vraman: Vibrational Raman ("Water Ring")
  !     comm:   Common mode
  !     resid:  Residual (like a common mode, or vice versa)
  !     pseudo1:Pseudo Absorber (no real use identified yet)
  !     pseudo2:Pseudo Absorber (no real use identified yet)
  !     pseudo3:Pseudo Absorber (no real use identified yet)
  !     polcor: Polarization correction
  !     us1:    First undersampling spectrum
  !     us2:    Second undersampling spectrum
  !     noname: Not yet determined, dummy placeholder
  ! ---------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       icf_idx    =  0, solar_idx  =  1, ring_idx   =  2, o3_t1_idx  =  3, &
       o3_t2_idx  =  4, o3_t3_idx  =  5, no2_t1_idx =  6, no2_t2_idx =  7, &
       o2o2_idx   =  8, so2_idx    =  9, bro_idx    = 10, oclo_idx   = 11, &
       hcho_idx   = 12, h2o_idx    = 13, lqh2o_idx  = 14, glyox_idx  = 15, &
       io_idx     = 16, hono_idx   = 17, vraman_idx = 18, comm_idx   = 19, &
       resid_idx  = 20, pabs1_idx  = 21, pabs2_idx  = 22, pabs3_idx  = 23, &
       polcor_idx = 24, us1_idx    = 25, us2_idx    = 26, noname_idx = 27

  ! ----------------------------------------------------------
  ! The minimum and maximum indices of reference spectra (rs).
  ! ----------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: min_rs_idx = solar_idx, max_rs_idx = noname_idx

  ! --------------------------------------------------------------
  ! The identification strings associated with the fitting indices
  ! --------------------------------------------------------------
  CHARACTER (LEN=7), DIMENSION (min_rs_idx:max_rs_idx), PARAMETER ::         &
       refspec_strings = (/                                                  &
       'solar  ', 'ring   ', 'o3_t1  ', 'o3_t2  ', 'o3_t3  ', 'no2_t1 ', 'no2_t2 ', &
       'o2o2   ', 'so2    ', 'bro    ', 'oclo   ', 'hcho   ', 'h2o    ', 'lqh2o  ', &
       'glyox  ', 'io     ', 'hono   ', 'vraman ', 'commod ', 'resid  ', 'pseudo1', &
       'pseudo2', 'pseudo3', 'polcor ', 'usamp1 ', 'usamp2 ', 'noname '            /)

  CHARACTER (LEN=24), DIMENSION (min_rs_idx:max_rs_idx), PARAMETER :: &
       refspec_titles = (/ &
       'Solar Reference Spectrum', &
       'Molecular Ring          ', &
       'O3 at Temperature T1    ', &
       'O3 at Temperature T2    ', &
       'O3 at Temperature T3    ', &
       'NO2 at Temperature T1   ', &
       'NO2 at Temperature T2   ', &
       'O2-O2                   ', &
       'SO2                     ', &
       'BrO                     ', &
       'OClO                    ', &
       'HCHO                    ', &
       'H2O                     ', &
       'Liquid H2O              ', &
       'CHOCHO (Glyoxal)        ', &
       'IO                      ', &
       'HONO                    ', &
       'Vibrational Raman       ', &
       'Common Mode Spectrum    ', &
       'Fitting Residuals       ', &
       'Pseudo Absorber 1       ', &
       'Pseudo Absorber 2       ', &
       'Pseudo Absorber 3       ', &
       'Polarization Correction ', &
       'Undersampling Phase 1   ', &
       'Undersampling Phase 2   ', &
       'NaN                     '  /)

  ! ==============================================
  ! Now we define the specific fitting parameters.
  ! ==============================================
  ! -----------------------------------------------------------------
  ! Particular fitting parameters: Solar fit and radiance calibration
  ! -----------------------------------------------------------------
  !     bl0: baseline, 0 order
  !     bl1: baseline, 1 order
  !     bl2: baseline, 2 order
  !     bl3: baseline, 3 order
  !     bl4: baseline, 4 order
  !     bl5: baseline, 5 order
  !     sc0: scaling,  0 order
  !     sc1: scaling,  1 order
  !     sc2: scaling,  2 order
  !     sc3: scaling,  3 order
  !     sc4: scaling,  4 order
  !     sc5: scaling,  5 order
  !     sin: solar intensity
  !     hwe: slit width at 1/e
  !     asy: slit function asymmetry
  !     shi: spectral shift
  !     squ: spectral squeeze
  ! --------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       bl0_idx =  1, bl1_idx =  2, bl2_idx =  3, bl3_idx =  4, bl4_idx =  5, bl5_idx =  6, &
       sc0_idx =  7, sc1_idx =  8, sc2_idx =  9, sc3_idx = 10, sc4_idx = 11, sc5_idx = 12, &
       sin_idx = 13, hwe_idx = 14, asy_idx = 15, sha_idx = 16, shi_idx = 17, squ_idx = 18, &
       max_calfit_idx = squ_idx

  CHARACTER (LEN=3), DIMENSION (max_calfit_idx), PARAMETER :: calfit_strings = (/ &
       'bl0', 'bl1', 'bl2', 'bl3', 'bl4', 'bl5', 'sc0', 'sc1', 'sc2', 'sc3', 'sc4', &
       'sc5', 'sin', 'hwe', 'asy', 'sha', 'shi', 'squ' /)

  CHARACTER (LEN=29), DIMENSION (max_calfit_idx), PARAMETER :: calfit_titles = (/ &
       'Baseline Polynomial 0th Order', &
       'Baseline Polynomial 1st Order', &
       'Baseline Polynomial 2nd Order', &
       'Baseline Polynomial 3rd Order', &
       'Baseline Polynomial 4th Order', &
       'Baseline Polynomial 5th Order', &
       'Scaling Polynomial 0th Order ', &
       'Scaling Polynomial 1st Order ', &
       'Scaling Polynomial 2nd Order ', &
       'Scaling Polynomial 3rd Order ', &
       'Scaling Polynomial 4th Order ', &
       'Scaling Polynomial 5th Order ', &
       'Solar Intensity Contribution ', &
       'Slit Function Half Width @1/e', &
       'Slit Function Asymmetry      ', &
       'Slit Function Shape Parameter', &
       'Wavelength Shift             ', &
       'Wavelength Squeeze           ' /)

  ! --------------------------------------------------------
  ! Particular fitting parameters: Radiance spectral fitting
  ! --------------------------------------------------------
  !     ad1: first added contribution
  !     lbe: Lambert-Beer terms
  !     ad2: second added contribution
  !     mns: minimum radiance fitting sub-index
  !     mxs: maximum radiance fitting sub-index
  ! -------------------------------------------
  INTEGER   (KIND=i4),                    PARAMETER :: ad1_idx = 1, lbe_idx = 2, ad2_idx = 3
  INTEGER   (KIND=i4),                    PARAMETER :: mns_idx = ad1_idx, mxs_idx = ad2_idx
  CHARACTER (LEN=3), DIMENSION (mxs_idx), PARAMETER :: &
       radfit_strings = (/ 'ad1', 'lbe', 'ad2' /)
  CHARACTER (LEN=14), DIMENSION (mxs_idx), PARAMETER :: &
       radfit_titles = (/ '(added first) ', '(Lambert-Beer)', '(added last)  ' /)

  ! -----------------------------------------------------------------
  ! Total number of (radiance) fitting parameters: The MAX_CALFIT_IDX
  ! calibration parameters, plus MAX_RS_IDX*AD2_IDX parameters
  ! associated with external or online computed reference spectra.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_max_fitpars = max_calfit_idx + mxs_idx*max_rs_idx

  ! --------------------------------
  ! ELSUNC fitting constraints index
  ! --------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       elsunc_unconstrained = 0, elsunc_same_lower = 1, elsunc_userdef = 2

  ! --------------------------------------------------------------------
  ! Indices for fitting wavelengths, spectrum, weights, and CCD detector
  ! --------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: wvl_idx = 1, spc_idx = 2, sig_idx = 3, ccd_idx = 4

  ! ---------------------------------------------------------------
  ! Indices for Solar and Radiance Wavelength Calibration, Radiance 
  ! Reference Fit, and reguar Radiance Fits
  ! ---------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       solcal_idx = 1, radcal_idx = 2, radref_idx = 3, radfit_idx = 4

  ! ------------------------------------------------------------------
  ! Quality flag indices. Used for flagging any CCD detector pixel for
  ! exclusion from the fitting process. For details on the definition
  ! of these flags see the Dutch Space L1B produce document. <--FIXME Need to be adapted to OMPS description
  ! ------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       qflg_mis_idx = 0, qflg_bad_idx = 1, qflg_err_idx = 2, qflg_tra_idx = 3, &
       qflg_rts_idx = 4, qflg_sat_idx = 5, qflg_noi_idx = 6, qflg_drk_idx = 7, &
       qflg_off_idx = 8, qflg_exp_idx = 9, qflg_str_idx = 10
  
  ! =================================================================
  ! This module defines a range of indices for the OMI SAO PGE fitting
  ! processes. These are mostly associated with file unit numbers as
  ! they are required by the PGE Process Control File (PCF). But they
  ! also include the number of SAO PGEs, the official PGE numbers for
  ! the molecules, etc.
  ! =================================================================

  ! ------------------------------------------------------------------------
  ! First define the number of SAO PGEs. We are fitting three molecules, so:
  ! ------------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_sao_pge = 4
  
  ! -----------------------------------------------------------------
  ! The operational OMI environment assigns reference numbers to each
  ! PGE. Here we adopt the official setting for future indexing.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       pge_oclo_idx = 10, pge_bro_idx  = 11, pge_hcho_idx  = 12, pge_o3_idx    = 13, &
       pge_no2_idx  = 14, pge_so2_idx  = 15, pge_gly_idx   = 16, pge_io_idx    = 17, &
       pge_h2o_idx  = 18, pge_hono_idx = 19, pge_o2o2_idx  = 20, pge_lqh2o_idx = 21, &
       pge_no2d_idx = 22
  INTEGER (KIND=i4), PARAMETER :: sao_pge_min_idx = pge_oclo_idx, sao_pge_max_idx = pge_no2d_idx

  ! ------------------------------------------------
  ! Identifications strings and indices for SAO PGEs
  ! ------------------------------------------------
  CHARACTER (LEN=8), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       sao_pge_names      = (/                                                 &
       'OMOCLO  ', 'OMBRO   ', 'OMHCHO  ', 'OMSAO3  ', 'OMSAONO2',             &
       'OMSAOSO2', 'OMCHOCHO', 'OMSAOIO ', 'OMH2O   ', 'OMHONO  ',             &
       'OMO2O2  ', 'OMLQH2O ', 'OMNO2D  '                                      /)
  CHARACTER (LEN=6), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       sao_molecule_names = (/                                                 &
       'OClO  ',   'BrO   ',   'HCHO  ',   'O3    ', 'NO2   ',                 &
       'SO2   ',   'CHOCHO',   'IO    ',   'H2O   ', 'HONO  ',                 &
       'O2O2  ',   'LqH2O ',   'NO2   '                                        /)

  ! ----------------------------------------------------------------------
  ! What follows are the input file unit numbers required for uniquely 
  ! identifying files in the PCF. Remember that each PGE has a separate
  ! OMI PGE number, and thus we have to define three plus one sets of 
  ! file units. The fourth set is that for O3, which requires a full PGE
  ! set of inputs without having the benefit of being an official SAO PGE.
  ! The easiest way to do this is to assign two full sets of input numbers
  ! for the HCHO PGE (#12).
  !
  ! Note on RESHAPE: 
  !   This array functions reshapes a 1-dim array into a multi-dimensional
  !   one. The first (new) array dimension varies fastest, followed by the
  !   second, and so on.
  ! ----------------------------------------------------------------------

  ! -------------------------------------------------------------
  ! Logical Unit Numbers (LUNs) in PCF. Some registered, some not.
  ! -------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       granule_s_lun           =  10258, & ! LUN for GranuleStartTime
       granule_e_lun           =  10259, & ! LUN for GranuleStartTime
       verbosity_lun           = 200100, & ! LUN for PGE verbosity threshold
       version_lun             = 200105, & ! LUN for PGEVERSION
       proccenter_lun          = 200110, & ! LUN for ProcessingCenter
       prochost_lun            = 200115, & ! LUN for ProcessingHost
       reproc_actual_lun       = 200135, & ! LUN for ReprocessingActual
       proclevel_lun           = 200170, & ! LUN for ProcessLevel
       instrument_name_lun     = 200175, & ! LUN for InstrumentName
       opmode_lun              = 200180, & ! LUN for OperationMode
       authoraffil_lun         = 200185, & ! LUN for Author Affiliation
       authorname_lun          = 200190, & ! LUN for Author Name
       orbitnumber_lun         = 200200, & ! LUN for orbit number
       swathname_lun           = 200210, & ! LUN for swath name
       versionid_lun           = 200220, & ! LUN for VersionID (formerly in MCF)
       molid_lun               = 700000, & ! LUN for PGE molecule ID; registered
       mcf_lun                 = 700001, & ! LUN for MCF file
       l1b_radiance_lun        = 700010, & ! LUN L1B radiance file (both UV or VIS)
       l1b_radianceref_lun     = 700015, & ! LUN L1B Radiance Reference file
       l1b_irradiance_lun      = 700020, & ! LUN L1B Solar irradiance  file
       slitfunc_lun            = 700050, & ! LUN for slit function data tile
       amf_table_lun           = 700210, & ! LUN AMF LUT file
       cld_climatology_lun     = 700230, & ! LUN cloud climatology file
       cld_lun                 = 700240, & ! LUN L2 cloud file
       to3_lun                 = 700250, & ! LUN L2 Total Ozone file
       albedo_lun              = 700280, & ! LUN for albedo file
       prefit_lun              = 700301, & ! LUN for prefit columns file
       solcomp_lun             = 700400, & ! LUN for solar composite spectrum file
       solmonthave_lun         = 700500, & ! LUN for solar monthly mean spectrum file
       refsec_lun              = 700600, & ! LUN for reference sector file
       refsec_cld_lun          = 700615, & ! LUN for L2 reference sector cloud file
       l2_output_lun           = 700999, & ! LUN for L2 output file
       format_lun              = 200221, & ! LUN for L2 file format
       localgranuleID_lun      = 200222, & ! LUN for local granule ID
       longname_lun            = 200223, & ! LUN for long name
       pgeversion_lun          = 200224, & ! LUN for pgeversion
       shortname_lun           = 200225, & ! LUN for short name
       creator_email_lun       = 200226, & ! LUN for creator email
       creator_institution_lun = 200227, & ! LUN for creator institution
       creator_name            = 200228, & ! LUN for creator name
       keywords_lun            = 200229, & ! LUN for keywords
       processing_level_lun    = 200230, & ! LUN for processing level
       title_lun               = 200231, & ! LUN for title
       locality_values_lun     = 200232, & ! LUN for locality values
       molecule_id_lun         = 200233    ! LUN for molecule_id


  ! -----------------------------------------------------------------
  ! Place any of the above LUNs that are associated with Config data
  ! into a single array, which is then used for reading all LUNs in 
  ! one loop. The CHARACTER arrays below hold the LUN explanation 
  ! (mainly for error handling), and the values that are read from 
  ! the PCF.
  ! 
  ! "Config" data are read using the PGS_PC_GetConfigData function.
  ! Note that this _excludes_ any file names, which must be read
  ! with the PGS_PC_GetReference function.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_config_luns = 29
  INTEGER (KIND=i4), DIMENSION (n_config_luns), PARAMETER :: config_lun_array = (/  &
       granule_s_lun, granule_e_lun, verbosity_lun, version_lun,&
       proccenter_lun, prochost_lun, reproc_actual_lun, proclevel_lun,&
       instrument_name_lun, opmode_lun, authoraffil_lun, authorname_lun,&
       orbitnumber_lun, swathname_lun, versionid_lun, molid_lun, &
       format_lun, localgranuleID_lun, longname_lun, pgeversion_lun, &
       shortname_lun, creator_email_lun, creator_institution_lun, creator_name, &
       keywords_lun, processing_level_lun, title_lun, locality_values_lun, &
       molecule_id_lun /)

  CHARACTER (LEN=29), DIMENSION (n_config_luns), PARAMETER :: config_lun_strings = (/  &
       "Granule Start Time           ", "Granule End Time             ", &
       "PGE verbosity threshold      ", "PGEVERSION                   ", &
       "ProcessingCenter             ", "ProcessingHost               ", &
       "REPROCESSINGACTUAL           ", "ProcessLevel                 ", &
       "InstrumentName               ", "OPERATIONMODE                ", &
       "AuthorAffiliation            ", "AuthorName                   ", &
       "OrbitNumber                  ", "SwathName                    ", &
       "VERSIONID                    ", "PGE Molecule ID              ", &
       "Format                       ", "LocalGranuleID               ", &
       "LongName                     ", "PGEVersion                   ", &
       "ShortName                    ", "creator_email                ", &
       "creator_institution          ", "creator_name                 ", &
       "keywords                     ", "processing_level             ", &
       "title                        ", "locality_values              ", &
       "molecule_id                  " /)

  CHARACTER (LEN=3), DIMENSION (n_config_luns), PARAMETER :: config_lun_autocopy = (/  &
       "no ", "no ", "no ", "yes", &
       "yes", "yes", "no ", "yes", &
       "yes", "yes", "yes", "yes", &
       "no ", "no ", "no ", "no ", &
       "no ", "no ", "no ", "no ", &
       "no ", "no ", "no ", "no ", &
       "no ", "no ", "no ", "no ", &
       "no " /)

  LOGICAL, DIMENSION (n_config_luns), PARAMETER :: yn_config_lun_autocopy = (/  &
       .FALSE., .FALSE., .FALSE., .TRUE.,  &
       .TRUE.,  .TRUE.,  .FALSE., .TRUE.,  &
       .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  &
       .FALSE., .FALSE., .FALSE., .FALSE., &
       .FALSE., .FALSE., .FALSE., .FALSE., &
       .FALSE., .FALSE., .FALSE., .FALSE., &
       .FALSE., .FALSE., .FALSE., .FALSE., &
       .FALSE. /)

  CHARACTER (LEN=maxchlen), DIMENSION (n_config_luns) :: config_lun_values
  ! ----------------------------------------------------------------------
  ! * Input file LUNs for reference spectra and algorithm control file.
  !   The first (ICF) is for static, non-reference spectra input files.
  ! ----------------------------------------------------------------------
  ! 700100: Input Control File
  ! 700101: Solar  spectrum
  ! 700102: Ring   spectrum
  ! 700103: O3     absorption cross sections at temperature T1
  ! 700104: O3     absorption cross sections at temperature T2
  ! 700105: O3     absorption cross sections at temperature T3
  ! 700106: NO2    absorption cross sections at temperature T1
  ! 700107: NO2    absorption cross sections at temperature T2
  ! 700108: O2-O2  absorption cross sections
  ! 700109: SO2    absorption cross sections
  ! 700110: BrO    absorption cross sections
  ! 700111: OClO   absorption cross sections
  ! 700112: HCHO   absorption cross sections
  ! 700113: H2O    absorption cross sections
  ! 700114: Liquid H2O absorption
  ! 700115: CHOCHO absorption cross sections                         
  ! 700116: IO     absorption cross sections                         
  ! 700117: HONO   absorption cross sections                         
  ! 700118: Vibrational Raman Scattering ("Water Ring")              
  ! 700119: Common mode spectrum                                     
  ! 700120: Fitting Residuals                                        
  ! 700121: Pseudo absorber (not yet used or implemented)            
  ! 700122: Polarization correction (not yet used or implemented)    
  ! 700123: Undersampling spectrum, Phase 1 (computed on-line)       
  ! 700124: Undersampling spectrum, Phase 2 (computed on-line)       
  ! 700125: BrO slant columns (for OMHCHO, pre-fitted from OMBRO)    
  ! 700126: O3  slant columns (for OMHCHO, pre-fitted from OMSAO3)   
  ! 700127: NoName dummy index (defines maximum index number)        
  ! ----------------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (icf_idx:max_rs_idx), PARAMETER :: &
       pge_static_input_luns = (/ &
       700100,  700101,  700102,  700103,  700104,  700105,  700106,  700107, &
       700108,  700109,  700110,  700111,  700112,  700113,  700114,  700115, &
       700116,  700117,  700118,  700119,  700120,  700121,  700122,  700123, &
       700124,  700125,  700126,  700127                                      /)

  ! -------------------------------------------------------
  ! Strings to search for in the fitting control input file
  ! -------------------------------------------------------
  TYPE ctr_strings
  ! ---------------------------------
  ! Processing mode selection strings
  ! ---------------------------------
     CHARACTER (LEN=10) :: procmode_diag = 'diagnostic'
     CHARACTER (LEN=21) :: comline_str     = 'Common mode iteration'
     CHARACTER (LEN=22) :: destriping_str  = 'Cross-Track Smoothing'
     CHARACTER (LEN=28) :: fitresconst_str = 'Fitting Residual Constraints'
     CHARACTER (LEN=26) :: genline_str     = 'General fitting parameters'
     CHARACTER (LEN=24) :: iofline_str     = 'Input/output data files'
     CHARACTER (LEN=18) :: molline_str     = 'Molecule(s) to fit'
     CHARACTER (LEN=26) :: maxgoodcol_str  = 'Maximum Good Column Amount'
     CHARACTER (LEN=22) :: nrmline_str     = 'Spectrum normalization'
     CHARACTER (LEN=37) :: o3amf_str       = 'Correct O3 AMF Wavelength Dependence?'
     CHARACTER (LEN=37) :: wfmod_amf_str   = 'Use weighting-function modified AMFs?'
     CHARACTER (LEN=15) :: procline_str    = 'Processing mode'
     CHARACTER (LEN=31) :: racline_str     = 'Radiance calibration parameters'
     CHARACTER (LEN=27) :: rafline_str     = 'Radiance fitting parameters'
     CHARACTER (LEN=26) :: rrsline_str     = 'Radiance reference setting'
     CHARACTER (LEN=24) :: rspline_str     = 'Input reference spectra'
     CHARACTER (LEN=24) :: scpline_str     = 'Solar composite spectrum'
     CHARACTER (LEN=28) :: socline_str     = 'Solar calibration parameters'
     CHARACTER (LEN=23) :: wavwindow_str   = 'Fitting window settings'
     CHARACTER (LEN=21) :: solmonthave_str = 'Solar monthly average'
     CHARACTER (LEN=27) :: refseccor_str   = 'Reference sector correction'
     CHARACTER (LEN=20) :: insowavcal_str  = 'Solar inverse wavcal'
     CHARACTER (LEN= 3) :: eoi3str = 'eoi'
  END type ctr_strings
  TYPE(ctr_strings) :: ctrstr

  ! -------------------
  ! End of input string
  ! -------------------
  CHARACTER (LEN=12) :: eoi_str = 'end_of_input'

  ! ----------------
  ! Metadata indices
  ! ----------------
  INTEGER (KIND=i4), PARAMETER :: md_inventory_idx = 2, md_archive_idx = 3

  ! -----------------------------------------
  ! Indices for diagnostic output control CCM
  ! -----------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       ccdpix_didx    =  1,& ! CCD Pixel Range
       commcnt_didx   =  2,& ! Common Mode Count
       commspc_didx   =  3,& ! Common Mode Spectrum
       commwvl_didx   =  4,& ! Common Mode Wavelengths
       xtrcor_didx    =  5,& ! Cross Track Stripe Correction
       corrcol_didx   =  6,& ! Fitting Parameter Columns
       corr_didx      =  7,& ! Fitting Parameter Correlations
       spcfit_didx    =  8,& ! Fitted Spectrum
       spcobs_didx    =  9,& ! Measured Spectrum 
       posobs_didx    = 10,& ! Measured Wavelength
       fitwt_didx     = 11,& ! Spectral Fit Weight
       correlm_didx   = 12,& ! Fitting Parameter Names
       correrr_didx   = 13,& ! Fitting Parameter Uncertainty
       spdata_didx    = 14,& ! Database Spec
       spdatw_didx    = 15,& ! Database Spec Wavel
       spnrmf_didx    = 16,& ! Database Norm Factor
       spname_didx    = 17,& ! Database Spec Names
       itnum_didx     = 18,& ! Iteration Count
       spcres_didx    = 19   ! Residual Spectrum

END MODULE OMSAO_indices_module

