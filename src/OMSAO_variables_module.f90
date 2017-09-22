MODULE OMSAO_variables_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
       o3_t1_idx, o3_t2_idx, o3_t3_idx,                                      &
       us1_idx, us2_idx, n_voc_amf_luns, ccd_idx, radfit_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, max_spec_pts, n_fit_winwav, max_mol_fit
  USE OMSAO_omidata_module, ONLY: nwavel_max, nxtrack_max, nlines_max

  IMPLICIT NONE

  ! -----------------------------------------------------
  ! Variables read from PCF file
  ! -----------------------------------------------------
  ! * Current PGE name and index
  ! * Verbosity threshold
  ! * Orbit number, Version ID, L1B radiance OPF version
  ! * Filenames: Slit function, Composite Solar Spectrum, 
  !              Monthly average Solar Spectrum, 
  !              radiance reference
  ! -----------------------------------------------------

  INTEGER (KIND=I4) :: l1br_opf_version

  TYPE pcf_variables
     INTEGER (KIND=I4) :: pge_idx
     CHARACTER (LEN=6) :: pge_name
     INTEGER (KIND=I4) :: verb_thresh_lev
     INTEGER (KIND=I4) :: orbit_number, ecs_version_id
     CHARACTER (LEN=maxchlen), DIMENSION (icf_idx:max_rs_idx) :: static_input_fnames
     CHARACTER (LEN=maxchlen) :: slitfunc_fname, solcomp_filename, solmonthave_filename, &
          l1b_radref_filename, l1b_rad_filename, l1b_irrad_filename, l2_filename, &
          amf_table_filename, climatology_filename, refsec_filename, refsec_cld_filename

  END type pcf_variables
  TYPE(pcf_variables):: pcfvar

  ! ----------------------------------------------------------------------
  ! * Index for the fitting parameters carrying the fitted column value.
  !  N_MOL_FIT:   Number of "molecules" that carry the final column; this
  !               can be one molecule at different temperatures.
  !  FITCOL_IDX:  The main molecule indices, corresponding to the list of
  !               reference spectra.
  !  FINCOL_IDX:  For the final summation of the fitted column: The total
  !               number is the number of different molecules times the
  !               allowed sub-indices. The second dimension is for the
  !               reference spectrum index - this eases the final sum over
  !               he fitted columns (see RADIANCE_FIT subroutine).
  !               includes the subindices, hence the dimension.
  ! N_FINCOL_IDX: Number of final column indices.
  ! * Logicals for processing mode (diagnostic or production), smooth, doas,
  !   normalization, solar composite, common mode, solar_monthly_average, 
  !   radiance reference, radiance reference remove target, use laboratory 
  !   slit function, i0 correction, wavelength dependent AMF, correct O3 AMF
  !   wavelength dependence, run destriping, remove bias, reference sector
  !   correction
  ! * Pixel number limits:
  !    1,2: First and last scan line number
  !    3,4: First and last cross-track pixel
  ! * Latitude limits: Orbit processing, common mode, radiance reference
  ! * Variables connected with ELSUNC numerical precision/convergence criteria
  ! * Solar composite type, Index position of Common Mode add-on, radiance
  !    reference polynomial order, solar fit maximum iteration number, radiance
  !    fit maximum iteration number, radiance wavelength calibration frequency,
  !    wavelength dependent AMF type index
  ! * Array for common mode initial fitting variables
  ! * Solar calibration fitting variables initial value, lower and upped limits
  ! * Solar calibration fitting variables strings
  ! * Logicals for prefitted columns
  ! * Undersampling spectra phase, maximum solar zenith angle, fix land albedo,
  !    fix snow albedo, AMF wavlength, AMF end wavelength, fix cloud albedo,
  !    AMF maximum solar zenith angle, minimum fitting wavelength, maximum
  !    fitting wavelength, maximum good column amount
  ! * Top of Atmosphere [km]
  ! * Undersampling
  ! * Fitting window limits
  ! * Excluded within the fitting window
  ! * Contstraints on the fitting residual: Window and Number of Iterations
  ! ---------------------------------------------------------------------------
  TYPE ctr_variables
     LOGICAL :: yn_diagnostic_run, yn_smooth, yn_doas, yn_spectrum_norm, &
          yn_solar_comp, yn_common_iter, yn_solmonthave, yn_radiance_reference, &
          yn_remove_target, yn_use_labslitfunc, yn_solar_i0, yn_amf_wfmod, &
          yn_o3amf_cor, yn_run_destriping, yn_remove_ctrbias, yn_newshift, &
          yn_refseccor

     LOGICAL, DIMENSION (2) :: yn_bro_prefit, yn_o3_prefit, yn_lqh2o_prefit
     LOGICAL, DIMENSION (us1_idx:us2_idx) :: have_undersampling

     INTEGER (KIND=I4) :: n_mol_fit, n_fincol_idx, solar_comp_typ, common_fitpos, &
          target_npol, max_itnum_sol, max_itnum_rad, radwavcal_freq, amf_wfmod_idx, &
          ctr_bias_pol, ctr_pol_base, ctr_pol_scal, ctr_pol_patt, ctr_nloop, ctr_nblocks, &
          ctr_fitfunc_calls

     INTEGER (KIND=I4), DIMENSION (max_mol_fit) :: fitcol_idx
     INTEGER (KIND=I4), DIMENSION (2,max_mol_fit*mxs_idx) :: fincol_idx
     INTEGER (KIND=I4), DIMENSION (4) :: pixnum_lim
     INTEGER (KIND=i4), DIMENSION (radfit_idx)  :: n_fitres_loop, fitres_range

     REAL (KIND=r4) :: zatmos

     REAL (KIND=r8) :: phase, szamax, amf_alb_lnd, amf_alb_sno, amf_wvl, amf_wvl2, amf_alb_cld, &
          amf_max_sza, winwav_min, winwav_max, max_good_col, tol,  epsrel,  epsabs,  epsx, ctr_maxcol, &
          pm_one

     REAL (KIND=r4), DIMENSION (2) :: radfit_latrange, common_latrange, radref_latrange, &
          fit_winexc_lim, ctrdst_latrange
     REAL (KIND=r8), DIMENSION (3) :: common_fitvar
     REAL (KIND=r8), DIMENSION (max_calfit_idx) :: fitvar_sol_init, lo_sunbnd, up_sunbnd
     REAL (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad_init, lo_radbnd, up_radbnd
     REAL (KIND=r8), DIMENSION (n_fit_winwav) :: fit_winwav_lim

     CHARACTER (LEN=6), DIMENSION (max_calfit_idx) :: fitvar_sol_str
     CHARACTER (LEN=6), DIMENSION (n_max_fitpars)  :: fitvar_rad_str
  END type ctr_variables
  TYPE(ctr_variables) :: ctrvar

  ! -------------------------------------------------
  ! Variables defined in preamble of original program
  ! -------------------------------------------------
  INTEGER (KIND=I4), DIMENSION (n_max_fitpars)  :: mask_fitvar_rad, mask_fitvar_cal, all_radfit_idx
  INTEGER (KIND=I4)                             :: n_fitvar_rad, n_fitvar_cal
  REAL    (KIND=r8), DIMENSION (max_calfit_idx) :: fitvar_cal, fitvar_cal_saved
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad, fitvar_rad_saved, lobnd, upbnd
  REAL    (KIND=r8), DIMENSION (max_rs_idx, nwavel_max) :: database

  ! -----------------------------
  ! Previously IMPLICIT variables
  ! -----------------------------
  REAL (KIND=r8) :: chisq, sol_wav_avg, rad_wav_avg
  REAL (KIND=r8) :: hw1e, e_asym

  ! -----------------------------------------------------------
  ! Variables related to reference spectra
  !  * number of reference spectra:       N_REfSPEC
  !  * indentification strings:           FITPAR_IDXNAME
  !  * file names with reference spectra: REFSPEC_FNAME
  !  * original (uniterpolated) data:     REFSPEC_ORIG_DATA
  !  * number of spectral points:         N_REFSPEC_PTS
  !  * first and last wavelenghts:        REFSPEC_FIRSTLAST_WAV
  ! -----------------------------------------------------------
  INTEGER   (KIND=I4) :: n_refspec

  ! -------------------------------------
  ! TYPE declaration for Reference Specta
  ! -------------------------------------
  TYPE, PUBLIC :: ReferenceSpectrum
     CHARACTER (LEN=maxchlen)                      :: Title, Units
     CHARACTER (LEN=maxchlen)                      :: FileName
     CHARACTER (LEN=maxchlen)                      :: FittingIdxName
     INTEGER   (KIND=I4)                           :: nPoints
     REAL      (KIND=r8)                           :: NormFactor, Temperature
     REAL      (KIND=r8), DIMENSION (2)            :: FirstLastWav
     REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecWavs
     REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecData     
  END TYPE ReferenceSpectrum

  ! -----------------------------------------
  ! TYPE declaration for Common Mode Spectrum
  ! -----------------------------------------
  TYPE, PUBLIC :: CommonModeSpectrum
     CHARACTER (LEN=maxchlen)                                  :: Title, Units
     CHARACTER (LEN=maxchlen)                                  :: FileName
     CHARACTER (LEN=maxchlen)                                  :: FittingIdxName
     INTEGER   (KIND=I4)                                       :: nPoints
     REAL      (KIND=r8)                                       :: NormFactor, Temperature
     REAL      (KIND=r8), DIMENSION (2)                        :: FirstLastWav
     INTEGER   (KIND=I2), DIMENSION (nxtrack_max,2)            :: CCDPixel
     INTEGER   (KIND=I4), DIMENSION (nxtrack_max)              :: RefSpecCount
     REAL      (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts) :: RefSpecWavs
     REAL      (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts) :: RefSpecData     
  END TYPE CommonModeSpectrum

  ! -------------------------------
  ! Array for all Reference Spectra
  ! -------------------------------
  TYPE (ReferenceSpectrum),  DIMENSION (max_rs_idx) :: refspecs_original
  TYPE (CommonModeSpectrum)                         :: common_mode_spec


  REAL (KIND=r8), DIMENSION (max_spec_pts) :: cubic_x, cubic_y, cubic_w

  ! ------------------------------------------------------------------
  ! Indices of fitting window defining wavelengths in current spectrum
  ! ------------------------------------------------------------------
  INTEGER (KIND=I4), DIMENSION (n_fit_winwav) :: fit_winwav_idx

  ! --------------------------------------------------------------------------
  ! The current solar and radiance spectrum, including wavelengths and weights
  ! --------------------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (2)                      :: radiance_wavcal_lnums
  REAL    (KIND=r8), DIMENSION (ccd_idx, nwavel_max)    :: curr_rad_spec
  REAL    (KIND=r8), DIMENSION (ccd_idx, nwavel_max)    :: curr_sol_spec
  REAL    (KIND=r8), DIMENSION (nwavel_max)             :: fitwavs, fitweights, currspec
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: rad_spec_wavcal, rad_wght_wavcal
  
  ! ----------------------------------------------------------------------------
  ! Number of calls to fitting function. This is counted in the fitting function
  ! itself, in an attempt to catch infinite loops of the ELSUNC routine, which
  ! occur on occasion for reasons not yet known.
  ! ----------------------------------------------------------------------------
  INTEGER (KIND=i4) :: num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls

  ! -----------------------------------------------------------------
  ! Generic dimension variables (initialized from either GOME or OMI)
  ! -----------------------------------------------------------------
  INTEGER (KIND=I4) :: n_sol_wvl, n_rad_wvl, n_database_wvl

  ! ---------------------------------------------------------------------------
  ! And this is the convolved solar spectrum. It is (re)initialized in the
  ! solar fit routines only if the above shift&squeeze parameters have changed.
  ! ---------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (max_spec_pts) :: solar_spec_convolved


  ! ---------------------------------------------------------
  ! Filenames specific for the AMF scheme in OMBRO and OMHCHO
  ! ---------------------------------------------------------
  CHARACTER (LEN=maxchlen), DIMENSION (n_voc_amf_luns)   :: voc_amf_filenames

  ! -----------------------------------------------------------------
  ! Logical for Scattering Weights, Gas Profile and Averaging Kernels
  ! Also filename
  ! -----------------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: OMSAO_OMLER_filename

END MODULE OMSAO_variables_module
