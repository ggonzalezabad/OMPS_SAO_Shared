MODULE OMSAO_slitfunction_module

  ! =================================================================
  !
  ! This module defines variables associated with error handling. It
  ! also loads/includes all (SDPTK) files that define error messages
  ! and generally deal with error handling.
  !
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_omidata_module,    ONLY: nxtrack_max
  USE OMSAO_OMPS_SLIT_READER
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! The following two quantities are used to determine whether we need to
  ! reconvolve the solar spectrum. Only if either SHIFT or SQUEEZE have 
  ! changed from one iteration to the other is a reconvolution necessary.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8) :: saved_shift = -1.0E+30_r8, saved_squeeze = -1.0E+30_r8


CONTAINS

  SUBROUTINE omi_slitfunc_convolve ( xtrack_pix, nwvl, wvl, spec_in, specmod, errstat )

    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    xtrack_pix ........... current cross-track pixel number
    !    nwvl ................. number of points in current spectrum
    !    wvl .................. wavelengths of current spectrum
    !    spec_in .............. current spectrum
    !    specmod .............. convolved spectrum
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------

    ! ---------------------------------------
    ! gga, several changes to fix convolution
    ! ---------------------------------------

    CHARACTER (LEN=21), PARAMETER :: modulename = 'omi_slitfunc_convolve'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                  INTENT (IN) :: nwvl, xtrack_pix
    REAL    (KIND=r8), DIMENSION(nwvl), INTENT (IN) :: wvl, spec_in
    
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION(nwvl), INTENT (OUT) :: specmod

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    INTEGER (KIND=i4) :: l, j, j1, j2, k1, k2
    REAL    (KIND=r8) :: sf_area
    LOGICAL           :: yn_full_range
    INTEGER (KIND=i4) :: ntmp
    REAL    (KIND=r8), DIMENSION (nwvl) :: convtmp, spectro

    ! ------------------
    ! OMPS slit function
    ! ------------------
    INTEGER (KIND=i4), PARAMETER       :: maxbands  = 300
    INTEGER (KIND=i4), PARAMETER       :: maxpoints = 501
    REAL    (KIND=r8), PARAMETER       :: startpoint = -2.5, endpoint = 2.5
    REAL    (KIND=r8)                  :: delvar
    LOGICAL                            :: LFAIL
    CHARACTER(LEN=256)                 :: msg
    REAL(KIND=8), DIMENSION(2)         :: wlr
    INTEGER                            :: iXt, ixa, iband
    INTEGER                            :: nwls
    REAL(KIND=8), DIMENSION(maxbands)  :: wls
    INTEGER,      DIMENSION(maxbands)  :: bandIdx
    REAL(KIND=8), DIMENSION(maxpoints) :: xa
    !! result
    REAL(KIND=r8), DIMENSION (nwvl, maxpoints) :: sf_wvals, sf_profiles
    ! The MAX(,) dimension prevents array bound problems for small values of NWVL
    REAL(KIND=r8), DIMENSION (maxpoints)       :: sfwvl_tmp, sfpro_tmp
    

    ! ---------------------------------------------------
    ! Initialize important output and temporary variables
    ! ---------------------------------------------------
    spectro(1:nwvl) = spec_in(1:nwvl)
    locerrstat      = pge_errstat_ok
    specmod(1:nwvl) = 0.0_r8
    sf_area         = 0.0_r8 ; delvar  = 0.0_r8
    sf_wvals        = 0.0_r8 ; sf_profiles = 0.0_r8


    ! ------------------------
    ! Read OMPS Slit functions
    ! ------------------------
    CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
         nbands, nXts, nSiltPts, slitMacroTab,&
         wlMacroTab, wlRelTab, wlRelminTab,   &
         wlRelmaxTab, wlStepTab, LFAIL, msg )

    ! --------------------
    ! Getting band indices
    ! --------------------
    iXt = xtrack_pix
    delvar = wvl(2) - wvl(1)
    wlr(1) = wvl(1) ; wlr(2) = wvl(nwvl)
    CALL OMPS_wlr2bands( wlr, iXt, wls, bandIdx, nwls, LFAIL, msg )

    ! -----------------------------------------------------------
    ! Working out the slit function profiles we are interested in
    ! for each point +- 1.99nm with a resolution of ~0.001. I may
    ! adjust this to improve performace vs accuracy
    ! -----------------------------------------------------------
    delvar = (endpoint - startpoint) / REAL (maxpoints-1, KIND=r8)
    DO ixa = 1, maxpoints
     xa(ixa) = startpoint + (ixa-1.0) * delvar
    END DO
    Do ixa = 1, nwvl
       ! ------------------------------------------------------
       ! Working out the wavelenghts for one band
       ! Find which band is closest to this particular wvl(ixa)
       ! ---------------------------------------
       ! Find the closed band center to wvl(ixa)
       ! ---------------------------------------
       iband = MINLOC(ABS(wls(1:nwls) - wvl(ixa)), 1)
       sf_wvals(ixa,1:maxpoints)  = wvl(ixa) + xa(1:maxpoints)
       sf_profiles(ixa,1:maxpoints) = OMPS_SlitA( xa(1:maxpoints), bandIdx(iband), iXt)
    End Do

    ! -----------------------------------------------------------------------------
    ! In the following loop over the wavelengths in the input spectrum, we perform
    ! the following steps:
    ! (1) Interpolation of the slit function to the spectrum wavelength:
    !     (a) Find indices J1 and J2 in the slit function wavelength array SF_WVALS
    !         bounding the non-zero part of the slit function from left and right;
    !     (b) Find indices K1 and K2 in the spectrum wavelength array WVL such that
    !         SF_WVALS(J1) <= WVL(K1)  < WVL(K2) <= SF_WAVLS(J2)
    !     (c) Interpolation.
    ! (2) Normalize the interpolated slit function to AREA=1.
    ! (3) Convolve the spectrum.
    ! -----------------------------------------------------------------------------
    ! gga

    DO l = 1, nwvl

       convtmp = 0.0_r8
    
       ! ----------------------------------------------
       ! Find indices J1 and J2; add 1 at either end to
       ! assure that we include the ZERO values.
       ! -----------------------------------------------------------------
       ! J1 is the minimum position of SF > 0 minus 1, but not less than 1
       ! -----------------------------------------------------------------
       j1  = 1
       get_sf_start: DO j = 2, maxpoints
          IF ( sf_profiles(l,j) > 0.0_r8 ) THEN
             j1 = j+1 ; EXIT get_sf_start
          END IF
       END DO get_sf_start
       ! ------------------------------------------------------------------------
       ! J2 is the maximum position of SF > 0 plus 1, but not more than N_SF_NPTS
       ! ------------------------------------------------------------------------
       j2  = maxpoints
       get_sf_end: DO j = maxpoints-1, 1, -1
          IF ( sf_profiles(l,j) > 0.0_r8 ) THEN
             j2 = j+1 ; EXIT get_sf_end
          END IF
       END DO get_sf_end
       IF ( j2 <= 0 ) j2 = maxpoints

       ! ----------------------------------------------------------------------------
       ! Now find the bounding indices K1 and K2 in the extended spectrum wavelengths
       ! ----------------------------------------------------------------------------
       CALL array_locate_r8 ( nwvl, wvl(1:nwvl), sf_wvals(l,j1), 'GE', k1 )
       CALL array_locate_r8 ( nwvl, wvl(1:nwvl), sf_wvals(l,j2), 'LE', k2 )
       IF (k1 == k2) THEN
          CYCLE
       END IF

       ! -----------------------------------------------------------------------------
       ! Now interpolate the slit function to the spectrum wavlengths. We are doing it
       ! this way disregard of the number of spectral points in either array (i.e., we
       ! do not interpolate the lower resolution grid to the higher resolution one)
       ! because we don't want to worry about introducing any additional undersampling
       ! effects. This may or may not be an issue but better safe than sorry.
       ! -----------------------------------------------------------------------------
       ! Copy to temporary arrays to avoid creation of same in call to subroutine
       ! ------------------------------------------------------------------------
       ntmp              = j2-j1+1
       sfwvl_tmp(1:ntmp) = sf_wvals   (l,j1:j2)
       sfpro_tmp(1:ntmp) = sf_profiles(l,j1:j2)
       
       ! -------------------------------------------------------
       ! Interpolation of the extended SF to obtain the area for
       ! normalization applied below
       ! -------------------------------------------------------
       CALL interpolation ( &
            modulename, ntmp, sfwvl_tmp(1:ntmp), sfpro_tmp(1:ntmp),                           &
            k2-k1+1, wvl(k1:k2), convtmp(k1:k2), 'endpoints', 0.0_r8, &
            yn_full_range, locerrstat )
       
       ! ------------------------------------------------
       ! If the SF doesn't contribute between nwvl(1) and 
       ! nwvl(nwvl) then we can skip this l.
       ! ------------------------------------------------
       IF (wvl(k2) < wvl(1) .OR. wvl(k1) > wvl(nwvl)) CYCLE
       
       ! --------------------------------------------------------
       ! Renormalize the slit function to AREA = 1, then convolve
       ! --------------------------------------------------------
       sf_area = SUM(convtmp)
       IF ( sf_area == 0.0_r8 ) sf_area = 1.0_r8
       
       ! -----------------------
       ! Save to temporal arrays
       ! -----------------------
       convtmp(k1:k2) = convtmp(k1:k2) / sf_area
              
       ! ----------------------------------------------------------------------------
       ! Done normalization, we perform the convolution of the spectrum with the slit
       ! function.
       ! ----------------------------------------------------------------------------
       specmod(k1:k2) = specmod(k1:k2) + spectro(l) * convtmp(k1:k2)

    END DO
    ! gga 
    
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE omi_slitfunc_convolve

  SUBROUTINE asymmetric_gaussian_sf ( npoints, hw1e, e_asym, wvlarr, specarr, specmod)

    ! =========================================================================
    !
    ! Convolves input spectrum with an asymmetric Gaussian slit function of
    ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                   _
    !                  |               x^2                   |
    !      g(x) =  EXP | - --------------------------------- |
    !                  |_   (hw1e * (1 + SIGN(x)*e_asym))^2 _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================


    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints
    REAL    (KIND=r8),                      INTENT (IN) :: hw1e, e_asym
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr
    
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: specmod

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                        :: i, j, nslit, sslit, eslit
    REAL    (KIND=r8)                        :: slitsum, sliterr, cwvl, lwvl, rwvl
    REAL    (KIND=r8), DIMENSION (3*npoints) :: spc_temp, wvl_temp, sf_val, xtmp, ytmp

    REAL (KIND=r8) :: signdp
    EXTERNAL signdp

    ! --------------------------------------------------------
    ! Initialize output variable (default for "no convolution"
    ! --------------------------------------------------------
    specmod(1:npoints) = specarr(1:npoints)
  
    ! -----------------------------------------------
    ! No Gaussian convolution if Halfwidth @ 1/e is 0
    ! -----------------------------------------------
    IF ( hw1e == 0.0_r8 ) RETURN


    ! ------------------------------------------------------------------------
    ! One temporary variable is SPC_TEMP, which is three times the size of
    ! SPEC. For the convolution routine to work (hopefully) in each and
    ! every case, we reflect the spectrum at its end points to always have a
    ! fully filled slit function. But this causes some real index headaches
    ! when the slit function wraps around at the ends. Performing the mirror
    ! imaging before we get to the convolution helps to keep things a little
    ! more simple.
    !
    ! Note that this approach is the same as for the pre-tabulated OMI lab
    ! slit function. It is adopted here because now we not only convolve the
    ! solar spectrum, but also any higher resolution reference cross sections,
    ! and these may not necessarily be equidistant in wavelength.
    ! ------------------------------------------------------------------------
    spc_temp(npoints+1:2*npoints) = specarr(1:npoints)
    wvl_temp(npoints+1:2*npoints) = wvlarr (1:npoints)
    DO i = 1, npoints
       spc_temp(npoints+1-i) = specarr(i)
       wvl_temp(npoints+1-i) = 2.0_r8*wvlarr(1)-wvlarr(i) -0.001_r8
       spc_temp(2*npoints+i) = specarr(npoints+1-i)
       wvl_temp(2*npoints+i) = 2.0_r8*wvlarr(npoints)-wvlarr(npoints+1-i) +0.001_r8
    END DO

    ! ------------------------------------------------------------------------
    ! We now compute the asymmetric Gaussian for every point in the spectrum.
    ! Starting from the center point, we go outwards and stop accumulating
    ! points when both sides are less than 0.001 of the maximum slit function.
    ! Since we are starting at the center wavelength, this can be set to 1.0.
    ! Remember that the original wavelength array is now located at indices
    ! NPOINTS+1:2*NPOINTS
    ! ------------------------------------------------------------------------
    DO i = 1, npoints
       sf_val = 0.0_r8
       cwvl = wvl_temp(npoints+i)

       sf_val(npoints+i) = 1.0_r8
       getslit: DO j = 1, npoints
          sslit = npoints+i-j ; lwvl = cwvl - wvl_temp(sslit)
          eslit = npoints+i+j ; rwvl = cwvl - wvl_temp(eslit)

          sf_val(sslit) = EXP(-lwvl**2 / ( hw1e * (1.0_r8 + signdp(lwvl)*e_asym) )**2)
          sf_val(eslit) = EXP(-rwvl**2 / ( hw1e * (1.0_r8 + signdp(rwvl)*e_asym) )**2)
          IF ( sf_val(sslit) < 0.0005_r8 .AND. sf_val(sslit) < 0.0005_r8 ) EXIT getslit
       END DO getslit

       ! ----------------------------------
       ! The number of slit function points
       ! ----------------------------------
       nslit = eslit - sslit + 1
       ! ----------------------------------------------------------------
       ! Compute the norm of the slitfunction. It should be close to 1
       ! already, but making sure doesn't hurt.
       ! ----------------------------------------------------------------
       xtmp(1:nslit) = wvl_temp(sslit:eslit)-cwvl
       ytmp(1:nslit) = sf_val  (sslit:eslit)
       CALL cubint ( &
            nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, slitsum, sliterr)

       IF ( slitsum > 0.0_r8 ) sf_val(sslit:eslit) = sf_val(sslit:eslit) / slitsum

       ! ---------------------------------------------------------------------
       ! Prepare array for integration: Multiply slit function values with the
       ! spectrum array to be convolved.
       ! ---------------------------------------------------------------------
       ytmp(1:nslit) = sf_val(sslit:eslit) * spc_temp(sslit:eslit)

       ! ----------------------------------------------------------
       ! Folding (a.k.a. integration) of spectrum and slit function
       ! ----------------------------------------------------------
       CALL cubint ( &
            nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, specmod(i), sliterr)
    END DO

    RETURN
  END SUBROUTINE asymmetric_gaussian_sf

END MODULE OMSAO_slitfunction_module
