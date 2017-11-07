MODULE OMSAO_slitfunction_module

  ! =================================================================
  !
  ! This module defines variables associated with error handling. It
  ! also loads/includes all (SDPTK) files that define error messages
  ! and generally deal with error handling.
  !
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_data_module,    ONLY: nxtrack_max
  USE OMPS_macroSlitNOAA_m
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! The following two quantities are used to determine whether we need to
  ! reconvolve the solar spectrum. Only if either SHIFT or SQUEEZE have 
  ! changed from one iteration to the other is a reconvolution necessary.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8) :: saved_shift = -1.0E+30_r8, saved_squeeze = -1.0E+30_r8

CONTAINS

  SUBROUTINE omi_slitfunc_convolve ( xtrack_pix, nwvl, wvl, spec_in, specmod, stretch, errstat )

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

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                  INTENT (IN) :: nwvl, xtrack_pix
    REAL    (KIND=r8), DIMENSION(nwvl), INTENT (IN) :: wvl, spec_in
    REAL    (KIND=r8),                  INTENT (IN) :: stretch
    
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
    INTEGER (KIND=i4) :: j

    ! ------------------
    ! OMPS slit function
    ! ------------------
    REAL (KIND=r8), PARAMETER :: dhalf = 2.0
    REAL (KIND=r8) :: delwvl, slitsum
    INTEGER :: fidx, lidx, nhalf, i
    INTEGER :: sidx, eidx, nslit, fpnt, lpnt
    REAL(KIND=8), DIMENSION(nwvl) :: locsli
    INTEGER, DIMENSION(nwvl) :: idxs    

    ! ---------------------------------------------------
    ! Initialize important output and temporary variables
    ! ---------------------------------------------------
    locerrstat = pge_errstat_ok
    specmod(1:nwvl) = spec_in(1:nwvl)
    fidx = 1; lidx = nwvl
    delwvl = wvl(fidx+1)-wvl(fidx)
    nhalf = CEILING(dhalf / delwvl)
    nslit = nhalf * 2 + 1
    sidx = MAX(fidx, nhalf + 1)
    eidx = MIN(lidx, nwvl - nhalf)

    ! Safe convolution
    DO i = sidx, eidx
       CALL compute_slitprofile (xtrack_pix, stretch, wvl(i),nslit, wvl(i-nhalf:i+nhalf), &
            locsli(1:nslit))
       slitsum = SUM(locsli(1:nslit))
       specmod(i) = DOT_PRODUCT(locsli(1:nslit), spec_in(i-nhalf:i+nhalf)) / slitsum
    END DO

    CALL compute_slitprofile (xtrack_pix, stretch, wvl(sidx),nslit, wvl(sidx-nhalf:sidx+nhalf), &
         locsli(1:nslit))
    slitsum = SUM(locsli(1:nslit))

    ! Left wrap mode
    DO i = fidx, sidx - 1
       idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
       lpnt = nhalf - i + 1
       idxs(1 : lpnt) = ABS(idxs(1 : lpnt)) + 2
       specmod(i) = DOT_PRODUCT(locsli(1:nslit), spec_in(idxs(1:nslit))) / slitsum
    ENDDO

    CALL compute_slitprofile (xtrack_pix, stretch, wvl(eidx),nslit, wvl(eidx-nhalf:eidx+nhalf), &
         locsli(1:nslit))
    slitsum = SUM(locsli(1:nslit))

    ! Right wrap mode
    DO i = eidx + 1, lidx
       idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
       fpnt = nhalf + (nwvl - i) + 1
       idxs(fpnt:nslit) =  nwvl * 2 - idxs(fpnt:nslit)
       specmod(i) = DOT_PRODUCT(locsli(1:nslit), spec_in(idxs(1:nslit))) / slitsum
    ENDDO
    
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE omi_slitfunc_convolve

  SUBROUTINE OMPS_Slit_DRIVER (ixt, wlc, scaling,  nslit, xa, slit )
    USE OMPS_macroSlitNOAA_m, ONLY : OMPS_SlitF, OMPS_SlitA, &
         OMPS_wlr2bands, OMPS_wl2band
    
    IMPLICIT NONE
    ! INPUT VARIABLES
    INTEGER, INTENT (IN)      :: iXt     ! idx of cross-track position
    REAL (KIND=8), INTENT(IN) :: wlc     ! band center
    REAL (KIND=8), INTENT(IN) :: scaling ! shape parameter for broadening/squeezing the slit function
    INTEGER, INTENT (IN)      :: nslit   ! the number of slit    
    REAL (KIND=8), INTENT(IN), DIMENSION (nslit) :: xa ! slit position relative to band center
    ! OUTPUT VARIABLES
    REAL(KIND=8), INTENT(OUT), DIMENSION(nslit) :: slit ! slitfunction
    ! LOCAL VARIABLES
    INTEGER :: ii, ib
    LOGICAL :: LFAIL
    CHARACTER(LEN=255) :: msg
    
    !! OMPS slit is tabulated according to cross track 1 to 36,
    IF ( ixt < 1 .or. ixt > 36 ) THEN
       write (*,'(A)') ' check the pixel range in OMPS Slit Deriver'
       stop
    ENDIF
    !! band index (1 to 196)
    !! Given a band center wavelength index, this function
    !! returns the band index to be used
    ib =  OMPS_wl2band( wlc, iXt, LFAIL, msg )
    IF ( ib < 1 .or. ib > 196 ) THEN
       write (*,'(A)') ' check the wave range in OMPS Slit Deriver'
       stop
    ENDIF
    !! scaling = 1.0, OMPS_SlitF is equal to database.
    !! scaling > 1.0, slitF is broadened.
    !! scaling < 1.0, slitF is squeezed.
    !! Normally sclaing should not go below 0.5 or above 1.5
    DO ii = 1, nslit
       slit(ii) = OMPS_SlitF( xa(ii), ib, iXt, scaling )
    ENDDO
  END SUBROUTINE OMPS_Slit_DRIVER
  
  SUBROUTINE compute_slitprofile (currpix, scailing, bandcenter, nslit, slitwave, slit)
    
    IMPLICIT NONE
    ! Input/Output variables   
    INTEGER, INTENT(IN) :: nslit, currpix
    REAL (KIND=8), INTENT(IN) :: bandcenter, scailing
    REAL (KIND=8), DIMENSION(nslit), INTENT(IN)  :: slitwave
    REAL (KIND=8), DIMENSION(nslit), INTENT(OUT) :: slit
    ! Local variables
    REAL(KIND=8), DIMENSION(nslit) :: xa
    
    IF ( scailing < 0.1 .and. scailing > 2.0 ) THEN
       print * , 'check omps_slit_scaling in omps_slit_module'
       stop
    ENDIF
    xa(1:nslit) = slitwave(1:nslit) - bandcenter
    call omps_slit_driver (currpix, bandcenter, scailing, nslit, xa, slit)
    RETURN
  END SUBROUTINE compute_slitprofile

  SUBROUTINE super_gaussian_sf ( npoints, hw1e, e_asym, g_shap, wvlarr, specarr, specmod)

    ! =========================================================================
    !
    ! Convolves input spectrum with an asymmetric Gaussian slit function of
    ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                            _
    !                  |   |            x^2                  |^g_shap |
    !      g(x) =  EXP | - |---------------------------------|        |
    !                  |_  | (hw1e * (1 + SIGN(x)*e_asym))   |       _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================


    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints
    REAL    (KIND=r8),                      INTENT (IN) :: hw1e, e_asym, g_shap
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
          sslit = npoints+i-j ; lwvl = - cwvl + wvl_temp(sslit)
          eslit = npoints+i+j ; rwvl = - cwvl + wvl_temp(eslit)
          sf_val(sslit) = EXP(-(ABS(lwvl / ( hw1e + signdp(lwvl)*e_asym ) ) )**g_shap )
          sf_val(eslit) = EXP(-(ABS(rwvl / ( hw1e + signdp(rwvl)*e_asym ) ) )**g_shap ) 
          IF ( sf_val(eslit) < 0.0005_r8 .AND. sf_val(sslit) < 0.0005_r8 ) EXIT getslit
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
  END SUBROUTINE super_gaussian_sf

END MODULE OMSAO_slitfunction_module
