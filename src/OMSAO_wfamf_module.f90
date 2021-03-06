MODULE OMSAO_wfamf_module

  ! ====================================================================
  ! This module defines variables associated with the wavelength depende
  ! nt AMF calculations and contains necessary subroutines to read files
  ! and calculate them
  ! ====================================================================
  USE OMSAO_data_module, ONLY: n_roff_dig
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, &
       pge_errstat_error, pge_errstat_fatal, &
       vb_lev_default, vb_lev_screen, he5_stat_fail, omi_s_success, &
       omsao_e_he5swopen, omsao_e_prefitcol, omsao_e_he5swlocate, &
       omsao_e_he5swattach, omsao_e_he5gdattach, omsao_e_he5gdopen, &
       omsao_w_he5gdclose, f_sep, error_check
  USE OMSAO_he5_module
  USE OMSAO_variables_module, ONLY: ctrvar, pcfvar

  IMPLICIT NONE

  ! ---------------------------------------
  ! Data obtained from the climatology file
  ! ---------------------------------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: latvals, lonvals, Ap, Bp
  REAL(KIND=r4), DIMENSION(:,:), ALLOCATABLE :: Psurface
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: Temperature, Gas_profiles

  ! ---------------------------------------
  ! Data obtained from Vlidort lookup table
  ! ---------------------------------------
  ! --------------------------------------------------
  ! Parameter for the definition of the vlidort arrays
  ! --------------------------------------------------
  ! ------------------------
  ! Cross sections variables
  ! ------------------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_OzC0, vl_OzC1, vl_OzC2
  
  ! --------------
  ! Grid variables
  ! --------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_pre
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_sza
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_vza
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_wav
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: vl_toms

  ! ------------------
  ! Profiles variables
  ! ------------------
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: vl_air, vl_alt, vl_ozo, vl_tem  
  
  ! --------------------------
  ! Parameterization variables
  ! --------------------------
  REAL(KIND=r4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: vl_I0, vl_I1, vl_I2, vl_Ir
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: vl_Sb  
  REAL(KIND=r4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: vl_dI0, vl_dI1, vl_dI2, vl_dIr
  REAL(KIND=r4)  :: vl_Factor

  ! -------------------
  ! Dimension variables
  ! -------------------
  INTEGER(KIND=i4) :: vl_nozo, vl_ncld, vl_nsza, vl_nvza, vl_nwav, vl_nalt
  
  ! ---------
  ! PCF stuff
  ! ---------
  INTEGER(KIND=i4), PARAMETER :: climatology_lun = 700270

  ! -----------------------------
  ! Dimensions of the climatology
  ! -----------------------------
  INTEGER (KIND=i4) :: Cmlat, Cmlon, CmETA, CmEp1

  ! -----------------------
  ! To find the right swath
  ! -----------------------
  INTEGER   (KIND=i4), PARAMETER, PRIVATE :: nmonths = 12
  CHARACTER (LEN=9), DIMENSION (nmonths), PARAMETER, PRIVATE :: months = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December '    /)

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER, PRIVATE :: zerocl = 0, onecl = 1

  ! -----------
  ! ISCCP stuff
  ! -----------
  ! --------------------------------------------
  ! TYPE declaration for ISCCP cloud climatology
  ! --------------------------------------------
  INTEGER (KIND=i4), PARAMETER, PRIVATE :: nlat_isccp=72, nlon_isccp=6596
  TYPE, PUBLIC :: CloudClimatology
     REAL (KIND=r8) :: scale_ctp, scale_cfr, delta_lat, &
          missval_cfr, missval_ctp
     REAL (KIND=r4), DIMENSION (nlat_isccp) :: latvals, delta_lon
     INTEGER (KIND=i4), DIMENSION (nlat_isccp) :: n_lonvals
     REAL (KIND=r4), DIMENSION (nlon_isccp) :: lonvals
     REAL (KIND=r8), DIMENSION (nlon_isccp) :: cfr, ctp
  END TYPE CloudClimatology

  ! --------------------------
  !(3) ISCCP Cloud Climatology
  ! --------------------------
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lat_field  = 'ISCCP_Lats'
  CHARACTER (LEN=15), PARAMETER, PRIVATE :: isccp_dlon_field = 'ISCCP_DeltaLons'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: isccp_nlon_field = 'ISCCP_NumLons'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lon_field  = 'ISCCP_Lons'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mcfr_field = 'ISCCP_MonthlyAVG_CloudFraction'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mctp_field = 'ISCCP_MonthlyAVG_CloudPressure'


CONTAINS

  SUBROUTINE amf_calculation_bis (               &
       pge_idx, cloud_filename,                  &
       nt, nx, lat, lon, sza, vza,               &
       xtrange, saocol, saodco, saoamf, amfdiag, &
       terrain_height, yn_write, errstat )

    ! =================================================================
    ! This subroutine computes the AMF factor using the following eleme
    ! nts:
    !     - Kleipool OMLER database
    !     - GEOS Chem climatology
    !     - VLIDORT calculated scattering weights
    ! =================================================================
    USE OMSAO_data_module, ONLY: albedo
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx, pge_idx
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza, &
         terrain_height
    LOGICAL :: yn_write
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
    CHARACTER (LEN=*), INTENT (IN) :: cloud_filename

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saocol, saodco, saoamf
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, itt, spixx, epixx
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1) :: amfgeo
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1) :: l2cfr, l2ctp
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1) :: local_albedo, cli_psurface
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: climatology, cli_temperature, &
         cli_heights
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: scattw

    locerrstat  = pge_errstat_ok

    ! ------------------------------------
    ! Initialize variables that are output
    ! ------------------------------------
    local_albedo       = r8_missval
    climatology  = r8_missval
    cli_heights  = r8_missval
    cli_psurface = r8_missval
    scattw       = r8_missval
    saoamf       = r8_missval
    amfgeo       = r8_missval
    amfdiag      = i2_missval

    ! -----------------------------------------
    ! If amf_wvl < 0.0 then the slant column is
    ! reported and AMFs equal to 1 with scattw, 
    ! akernels and climatology set to missval
    ! -----------------------------------------
    IF (ctrvar%amf_wvl .LT. 0.0) THEN

       DO itt = 0, nt-1
          spixx = xtrange(itt,1) ; epixx = xtrange(itt,2)
          saoamf(spixx:epixx,itt) = 1.0_r8
       END DO

    ELSE

       ! -------------------------
       ! Compute the geometric AMF
       ! -------------------------
       CALL compute_geometric_amf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )
       
       ! -------------------------------------------------------
       ! Initialize molecular AMF with geometric AMF. Subsequent
       ! subroutines will replace any entries where the true
       ! molecular AMF can be computed.
       ! -------------------------------------------------------
       saoamf = amfgeo
          
       ! ---------------------------------------
       ! Write the albedo to the output file he5
       ! ---------------------------------------
       local_albedo(1:nx,0:nt-1) = albedo(1:nx,0:nt-1)
       IF (yn_write) CALL write_albedo_he5 ( local_albedo, nt, nx, locerrstat)
       
       ! ------------------------------
       ! Read the OMPS L2 cloud product
       ! ------------------------------
       locerrstat = pge_errstat_ok
       CALL amf_read_ompsclouds ( cloud_filename, nt, nx, l2cfr, l2ctp, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       IF ( locerrstat >= pge_errstat_error ) THEN
          l2cfr = r8_missval
          l2ctp = r8_missval
       END IF

       ! ---------------------------------------------------------------------
       ! The climatology has already been read, inside omi_pge_fitting_process
       ! Now it is only needed to interpolate to the pixels of the granule.
       ! It was read there to obtain the dimensions of the number of levels.
       ! ---------------------------------------------------------------------
       CALL omi_climatology (climatology, cli_heights, cli_psurface, cli_temperature, lat, lon, &
            nt, nx, xtrange, locerrstat)
       
       ! -------------------------------------
       ! Write the climatology to the he5 file
       ! -------------------------------------
       IF (yn_write) CALL write_climatology_he5 (climatology, cli_heights, nt, nx, CmETA, locerrstat)

       ! ------------------------------------------------------------------
       ! Read VLIDORT look up table. Variables are declared at module level
       ! ------------------------------------------------------------------
       IF (yn_write) CALL read_vlidort (locerrstat)

       ! ----------------------------------------------------------------------
       ! amfdiag is used to keep track of the pixels were enough information is
       ! available to carry on the AMFs calculation.
       ! ----------------------------------------------------------------------
       CALL amf_diagnostic ( nt, nx, sza, vza, xtrange, &
            l2cfr, l2ctp, amfdiag  )

       ! --------------------------------------------------------
       ! Compute Scattering weights in the look up table grid but
       ! with the correct albedo. amfdiag is used to skip pixel
       ! ---------------------------------------------------------
       CALL compute_scatt ( nt, nx, local_albedo, sza, vza, l2ctp, l2cfr, terrain_height, cli_heights, amfdiag, &
            scattw)

       ! ----------------------------
       ! Deallocate Vlidort variables
       ! ----------------------------
       IF (.NOT. yn_write) CALL vlidort_allocate ("d", vl_nozo, vl_ncld, vl_nsza, vl_nvza, vl_nwav, vl_nalt, &
            errstat)

       ! -----------------------------------------------------------------
       ! Work out the AMF using the scattering weights and the climatology
       ! Work out Averaging Kernels
       ! -----------------------------------------------------------------
       CALL compute_amf ( nt, nx, CmETA, climatology, scattw, saoamf, amfdiag)
       
       ! -----------------------------------------------------------------
       ! Write out scattering weights, altitude grid and averaging kernels
       ! -----------------------------------------------------------------
       IF (yn_write) CALL write_scatt_he5 (scattw, nt, nx, CmETA, locerrstat)
       
    END IF

    ! --------------------------
    ! Apply the air mass factors
    ! --------------------------
    WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval ) 
       saocol = saocol / saoamf
       saodco = saodco / saoamf
    END WHERE

    ! -----------------------------------------------
    ! Write AMFs, AMF diagnosting, and AMF-adjusted
    ! columns and column uncertainties to output file
    ! -----------------------------------------------
    IF (yn_write) CALL he5_amf_write ( pge_idx, nx, nt, saocol, saodco, saoamf, &
         amfgeo, amfdiag, l2cfr, l2ctp, locerrstat )
    
    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE amf_calculation_bis
  
  SUBROUTINE omi_climatology (climatology, local_heights, local_psurf, local_temperature, &
                              lat, lon, nt, nx, xtrange, locerrstat)

    ! =========================================
    ! Extract Gas climatology to granule pixels
    ! No interpolation or something like that,
    ! Just pick the closest model grid
    ! =========================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: locerrstat
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1, CmETA), INTENT (INOUT) :: &
         climatology, local_temperature, local_heights
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (INOUT) :: local_psurf  

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itimes, ixtrack, spix, epix, idx_lat, idx_lon, ilevel, n, n1
    REAL (KIND=r8) :: rho1, rho2, grad, aircolumn
    REAL (KIND=r8), DIMENSION (0:CmETA) :: lhgt, lpre, ltmp

    ! -----------------------
    ! Some physical constants
    ! -----------------------
    REAL (KIND=r8), PARAMETER ::                &
         rho_stand = 2.6867773e+19_r8,          & ! Loschmidt, in [cm^-3]
         pzero     = 1013.25_r8,                & ! P0 in [mb]
         tzero     = 273.15_r8,                 &
         rho_zero  = rho_stand * tzero / pzero

    ! -------------------------------
    ! Air density conversion constant
    ! -------------------------------
    REAL (KIND=r8), PARAMETER :: km2cm  = 1.0E+05_r8

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok
    
    ! ------------------------------------------------------------
    ! Find the Climatology corresponding to each lat and lon pixel
    ! ------------------------------------------------------------
    DO itimes = 0, nt-1

       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       DO ixtrack = spix, epix

          IF (lon(ixtrack,itimes) .LT. -180.0_r4 .OR. &
              lat(ixtrack,itimes) .LT.  -90.0_r4 .OR. &
              lon(ixtrack,itimes) .GT.  180.0_r4 .OR. &
              lat(ixtrack,itimes) .GT.  90.0_r4) CYCLE

          ! ----------------------------------------------------
          ! Just selecting the closest location to lat lon pixel
          ! ----------------------------------------------------
          idx_lat = MINVAL(MINLOC(ABS(latvals(1:Cmlat) - lat(ixtrack,itimes) )))
          idx_lon = MINVAL(MINLOC(ABS(lonvals(1:Cmlon) - lon(ixtrack,itimes) )))

          climatology(ixtrack,itimes,1:CmETA)       = Gas_profiles(idx_lon,idx_lat,1:CmETA)
          local_temperature(ixtrack,itimes,1:CmETA) = Temperature(idx_lon,idx_lat,1:CmETA)
          local_psurf(ixtrack,itimes)               = Psurface(idx_lon,idx_lat)
          DO ilevel = 1, CmETA
             local_heights(ixtrack,itimes,ilevel)   = ( ( Ap(ilevel) + ( local_psurf(ixtrack,itimes) * Bp(ilevel)   ) ) + &
                                                      ( Ap(ilevel+1) + ( local_psurf(ixtrack,itimes) * Bp(ilevel+1) ) ) ) / 2 !(in hPa)
          END DO

          ! -------------------
          ! Develop air density
          ! -------------------
          ! --------------------------------------------------------
          ! Some notes on Units
          ! -------------------
          ! Values from the atmospheric profile now have these units:
          !   pressures:    hPa = mb
          !   temperatures: K
          !   heights       km
          !
          ! CONST: 10^-1 K m^-3 hPa^-1
          !
          ! ==> aircolumns: 10^4 m^-2  or 10^10 km^-2  or  cm^-2
          ! --------------------------------------------------------
          lhgt              = 0.0_r8
          lpre              = 0.0_r8
          ltmp              = 0.0_r8
          grad              = 0.0_r8

          DO n = 0, CmETA-1
             lhgt(n) = -16.0_r8 * log10(local_heights(ixtrack,itimes,n+1) / 1013.0_r8)
          END DO

          ! Create a fake top of atmos 1 km above top altitude of the
          ! climatology
          lhgt(CmETA)     = local_heights(ixtrack,itimes,CmETA) / 1000.0_r8 + 1.0_r8
          lpre(0:CmETA-1) = local_heights(ixtrack,itimes,1:CmETA)
          ltmp(0:CmETA-1) = local_temperature(ixtrack,itimes,1:CmETA)

          ! --------------------------------------------
          ! extrapolated temperature and pressure at TOA
          ! --------------------------------------------
          grad        =  LOG ( lpre(CmETA-2)/lpre(CmETA-1) ) / (lhgt(CmETA-2) - lhgt(CmETA-1))
          lpre(CmETA) = EXP( LOG( lpre(CmETA-1)) + grad * (lhgt(CmETA) - lhgt(CmETA-1) ) )
          grad        = (ltmp(CmETA-2) - ltmp(CmETA-1)) / (lhgt(CmETA-2) - lhgt(CmETA-1))
          ltmp(CmETA) = ltmp(CmETA-1) + grad*(lhgt(CmETA) - lhgt(CmETA-1))
          
          DO n = 1, CmETA

             rho1              = 0.0_r8
             rho2              = 0.0_r8
             aircolumn         = 0.0_r8

             n1                   = n - 1
             rho1                 = lpre(n1) / ltmp(n1)
             rho2                 = lpre(n)  / ltmp(n)
             aircolumn            = 0.5_r8*rho_zero*(rho1+rho2)*(lhgt(n)-lhgt(n1))*km2cm

             ! -------------------------------------------------------------
             ! Divided by 1E09 to convert from ppb in the climatology to VMR
             ! -------------------------------------------------------------
             climatology(ixtrack,itimes,n) = aircolumn * climatology(ixtrack,itimes,n) / 1.0E09

          END DO

          !  Set non-physical entries to zero.
          WHERE ( climatology(ixtrack,itimes,1:CmETA) < 0.0_r8 )
             climatology(ixtrack,itimes,1:CmETA) = 0.0_r8
          END WHERE

       END DO
    END DO

  END SUBROUTINE omi_climatology

  SUBROUTINE read_climatology(errstat)
    
    ! ==========================================================
    ! This subroutine reads in the climatology from GEOS-Chem or
    ! other source. The climatology needs to be produced keeping
    ! the format so the reader will be flexible enough to move
    ! from one climatology to another one
    ! ==========================================================
    
    USE OMSAO_indices_module,   ONLY: sao_molecule_names
    IMPLICIT NONE    

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: nswath, locerrstat, swath_id, swath_file_id, swlen, he5stat, &
         ismonth, ndatafields
    INTEGER (KIND=i4), DIMENSION(10) :: datafield_rank, datafield_type
    INTEGER (KIND=C_LONG) :: nswathcl, Cmlatcl, Cmloncl, CmETAcl, CmEp1cl
    CHARACTER (LEN=maxchlen) :: swath_file, locswathname, gasdatafieldname, &
         datafield_name
    CHARACTER (LEN=10*maxchlen) :: swath_name

    CHARACTER (LEN= 9), PARAMETER :: cli_lat_field = 'Latitudes'
    CHARACTER (LEN=10), PARAMETER :: cli_lon_field = 'Longitudes'
    CHARACTER (LEN= 3), PARAMETER :: cli_Ap_field = 'A_p'
    CHARACTER (LEN= 3), PARAMETER :: cli_Bp_field = 'B_p'
    CHARACTER (LEN=15), PARAMETER :: cli_Psurf_field  = 'SurfacePressure'
    CHARACTER (LEN=18), PARAMETER :: cli_Temperature_field = 'TemperatureProfile'

    REAL (KIND=r4) :: scale_lat, scale_lon, scale_Ap, scale_Bp, scale_gas, scale_Psurf, &
         scale_temperature

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=16), PARAMETER :: modulename = 'read_climatology' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok
    
    swath_file = TRIM(ADJUSTL(pcfvar%climatology_fname))
    ismonth    = granule_month

    ! --------------------------------------------------------------
    ! Open he5 OMI climatology and check SWATH_FILE_ID (-1 if error)
    ! --------------------------------------------------------------
    swath_file_id = HE5_SWOPEN (swath_file, he5f_acc_rdonly)
    IF (swath_file_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF
    
    ! -----------------------------------------------------------
    ! Check for existing HE5 swathw and attach to the one we need
    ! -----------------------------------------------------------
    nswathcl = HE5_SWinqswath(TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    nswath   = INT(nswathcl, KIND=i4 )

    ! ----------------------------------------------------------------
    ! If there is only one swath in the file, we can attach to it but
    ! if there are more (NSWATH > 1), then we must find the swath that
    ! corresponds to the current month.
    ! ----------------------------------------------------------------
    IF (nswath > 1) THEN
       CALL extract_swathname(nswath, TRIM(ADJUSTL(swath_name)), &
            TRIM(ADJUSTL(months(granule_month))), locswathname)
       ! ---------------------------------------------------------------------------
       ! Check if we found the correct swath name. If not, report an error and exit.
       ! ---------------------------------------------------------------------------
       IF ( INDEX (TRIM(ADJUSTL(locswathname)),TRIM(ADJUSTL(months(ismonth)))) == 0 ) THEN
          CALL error_check ( &
               0, 1, pge_errstat_error, OMSAO_E_HE5SWLOCATE, modulename, &
               vb_lev_default, locerrstat )
          errstat = MAX ( errstat, locerrstat )
          RETURN
       END IF
    ELSE
       locswathname = TRIM(ADJUSTL(swath_name))
    END IF

    ! -----------------------------
    ! Attach to current month swath
    ! -----------------------------
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(locswathname)) )
    IF ( swath_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! ------------------------------------
    ! Read dimensions of Climatology swath
    ! ------------------------------------
    locerrstat = pge_errstat_ok
    CALL climatology_getdim ( swath_id, Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
    IF ( ismonth < 1 .OR. ismonth > nmonths .OR. locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       RETURN
    END IF

    ! ----------------------------------------------------------------------------
    ! Create KIND=4/KIND=8 variables. We have to use those dimensions a few times
    ! in this subroutine, so it saves some typing if we do the conversion once and
    ! save them in new variables.
    ! ----------------------------------------------------------------------------
    Cmlatcl = INT ( Cmlat, KIND=C_LONG )
    Cmloncl = INT ( Cmlon, KIND=C_LONG )
    CmETAcl = INT ( CmETA, KIND=C_LONG )
    CmEp1cl = INT ( CmEp1, KIND=C_LONG )

    ! ---------------------------
    ! Allocate Climatology arrays
    ! ---------------------------
    locerrstat = pge_errstat_ok
    CALL climatology_allocate ( "a", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       CALL climatology_allocate ( "d", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
       RETURN
    END IF

    ! -------------------------------
    ! Read dimension-defining arrays
    ! -------------------------------
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = Cmlatcl
    he5stat = HE5_SWrdfld ( swath_id, cli_lat_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, latvals(1:Cmlat) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = Cmloncl
    he5stat = HE5_SWrdfld ( swath_id, cli_lon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, lonvals(1:Cmlon) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = CmEp1cl
    he5stat = HE5_SWrdfld ( swath_id, cli_Ap_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, Ap(1:CmEp1) )
    he5stat = HE5_SWrdfld ( swath_id, cli_Bp_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, Bp(1:CmEp1) )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology arrays access failed.', vb_lev_default, errstat )

    ! -----------------------------------------------
    ! Read dimension-defining scale factor attributes
    ! -----------------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_lat_field, "ScaleFactor", scale_lat )
    he5stat = HE5_SWrdlattr ( swath_id, cli_lon_field, "ScaleFactor", scale_lon )
    he5stat = HE5_SWrdlattr ( swath_id, cli_Ap_field,  "ScaleFactor", scale_Ap )
    he5stat = HE5_SWrdlattr ( swath_id, cli_Bp_field,  "ScaleFactor", scale_Bp )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology attributes access failed.', vb_lev_default, errstat )

    ! -----------------------------------
    ! Apply scaling factors to geo fields
    ! -----------------------------------
    lonvals = lonvals * scale_lon
    latvals = latvals * scale_lat
    Ap      = Ap      * scale_Ap
    Bp      = Bp      * scale_Bp

    ! -----------------------------------------------
    ! Read the tables: Psurface, Heights, Temperature
    ! -----------------------------------------------
    he5_start_2d  = (/ zerocl, zerocl /)
    he5_stride_2d = (/  onecl,  onecl /)
    he5_edge_2d   = (/ Cmloncl, Cmlatcl /)

    he5_start_3d  = (/ zerocl, zerocl, zerocl/)
    he5_stride_3d = (/  onecl,  onecl,  onecl /)
    he5_edge_3d   = (/ Cmloncl, Cmlatcl, CmETAcl /)
    he5stat = HE5_SWrdfld (                                &
         swath_id, cli_Psurf_field,                        &
         he5_start_2d, he5_stride_2d, he5_edge_2d,         &
         Psurface(1:Cmlon,1:Cmlat) )
    he5stat = HE5_SWrdfld (                                &
         swath_id, cli_Temperature_field,                  &
         he5_start_3d, he5_stride_3d, he5_edge_3d,         &
         Temperature(1:Cmlon,1:Cmlat,1:CmETA) )

    ! --------------------------------------
    ! Read datafields scale factor attribute
    ! --------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_Psurf_field,       "ScaleFactor", scale_Psurf       )
    he5stat = HE5_SWrdlattr ( swath_id, cli_Temperature_field, "ScaleFactor", scale_Temperature )

    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology data fields access failed.', vb_lev_default, errstat )

    ! -----------------------------------------------------------------------
    ! Finding out the data field for the gas of interest (.eq. to target gas)
    ! -----------------------------------------------------------------------
    ndatafields = HE5_swinqdflds(swath_id, datafield_name, datafield_rank, datafield_type)
    CALL extract_swathname(nswath, TRIM(ADJUSTL(datafield_name)), &
         TRIM(ADJUSTL(sao_molecule_names(pcfvar%pge_idx))), gasdatafieldname)

    ! ---------------------------------------------------------------------------
    ! Check if we found the correct swath name. If not, report an error and exit.
    ! ---------------------------------------------------------------------------
    IF ( INDEX (TRIM(ADJUSTL(gasdatafieldname)),TRIM(ADJUSTL(sao_molecule_names(pcfvar%pge_idx)))) == 0 ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWLOCATE, modulename, &
            vb_lev_default, locerrstat )
       WRITE(*,*) "Climatology file does not contain data for ", sao_molecule_names(pcfvar%pge_idx)
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! -----------------------------
    ! Read data from this datafield
    ! -----------------------------
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(gasdatafieldname)),        &
         he5_start_3d, he5_stride_3d, he5_edge_3d,         &
         Gas_profiles(1:Cmlon,1:Cmlat,1:CmETA) )

    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(gasdatafieldname)),&
              "ScaleFactor", scale_gas       )

    ! ------------------------------------
    ! Apply scaling factors to data fields
    ! ------------------------------------
    Temperature  = Temperature  * scale_Temperature
    Psurface     = Psurface     * scale_Psurf
    Gas_profiles = Gas_profiles * scale_gas

  END SUBROUTINE read_climatology

  SUBROUTINE omi_omler_albedo( lat, lon, albedo, nt, nx,  &
                               errstat)

    ! ==================================================================
    ! This subroutine reads the OMLER albedo data base for the month of
    ! the orbit to processed. Then it interpolates the values for each
    ! one of the pixels of the orbit to be analyzed
    ! ==================================================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (INOUT) :: albedo

    ! ------------------------------------------------------------------
    ! Local variables, the variables to hold the OMLER data are going to
    ! be allocated and deallocated within this subroutine.
    ! ------------------------------------------------------------------
    REAL (KIND=r4), ALLOCATABLE, DIMENSION(:) :: OMLER_longitude, &
         OMLER_latitude, OMLER_wvl
    INTEGER (KIND=i2), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
         OMLER_monthly_albedo
    REAL (KIND=r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
         OMLER_albedo, OMLER_wvl_albedo

    ! --------------------
    ! More Local variables
    ! --------------------
    CHARACTER (LEN=34), PARAMETER :: grid_name = &
         'EarthSurfaceReflectanceClimatology'
    CHARACTER (LEN=maxchlen) :: grid_file

    INTEGER (KIND=i4), DIMENSION(1) :: minwvl, maxwvl, minlon,    &
         maxlon, minlat, maxlat
    INTEGER (KIND=i4), PARAMETER :: OMLER_n_latitudes   = 360, &
         OMLER_n_longitudes  = 720, &
         OMLER_n_wavelenghts =  23, &
         one                 =   1
    INTEGER (KIND=i4) :: ixtrack, itimes, ilon, ilat, nlon, nlat, OMnwvl
    INTEGER (KIND=i4) :: grid_id, grid_file_id, month


    REAL (KIND=r8), DIMENSION(1) :: plon, plat, midwvl
    REAL (KIND=r4) :: scale_factor, offset
    REAL (KIND=r8) :: lonp, latp

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: locerrstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=16), PARAMETER :: modulename = 'omi_omler_albedo' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok

    grid_file = TRIM(ADJUSTL(pcfvar%albedo_fname))

    ! -------------------------------------------------------------------------------
    ! Open he5 OMI OMLER grid file and check GRID_FILE_ID (-1 if error)
    ! -------------------------------------------------------------------------------
    grid_file_id = HE5_GDOPEN (grid_file, he5f_acc_rdonly)
    IF (grid_file_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5GDOPEN, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF

    ! ----------------------------------------------
    ! Attach to grid and check GRID_ID (-1 if error)
    ! ----------------------------------------------
    grid_id = HE5_GDattach (grid_file_id, grid_name)
    IF (grid_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5GDATTACH, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX( errstat, locerrstat)
    END IF

    ! -------------------------
    ! Read longitude data field
    ! -------------------------
    ALLOCATE (OMLER_longitude(OMLER_n_longitudes))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_longitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Longitude", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_longitude)
    errstat = MAX (errstat, locerrstat)

    ! ------------------------
    ! Read latitude data field
    ! ------------------------
    ALLOCATE (OMLER_latitude(OMLER_n_latitudes))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_latitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Latitude", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_latitude)
    errstat = MAX (errstat, locerrstat)

    ! ---------------------------
    ! Read wavelenghts data field
    ! ---------------------------
    ALLOCATE (OMLER_wvl(OMLER_n_wavelenghts))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_wavelenghts
    locerrstat = HE5_GDRDFLD(grid_id, "Wavelength", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_wvl)
    errstat = MAX (errstat, locerrstat)

    ! -----------------------------------------------
    ! Select the wavelenghts; finding array positions
    ! -----------------------------------------------
    midwvl = (ctrvar%amf_wvl + ctrvar%amf_wvl2) / 2.0_r8
    minwvl = MINLOC(OMLER_wvl, OMLER_wvl .GE. REAL(ctrvar%winwav_min,KIND=r4))
    maxwvl = MAXLOC(OMLER_wvl, OMLER_wvl .LE. REAL(ctrvar%winwav_max,KIND=r4))
    OMnwvl = maxwvl(1)-minwvl(1)+1

    ! ---------------------------
    ! Read the albedo data field:
    ! -Month
    ! -Selected wavelenghts
    ! ---------------------------
    ALLOCATE (OMLER_monthly_albedo(OMLER_n_longitudes, &
              OMLER_n_latitudes, OMnwvl,1))
    ALLOCATE (OMLER_wvl_albedo(OMLER_n_longitudes,     &
              OMLER_n_latitudes, OMnwvl,1))
    ALLOCATE (OMLER_albedo(OMLER_n_longitudes,         &
              OMLER_n_latitudes, 1,1))

    month = granule_month - 1

    he5_start_4d  = (/  0,  0, minwvl-1, month/)
    he5_stride_4d = (/  1,  1,        1,     1/)
    he5_edge_4d   = (/OMLER_n_longitudes,OMLER_n_latitudes, &
                      OMnwvl, 1/)
    locerrstat = HE5_GDRDFLD(grid_id, "MonthlySurfaceReflectance",  &
      he5_start_4d, he5_stride_4d, he5_edge_4d, OMLER_monthly_albedo)
    errstat = MAX (errstat, locerrstat)

    ! ----------------------------
    ! Read the albedo scale factor
    ! ----------------------------
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlySurfaceReflectance", &
         "ScaleFactor", scale_factor)
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlySurfaceReflectance", &
         "Offset", offset)

    ! --------------------
    ! Deattached from grid
    ! --------------------
    locerrstat = HE5_GDDETACH(grid_id)

    ! -------------------------------------------
    ! Close he5 OMI OMLER grid file (-1 if error)
    ! -------------------------------------------
    locerrstat = HE5_GDclose ( grid_file_id)
    IF ( locerrstat == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_W_HE5GDCLOSE, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF
  
    OMLER_wvl_albedo = REAL(offset, KIND = r8) +          &
                      REAL(scale_factor, KIND = r8)*      &
                      REAL(OMLER_monthly_albedo, KIND=r8)
 
    ! ----------------------------------------------------------------
    ! Interpolate for each pixel to one single wavelenght, just at the
    ! mid point of the fitting window: midwvl.
    ! ----------------------------------------------------------------
    DO ilon = 1, OMLER_n_longitudes
       DO ilat = 1, OMLER_n_latitudes
          CALL ezspline_1d_interpolation (                            &
               OMnwvl, REAL(OMLER_wvl(minwvl(1):maxwvl(1)), KIND=r8), &
               OMLER_wvl_albedo(ilon,ilat,1:OMnwvl,1),                &
               one, midwvl, OMLER_albedo(ilon,ilat,1,1), errstat )
       END DO
    END DO

    ! --------------------------------------------------
    ! Interpolate to the lat and longitude of each pixel
    ! --------------------------------------------------
    DO itimes = 0, nt-1

       DO ixtrack = 1, nx
          
          plon(1) = r8_missval
          plat(1) = r8_missval
          plon(1) = REAL(lon(ixtrack,itimes), KIND=r8)
          plat(1) = REAL(lat(ixtrack,itimes), KIND=r8)
          lonp = plon(1)
          latp = plat(1)

          IF (lon(ixtrack,itimes) .LT. -180.0_r4 .OR. &
              lat(ixtrack,itimes) .LT.  -90.0_r4 .OR. &
              lon(ixtrack,itimes) .GT.  180.0_r4 .OR. &
              lat(ixtrack,itimes) .GT.  90.0_r4) CYCLE

          ! -------------------------------------------------
          ! To speed up the interpolation a field of +- 2 deg
          ! rees on latitude and longitude will be used from
          ! OMLER centered around the OMI pixel.
          ! ----------------------------------------------------------
          ! Finding which OMLER pixels go in to the interpolation
          ! If we are too close to the boundaries I just move them to
          ! the interior. Good for the latitude, no so good for the 
          ! longitude, but over the Pacific it should make no much dif
          ! ference
          ! ----------------------------------------------------------
          IF (latp .LT. -88.0) latp = -89.0
          IF (latp .GT.  88.0) latp =  89.0
          IF (lonp .LT. -178.0) lonp = -179.0
          IF (lonp .GT.  178.0) lonp =  179.0

          minlon = MINLOC(OMLER_longitude, OMLER_longitude .GE. lonp-1.0)
          maxlon = MAXLOC(OMLER_longitude, OMLER_longitude .LE. lonp+1.0)
          minlat = MINLOC(OMLER_latitude,  OMLER_latitude  .GE. latp-1.0)
          maxlat = MAXLOC(OMLER_latitude,  OMLER_latitude  .LE. latp+1.0)
          nlon   = maxlon(1)-minlon(1)+1
          nlat   = maxlat(1)-minlat(1)+1                

          CALL ezspline_2d_interpolation ( nlon, nlat,             &
               REAL(OMLER_longitude(minlon(1):maxlon(1)),KIND=r8), &
               REAL(OMLER_latitude (minlat(1):maxlat(1)),KIND=r8), &
               OMLER_albedo(minlon(1):maxlon(1),                   &
                            minlat(1):maxlat(1), 1, 1),            &
               one, one, plon, plat, albedo(ixtrack,itimes),       &
               locerrstat )
       END DO
    END DO
    
    ! --------------------
    ! Deallocate variables
    ! --------------------
    DEALLOCATE (OMLER_monthly_albedo)
    DEALLOCATE (OMLER_wvl_albedo)
    DEALLOCATE (OMLER_albedo)
    DEALLOCATE (OMLER_longitude)
    DEALLOCATE (OMLER_latitude)
    DEALLOCATE (OMLER_wvl)
   
    errstat = MAX(errstat, locerrstat)
    
  END SUBROUTINE omi_omler_albedo
    
  SUBROUTINE extract_swathname ( nswath, multi_swath, swathstr, single_swath )

    ! ---------------------------------------------------------------------
    ! Extracts SINGLE_SWATH from MULTI_SWATH, based on presence of SWATHSTR
    ! ---------------------------------------------------------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nswath
    CHARACTER (LEN=*), INTENT (IN) :: multi_swath, swathstr

    ! ----------------
    ! Output variables
    ! ----------------
    CHARACTER (LEN=*), INTENT (OUT) :: single_swath

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4), DIMENSION(0:nswath) :: swsep
    INTEGER (KIND=i4) :: mslen, k, j1, j2, nsep
    CHARACTER (LEN=LEN(multi_swath)) :: tmpstr


    ! --------------------------
    ! Initialize output variable
    ! --------------------------
    single_swath = ''

    ! ---------------------------------
    ! Find length of MULTI_SWATH string
    ! ---------------------------------
    tmpstr = TRIM(ADJUSTL(multi_swath))
    mslen  = LEN_TRIM(ADJUSTL(tmpstr))

    ! --------------------------------------
    ! First find the number of "," in TMPSTR
    ! --------------------------------------
    nsep = 0 ; swsep(0:nswath) = 0
    DO k = 1, mslen
       IF ( tmpstr(k:k) == ',' ) THEN
          nsep = nsep + 1
          swsep(nsep) = k
       END IF
    END DO
    IF ( nsep == nswath-1 ) THEN
       nsep = nsep + 1 ; swsep(nsep) = mslen + 1
    END IF


    ! ----------------------------------------------------------------------
    ! Hangle along the positions of separators (commas, ",") between the
    ! concatinated Swath Name entries. The first Swath Name to contain 
    ! SWATHSTR is taken as the match - not a perfect rationale but simple
    ! enough if we have set up the AMF table file correctly.
    ! ----------------------------------------------------------------------    
    getswath: DO k = 1, nswath
       j1 = swsep(k-1)+1  ;  j2 = swsep(k)-1
       IF ( INDEX ( tmpstr(j1:j2), TRIM(ADJUSTL(swathstr)) ) > 0 ) THEN
          single_swath = TRIM(ADJUSTL(tmpstr(j1:j2)))
          EXIT getswath
       END IF
    END DO getswath

    RETURN
  END SUBROUTINE extract_swathname

  SUBROUTINE climatology_getdim ( &
       swath_id, Cmlat, Cmlon, CmETA, CmEp1, errstat )

    ! --------------------------------
    ! Return dimensions of Climatology
    ! --------------------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    INTEGER (KIND=i4), INTENT (OUT) :: Cmlat, Cmlon, CmETA, CmEp1

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4), PARAMETER :: maxdim = 100
    INTEGER (KIND=i4) :: ndim, nsep
    INTEGER (KIND=C_LONG) :: ndimcl
    INTEGER (KIND=i4) :: i, j, swlen, iend, istart
    INTEGER (KIND=i4), DIMENSION(0:maxdim) :: dim_array, dim_seps
    INTEGER (KIND=C_LONG), DIMENSION(0:maxdim) :: dim_arraycl
    CHARACTER (LEN=10*maxdim) :: dim_chars

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    Cmlat = -1 ; Cmlon = -1 ; CmETA = -1; CmEp1 = -1

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:maxdim) )
    ndim   = INT ( ndimcl, KIND=i4 )
    IF ( ndim <= 0 ) THEN
       errstat = MAX ( errstat, pge_errstat_error )
       RETURN
    END IF
    dim_array(0:maxdim) = INT ( dim_arraycl(0:maxdim), KIND=i4 )

    dim_chars =     TRIM(ADJUSTL(dim_chars))
    swlen     = LEN_TRIM(ADJUSTL(dim_chars))

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 ; dim_seps(0:ndim) = 0
    getseps: DO i = 1, swlen 
       IF ( dim_chars(i:i) == ',' ) THEN
          nsep           = nsep + 1
          dim_seps(nsep) = i
       END IF
    END DO getseps
    nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO j = 0, nsep-1
       istart = dim_seps(j)+1 ; iend = dim_seps(j+1)-1

       SELECT CASE ( dim_chars(istart:iend) )
       CASE ( "nLat" )
          Cmlat = dim_array(j)
       CASE ( "nLon" )
          Cmlon = dim_array(j)
       CASE ( "nETA" )
          CmETA = dim_array(j)
       CASE ( "nEp1")
          CmEp1 = dim_array(j)
       CASE DEFAULT
          ! Whatever. Nothing to be done here.
       END SELECT
    END DO getdims

    RETURN
  END SUBROUTINE climatology_getdim

  SUBROUTINE climatology_allocate ( ad, Cmlat, Cmlon, CmETA, CmEp1, errstat )

    USE OMSAO_casestring_module, ONLY: lower_case
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=1), INTENT (IN) :: ad
    INTEGER (KIND=i4), INTENT (IN) :: Cmlat, Cmlon, CmETA, CmEp1

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: estat
    CHARACTER (LEN=1) :: adlow


    estat = pge_errstat_ok

    ! ------------------------------------------
    ! Make sure AD ("a" or "d") is in lower case
    ! ------------------------------------------
    adlow = lower_case ( ad )

    SELECT CASE ( adlow )
    CASE ('a')
       ALLOCATE (latvals (Cmlat),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (lonvals (Cmlon),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Ap      (CmEp1),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Bp      (CmEp1),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Temperature (Cmlon, Cmlat, CmETA), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Gas_profiles(Cmlon, Cmlat, CmETA), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Psurface(Cmlon,Cmlat),            STAT=estat ) ; errstat = MAX ( errstat, estat )
    CASE ('d')
       IF ( ALLOCATED ( latvals      ) )  DEALLOCATE ( latvals      )
       IF ( ALLOCATED ( lonvals      ) )  DEALLOCATE ( lonvals      )
       IF ( ALLOCATED ( Ap           ) )  DEALLOCATE ( Ap           )
       IF ( ALLOCATED ( Bp           ) )  DEALLOCATE ( Bp           )
       IF ( ALLOCATED ( Temperature  ) )  DEALLOCATE ( Temperature  )
       IF ( ALLOCATED ( Gas_profiles ) )  DEALLOCATE ( Gas_profiles )
       IF ( ALLOCATED ( Psurface     ) )  DEALLOCATE ( Psurface     )
    CASE DEFAULT
       ! Whatever. Nothing to be done here
    END SELECT

    RETURN
  END SUBROUTINE climatology_allocate

  SUBROUTINE vlidort_allocate ( ad, anozo, ancld, ansza, anvza, anwav, analt, errstat )

    USE OMSAO_casestring_module, ONLY: lower_case
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=1), INTENT (IN) :: ad
    INTEGER (KIND=i4), INTENT (IN) :: anozo, ancld, ansza, anvza, anwav, analt

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: estat
    CHARACTER (LEN=1) :: adlow


    estat = pge_errstat_ok

    ! ------------------------------------------
    ! Make sure AD ("a" or "d") is in lower case
    ! ------------------------------------------
    adlow = lower_case ( ad )

    SELECT CASE ( adlow )
    CASE ('a')
       ALLOCATE (vl_OzC0(anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_OzC1(anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_OzC2(anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )

       ALLOCATE (vl_pre(ancld), STAT=estat )  ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_sza(ansza), STAT=estat )  ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_vza(anvza), STAT=estat )  ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_wav(anwav), STAT=estat )  ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_toms(anozo), STAT=estat ) ; errstat = MAX ( errstat, estat )

       ALLOCATE (vl_air(anozo,ancld,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_alt(anozo,ancld,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_ozo(anozo,ancld,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_tem(anozo,ancld,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )

       ALLOCATE (vl_I0(anozo,ancld,ansza,anvza,anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_I1(anozo,ancld,ansza,anvza,anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_I2(anozo,ancld,ansza,anvza,anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_Ir(anozo,ancld,ansza,anvza,anwav), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_Sb(anozo,ancld,anwav), STAT=estat )             ; errstat = MAX ( errstat, estat )

       ALLOCATE (vl_dI0(anozo,ancld,ansza,anvza,anwav,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_dI1(anozo,ancld,ansza,anvza,anwav,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_dI2(anozo,ancld,ansza,anvza,anwav,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vl_dIr(anozo,ancld,ansza,anvza,anwav,analt), STAT=estat ) ; errstat = MAX ( errstat, estat )

    CASE ('d')
       IF ( ALLOCATED ( vl_OzC0 ) )  DEALLOCATE ( vl_OzC0 )
       IF ( ALLOCATED ( vl_OzC1 ) )  DEALLOCATE ( vl_OzC1 )
       IF ( ALLOCATED ( vl_OzC2 ) )  DEALLOCATE ( vl_OzC2 )

       IF ( ALLOCATED ( vl_pre  ) )  DEALLOCATE ( vl_pre  )
       IF ( ALLOCATED ( vl_sza  ) )  DEALLOCATE ( vl_sza  )
       IF ( ALLOCATED ( vl_vza  ) )  DEALLOCATE ( vl_vza  )
       IF ( ALLOCATED ( vl_wav  ) )  DEALLOCATE ( vl_wav  )
       IF ( ALLOCATED ( vl_toms ) )  DEALLOCATE ( vl_toms )

       IF ( ALLOCATED ( vl_air ) )  DEALLOCATE ( vl_air )
       IF ( ALLOCATED ( vl_alt ) )  DEALLOCATE ( vl_alt )
       IF ( ALLOCATED ( vl_ozo ) )  DEALLOCATE ( vl_ozo )
       IF ( ALLOCATED ( vl_tem ) )  DEALLOCATE ( vl_tem )

       IF ( ALLOCATED ( vl_I0 ) )  DEALLOCATE ( vl_I0 )
       IF ( ALLOCATED ( vl_I1 ) )  DEALLOCATE ( vl_I1 )
       IF ( ALLOCATED ( vl_I2 ) )  DEALLOCATE ( vl_I2 )
       IF ( ALLOCATED ( vl_Ir ) )  DEALLOCATE ( vl_Ir )
       IF ( ALLOCATED ( vl_Sb ) )  DEALLOCATE ( vl_Sb )

       IF ( ALLOCATED ( vl_dI0 ) )  DEALLOCATE ( vl_dI0 )
       IF ( ALLOCATED ( vl_dI1 ) )  DEALLOCATE ( vl_dI1 )
       IF ( ALLOCATED ( vl_dI2 ) )  DEALLOCATE ( vl_dI2 )
       IF ( ALLOCATED ( vl_dIr ) )  DEALLOCATE ( vl_dIr )

    CASE DEFAULT
       ! Whatever. Nothing to be done here
    END SELECT

    RETURN
  END SUBROUTINE vlidort_allocate

  SUBROUTINE compute_geometric_amf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

    USE OMSAO_parameters_module, ONLY: deg2rad
    USE OMSAO_data_module, ONLY: geo_amf, oobview_amf

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nx, nt
    REAL (KIND=r4), DIMENSION (nx,0:nt-1), INTENT (IN) :: sza, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfgeo
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: it, spix, epix


    ! ---------------------------------------------
    ! Compute geometric AMF and set diagnostic flag
    ! ---------------------------------------------

    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with OMI_OOBVIEW_AMF.
    ! ----------------------------------------------------------------
    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)
       WHERE ( &
            sza(spix:epix,it) /= r4_missval .AND. &
            sza(spix:epix,it) >=     0.0_r4 .AND. &
            sza(spix:epix,it) <     90.0_r4 .AND. &
            vza(spix:epix,it) /= r4_missval .AND. &
            vza(spix:epix,it) >=     0.0_r4 .AND. &
            vza(spix:epix,it) <     90.0_r4         )
          amfgeo(spix:epix,it) = &
               1.0_r8 / COS ( REAL(sza(spix:epix,it),KIND=r8)*deg2rad ) + &
               1.0_r8 / COS ( REAL(vza(spix:epix,it),KIND=r8)*deg2rad )
          amfdiag(spix:epix,it) = geo_amf
       ELSEWHERE
          amfdiag(spix:epix,it) = oobview_amf
       ENDWHERE
    END DO

    RETURN
  END SUBROUTINE compute_geometric_amf

  SUBROUTINE read_vlidort (errstat)

    ! ====================================================
    ! This subroutine reads in the VLIDORT calculations to
    ! compute the Scattering Weights.
    ! It should check if the fitting window is included in
    ! the file, if not a warning should be printed and all
    ! the AMF diagnostic set to non computed.
    ! ====================================================
    IMPLICIT NONE    

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: hdferr

    INTEGER(HID_T) :: input_file_id                                  ! File identifier
    INTEGER(HID_T) :: OzC0_did, OzC1_did, OzC2_did,                & ! Dataset identifiers
                      llp_did, sza_did, toz_did, vza_did, wav_did, &
                      I0_did, I1_did, I2_did, Ir_did, Sb_did,      &
                      dI0_did, dI1_did, dI2_did, dIr_did, Fac_did, & 
                      air_did, alt_did, ozo_did, tem_did, dspace,  &
                      Tozo_datatype_id

    INTEGER(HSIZE_T), DIMENSION(1) :: hllp_dim, hllp_maxdim, &
                                      hsza_dim, hsza_maxdim, &
                                      htoz_dim, htoz_maxdim, &
                                      hvza_dim, hvza_maxdim, &
                                      hwav_dim, hwav_maxdim, &
                                      hfac_dim, hfac_maxdim
    INTEGER(HSIZE_T), DIMENSION(3) :: halt_dim, halt_maxdim, &
                                      hSb_dim,  hSb_maxdim
    INTEGER(HSIZE_T), DIMENSION(5) :: hI0_dim,  hI0_maxdim
    INTEGER(HSIZE_T), DIMENSION(6) :: hdI0_dim, hdI0_maxdim

    INTEGER(SIZE_T)                :: size
    LOGICAL, SAVE :: h5inited = .FALSE.

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    errstat = pge_errstat_ok

    ! ---------------------------------
    ! Initialize hdf5 FORTRAN Interface
    ! ---------------------------------
    if (.NOT.h5inited) then
      CALL h5open_f(hdferr)
      h5inited = .TRUE.
    endif 

    ! ------------------
    ! Dataset data types
    ! ------------------
    size = 4
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, Tozo_datatype_id, hdferr)
    CALL h5tset_size_f(Tozo_datatype_id, size, hdferr)

    ! ******************************************
    ! Find out the dimensions of the input file:
    !  # of pressure levels
    !  # of SZA
    !  # of Ozone profiles
    !  # of VZA
    !  # of wavelenghts
    !  # of altitude levels  
    ! ******************************************
    ! -------------------
    ! Opening input TABLE
    ! -------------------
    CALL h5fopen_f(TRIM(ADJUSTL(pcfvar%amf_table_fname)), H5F_ACC_RDONLY_F, &
                   input_file_id, hdferr)
    IF (hdferr .eq. -1) THEN
       WRITE(*,100) 'ERROR: Opening '//TRIM(ADJUSTL(pcfvar%amf_table_fname))
    END IF

    ! --------------------------------------------------------------------------
    ! Open ozone cross sections, grid, intensity, jacobians and profile datasets
    ! --------------------------------------------------------------------------
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C0', OzC0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C1', OzC1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C2', OzC2_did, hdferr)

    CALL h5dopen_f(input_file_id,'/Grid/Lower level pressure', llp_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/SZA', sza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/TOMS ozone', toz_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/VZA', vza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/Wavelength', wav_did,hdferr)
    
    CALL h5dopen_f(input_file_id,'/Intensity/I0', I0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/I1', I1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/I2', I2_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Ir', Ir_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Sb', Sb_did, hdferr)
    
    CALL h5dopen_f(input_file_id,'/Jacobians/Factor', Fac_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI0', dI0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI1', dI1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI2', dI2_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dIr', dIr_did, hdferr)
    
    CALL h5dopen_f(input_file_id,'/Profiles/Air profile', air_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Altitude', alt_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Ozone profile', ozo_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Temperature profile', tem_did, hdferr)
    
    ! -----------------------
    ! Find out the dimensions
    ! -----------------------
    CALL h5dget_space_f(llp_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hllp_dim, hllp_maxdim, hdferr)
    CALL h5dget_space_f(sza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hsza_dim, hsza_maxdim, hdferr)
    CALL h5dget_space_f(toz_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, htoz_dim, htoz_maxdim, hdferr)
    CALL h5dget_space_f(vza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hvza_dim, hvza_maxdim, hdferr)
    CALL h5dget_space_f(wav_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hwav_dim, hwav_maxdim, hdferr)
    CALL h5dget_space_f(I0_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hI0_dim, hI0_maxdim, hdferr)
    CALL h5dget_space_f(Sb_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hSb_dim, hSb_maxdim, hdferr)
    CALL h5dget_space_f(Fac_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hfac_dim, hfac_maxdim, hdferr)
    CALL h5dget_space_f(dI0_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hdI0_dim, hdI0_maxdim, hdferr)
    CALL h5dget_space_f(alt_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, halt_dim, halt_maxdim, hdferr)
    
    ! ---------------------------------------------------------------
    ! Allocate & initialize variables now that we have the dimensions
    ! ---------------------------------------------------------------
    vl_nozo = hdI0_dim(1)
    vl_ncld = hdI0_dim(2)
    vl_nsza = hdI0_dim(3)
    vl_nvza = hdI0_dim(4)
    vl_nwav = hdI0_dim(5)
    vl_nalt = hdI0_dim(6)

    CALL vlidort_allocate ("a", vl_nozo, vl_ncld, vl_nsza, vl_nvza, vl_nwav, vl_nalt, errstat)
  
    ! ----------------------------------------------------
    ! Read from the h5 file all these small size variables
    ! ----------------------------------------------------
    CALL h5dread_f(OzC0_did, H5T_NATIVE_REAL, vl_OzC0(1:vl_nwav), hwav_dim, hdferr)
    CALL h5dread_f(OzC1_did, H5T_NATIVE_REAL, vl_OzC1(1:vl_nwav), hwav_dim, hdferr)
    CALL h5dread_f(OzC2_did, H5T_NATIVE_REAL, vl_OzC2(1:vl_nwav), hwav_dim, hdferr)

    CALL h5dread_f(llp_did, H5T_NATIVE_REAL,  vl_pre(1:vl_ncld),  hllp_dim, hdferr)
    CALL h5dread_f(sza_did, H5T_NATIVE_REAL,  vl_sza(1:vl_nsza),  hsza_dim, hdferr)
    CALL h5dread_f(toz_did, Tozo_datatype_id, vl_toms(1:vl_nozo), htoz_dim, hdferr)
    CALL h5dread_f(vza_did, H5T_NATIVE_REAL,  vl_vza(1:vl_nvza),  hvza_dim, hdferr)
    CALL h5dread_f(wav_did, H5T_NATIVE_REAL,  vl_wav(1:vl_nwav),  hwav_dim, hdferr)
    
    CALL h5dread_f(I0_did, H5T_NATIVE_REAL,  &
         vl_I0(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav),  &
         hI0_dim, hdferr)
    CALL h5dread_f(I1_did, H5T_NATIVE_REAL,  &
         vl_I1(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav),  &
         hI0_dim, hdferr)
    CALL h5dread_f(I2_did, H5T_NATIVE_REAL,  &
         vl_I2(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav),  &
         hI0_dim, hdferr)
    CALL h5dread_f(Ir_did, H5T_NATIVE_REAL,  &
         vl_Ir(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav),  &
         hI0_dim, hdferr)
    CALL h5dread_f(Sb_did, H5T_NATIVE_REAL,  &
         vl_Sb(1:vl_nozo,1:vl_ncld,1:vl_nwav),  &
         hSb_dim, hdferr)
    
    CALL h5dread_f(Fac_did, H5T_NATIVE_REAL,  vl_Factor, hfac_dim, hdferr)
    CALL h5dread_f(dI0_did, H5T_NATIVE_REAL,  &
         vl_dI0(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav,1:vl_nalt), &
         hdI0_dim, hdferr)
    CALL h5dread_f(dI1_did, H5T_NATIVE_REAL,  &
         vl_dI1(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav,1:vl_nalt), &
         hdI0_dim, hdferr)
    CALL h5dread_f(dI2_did, H5T_NATIVE_REAL,  &
         vl_dI2(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav,1:vl_nalt), &
         hdI0_dim, hdferr)
    CALL h5dread_f(dIr_did, H5T_NATIVE_REAL,  &
         vl_dIr(1:vl_nozo,1:vl_ncld,1:vl_nsza,1:vl_nvza,1:vl_nwav,1:vl_nalt), &
         hdI0_dim, hdferr)
    
    CALL h5dread_f(air_did, H5T_NATIVE_REAL,  vl_air(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(alt_did, H5T_NATIVE_REAL,  vl_alt(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(ozo_did, H5T_NATIVE_REAL,  vl_ozo(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(tem_did, H5T_NATIVE_REAL,  vl_tem(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    
    ! --------------
    ! Close datasets
    ! --------------
    CALL h5dclose_f (OzC0_did, hdferr)
    CALL h5dclose_f (OzC1_did, hdferr)
    CALL h5dclose_f (OzC2_did, hdferr)
    
    CALL h5dclose_f (llp_did, hdferr)
    CALL h5dclose_f (sza_did, hdferr)
    CALL h5dclose_f (toz_did, hdferr)
    CALL h5dclose_f (vza_did, hdferr)
    CALL h5dclose_f (wav_did, hdferr)
    
    CALL h5dclose_f(I0_did, hdferr)
    CALL h5dclose_f(I1_did, hdferr)
    CALL h5dclose_f(I2_did, hdferr)
    CALL h5dclose_f(Ir_did, hdferr)
    CALL h5dclose_f(Sb_did, hdferr)
    
    CALL h5dclose_f(Fac_did, hdferr)
    CALL h5dclose_f(dI0_did, hdferr)
    CALL h5dclose_f(dI1_did, hdferr)
    CALL h5dclose_f(dI2_did, hdferr)
    CALL h5dclose_f(dIr_did, hdferr)
    
    CALL h5dclose_f (air_did, hdferr)
    CALL h5dclose_f (alt_did, hdferr)
    CALL h5dclose_f (ozo_did, hdferr)
    CALL h5dclose_f (tem_did, hdferr)
    
    ! ----------
    ! Close file
    ! ----------
    CALL h5fclose_f(input_file_id, hdferr)

    errstat = hdferr
    
100 FORMAT (A)

  END SUBROUTINE read_vlidort

  SUBROUTINE amf_read_ompsclouds ( cloud_filename, nt, nx, l2cfr, l2ctp, errstat )

    USE ReadH5dataset, ONLY: H5ReadDataset

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    CHARACTER (LEN=*), INTENT (IN) :: cloud_filename

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: l2cfr, l2ctp

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    REAL (KIND=r4), PARAMETER :: missval_cfr=-9999.0, missval_ctp=-9999.0
    REAL (KIND=r4), DIMENSION(:,:), POINTER :: cfr => NULL()
    REAL (KIND=r4), DIMENSION(:,:), POINTER :: ctp => NULL()

    locerrstat = pge_errstat_ok

    ! --------
    ! Allocate
    ! --------
    IF (ASSOCIATED(cfr)) DEALLOCATE(cfr)
    IF (ASSOCIATED(ctp)) DEALLOCATE(ctp)

    ALLOCATE(cfr(nx,nt), STAT = errstat)
    ALLOCATE(ctp(nx,nt), STAT = errstat)

    ! ----------------------------------------------------------------
    ! Read cloud fraction and cloud top pressure, and check for error.
    ! Eventually we may read the cloud uncertainties also, but for the
    ! first version we stick with just the basic cloud products.
    ! ----------------------------------------------------------------
    CALL H5ReadDataset(cloud_filename, &
         "/SCIENCE_DATA/CloudFractionforO3", cfr)
    CALL H5ReadDataset(cloud_filename, &
         "/SCIENCE_DATA/CloudPressureforO3", ctp)

    ! -----------------------------------------------------------
    ! Assign the cloud fraction array used in the AMF calculation
    ! -----------------------------------------------------------
    l2cfr(1:nx,0:nt-1) = REAL(cfr(1:nx,1:nt), KIND=r8)

    ! ---------------------------------------------------------------
    ! Assign the cloud top pressure array used in the AMF calculation
    ! ---------------------------------------------------------------
    l2ctp(1:nx,0:nt-1) = REAL(ctp(1:nx,1:nt), KIND=r8)

    ! ---------------------
    ! De-allocate cfr & ctp
    ! ---------------------
    DEALLOCATE(cfr)
    DEALLOCATE(ctp)

    ! ------------------------------------------------
    ! Force the cloud parameters into physical bounds.
    ! But make sure not to remove MissingValues.
    ! ------------------------------------------------
    WHERE ( l2cfr > REAL(missval_cfr, KIND=r8) .AND. l2cfr < 0.0_r8 )
       l2cfr = 0.0_r8
    ENDWHERE
    WHERE ( l2cfr > 1.0_r8 )
       l2cfr = 1.0_r8
    ENDWHERE
    WHERE ( l2ctp > REAL(missval_ctp, KIND=r8) .AND. l2ctp < 0.0_r8 )
       l2ctp = 0.0_r8
    ENDWHERE

    ! ------------------------------------------------------
    ! Replace cloud missing values by SAO PGE missing values
    ! ------------------------------------------------------
    WHERE ( l2cfr <= REAL(missval_cfr, KIND=r8) )
       l2cfr = r8_missval
    ENDWHERE
    WHERE ( l2ctp <= REAL(missval_ctp, KIND=r8) )
       l2ctp = r8_missval
    ENDWHERE
    
    RETURN
  END SUBROUTINE amf_read_ompsclouds

  
  SUBROUTINE amf_diagnostic ( nt, nx, sza, vza, xtrange, &
       l2cfr, l2ctp, amfdiag )

    USE OMSAO_data_module, ONLY: oobview_amf, bigsza_amf
    
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange

    ! ----------------
    ! Modified variabe
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: l2cfr, l2ctp

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: it, spix, epix

    ! -------------------------------------------------------------------
    ! AMFDIAG has already been set to "geometric" AMF where SZA and VZA
    ! information was available, and "missing"/"out of bounds" otherwise.
    ! Here we need to check for cloud information as well as for snow 
    ! and glint.
    ! -------------------------------------------------------------------

    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with real diagnostic flags.
    ! ----------------------------------------------------------------

    ! -------------------------------------------------------------------
    ! If the file privided for Scattering weights has no information with
    ! in the fitting window, then nothing is to be done here.
    ! The geometric AMFs flag will remain in amfdiag and no further calcu
    ! lation will be performed inside the calculate_amf and calculate_sca
    ! subroutines.
    ! -------------------------------------------------------------------
    ! ----------------------------------------------------------
    ! Check that the value in amf_wvl is within the range of the
    ! suplied scattering weights file.
    ! ----------------------------------------------------------
    IF ( &
         (ctrvar%amf_wvl  .LT. MINVAL(vl_wav) ) .OR. &
         (ctrvar%amf_wvl2 .GT. MAXVAL(vl_wav) ) ) THEN
       RETURN
    END IF

    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)

       ! ----------------
       ! Missing SZA, VZA
       ! ----------------
       WHERE (                                   &
            sza(spix:epix,it) <= r8_missval .OR. &
            vza(spix:epix,it) <= r8_missval        )
          amfdiag(spix:epix,it) = i2_missval
       END WHERE

       ! ----------------------------------------
       ! Out-of-Bound SZA, VZA (but not missing!)
       ! ----------------------------------------
       WHERE ( &
            ( sza(spix:epix,it) < MINVAL(vl_sza) ) .OR. &
            ( sza(spix:epix,it) > MAXVAL(vl_sza) ) .OR. &
            ( vza(spix:epix,it) < MINVAL(vl_vza) ) .OR. &
            ( vza(spix:epix,it) > MAXVAL(vl_vza) )      )
          amfdiag(spix:epix,it) = oobview_amf
       END WHERE

       ! -------------------------------------------------
       ! Out of bounds clouds (to high over land), make it
       ! the highest possible value in the look up table.
       ! Ask Xiong...
       ! -------------------------------------------------
       WHERE ( &
            ( l2ctp(spix:epix,it) < MINVAL(vl_pre(1:vl_ncld))*1013.0_r8 ) )
              l2ctp(spix:epix,it) = MINVAL(vl_pre(1:vl_ncld))*1013.0_r8
       END WHERE

       ! -------------------------------------------------
       ! Out of bounds clouds (to low over land), make it
       ! the highest possible value in the look up table.
       ! Ask Xiong...
       ! -------------------------------------------------
       WHERE ( &
            ( l2ctp(spix:epix,it) > MAXVAL(vl_pre(1:vl_ncld))*1013.0_r8 ) )
              l2ctp(spix:epix,it) = MAXVAL(vl_pre(1:vl_ncld))*1013.0_r8
       END WHERE

       ! ------------------------------------------------------
       ! For pixel without cloud information set amf to missval
       ! and flag to missval
       ! ------------------------------------------------------
       WHERE ( &
            ( l2cfr(spix:epix,it) .EQ. r8_missval ) .OR. &
            ( l2ctp(spix:epix,it) .EQ. r8_missval)       )
              amfdiag(spix:epix,it) = i2_missval
       END WHERE

       ! ------------------------------------------------------
       ! And AMFDIAG values > OOB must be good and are set to 0
       ! if we have "good" clouds
       ! ------------------------------------------------------
       WHERE ( &
            ( amfdiag(spix:epix,it) > oobview_amf ) .AND. &
            (l2cfr(spix:epix,it) >= 0.0_r8            ) .AND. &
            (l2ctp(spix:epix,it) >= 0.0_r8            )       )
             amfdiag(spix:epix,it) = 0_i2
       END WHERE

       ! --------------------------------------------------
       ! Angles above the top value set on the control file
       ! are calculated "using this maximum value".
       ! --------------------------------------------------
       WHERE ( &
            ( sza(spix:epix,it)     .GE. ctrvar%amf_max_sza) .AND. &
            ( amfdiag(spix:epix,it) .GT. oobview_amf ) )
              amfdiag(spix:epix,it) = bigsza_amf + amfdiag(spix:epix,it)
       END WHERE

    END DO

    RETURN
  END SUBROUTINE amf_diagnostic

  SUBROUTINE compute_scatt ( nt, nx, albedo, sza, vza, l2ctp, l2cfr, terrain_height, cli_heights, amfdiag, &
                             scattw)

    USE OMSAO_lininterpolation_module, ONLY: lininterpol
    USE OMSAO_variables_module, ONLY: pcfvar

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza, terrain_height
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: albedo, l2ctp, l2cfr
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT (IN) :: cli_heights
    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT (INOUT) :: scattw

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itime, ixtrack, ispre, isozo, isalt, iswav, issza, isvza, status, one
    INTEGER (KIND=i4), DIMENSION(1) :: iwavs, iwavf, index_thg, index_cld
    REAL (KIND=r8) :: temp, tempsquare, grad
    REAL (KIND=r8) :: ozo_abs, Intensity, Jacobian, Oz_xs, Intensity_cld, Jacobian_cld
    REAL (KIND=r8) :: crf, nwavs
    REAL (KIND=r8), DIMENSION(vl_ncld, vl_nsza, vl_nvza, vl_nalt) :: scattwe, scattwe_cld
    REAL (KIND=r8), DIMENSION(vl_ncld, vl_nsza, vl_nvza) :: Inte_clear, Inte_cloud
    REAL (KIND=r8), DIMENSION(vl_nalt) :: re_alt
    REAL (KIND=r8), DIMENSION(vl_ncld) :: re_pre
    REAL (KIND=r8), DIMENSION(vl_nsza) :: re_sza
    REAL (KIND=r8), DIMENSION(vl_nvza) :: re_vza
    REAL (KIND=r8) :: local_alb, local_sza, local_vza, local_thg, local_cld, local_cfr
    REAL (KIND=r8), DIMENSION(1) :: ezlocal_sza, ezlocal_vza
    REAL (KIND=r8), DIMENSION(1,1,1) :: cloud_scattw, clear_scattw
    REAL (KIND=r8), DIMENSION(1,1) :: Icr, Icl
    REAL (KIND=r8), DIMENSION(CmETA) :: local_chg
    REAL (kind=8),  PARAMETER :: d2r = 3.141592653589793d0/180.0  !! JED fix
    ! -----------------------------------
    ! Find look up table wavelength index
    ! No interpolation, closest available
    ! is selected.
    ! -----------------------------------
    one = 1_i4
    iwavs = MINLOC(ABS(vl_wav - REAL(ctrvar%amf_wvl,  KIND = r4) ))
    iwavf = MINLOC(ABS(vl_wav - REAL(ctrvar%amf_wvl2, KIND = r4) ))
    nwavs = REAL(iwavf(1), KIND=r8) - REAL(iwavs(1), KIND=r8) + 1.0_r8

    ! --------------------------------------------------------
    ! Re-order the pressure, altitude, sza, and vza dimensions
    ! from (look up table). Needed for interpolation. 
    ! This should be moved to the program that creates the
    ! look up tables
    ! ------------------------------- ------------------------
    DO issza = 1, vl_nsza
       re_sza(issza) = cos(d2r*REAL(vl_sza(vl_nsza+1-issza), KIND = r8))  ! JED fix
    END DO
    DO isvza = 1, vl_nvza
       re_vza(isvza) = cos(d2r*REAL(vl_vza(vl_nvza+1-isvza), KIND = r8)) ! JED fix
    END DO
    DO ispre = 1, vl_ncld
       re_pre(ispre) = REAL(vl_pre(ispre), KIND = r8) * 1013.0_r8
    END DO
    DO isalt = 1, vl_nalt
       re_alt(isalt) = 1013.0_r8 * (10.0_r8 ** ( REAL(vl_alt(1,1,isalt), KIND = r8) / (-16.0_r8)))
    END DO

    ! ---------------
    ! Loop over lines
    ! ---------------
    DO itime = 0, nt-1
       ! --------------------------
       ! Loop over xtrack positions
       ! --------------------------
       DO ixtrack = 1, nx

          IF (amfdiag(ixtrack,itime) .LT. 0) CYCLE

          ! ----------------------------------------------
          ! If this point is reached then scattw should be
          ! different from r8_missval and it needs to be
          ! initialized to 0.0 to work out the average
          ! ----------------------------------------------
          scattw(ixtrack,itime,:) = 0.0_r8

          ! ----------------------------------------------
          ! If sza > amf_max_sza set it for calculation to
          ! amf_max_sza
          ! ----------------------------------------------
          local_sza = REAL(sza(ixtrack,itime), KIND = r8)
          IF (local_sza .GT. ctrvar%amf_max_sza) local_sza = ctrvar%amf_max_sza

          ! ---------------------
          ! Albedo for this pixel
          ! ---------------------
          local_alb    = albedo(ixtrack,itime)
          local_cld    = l2ctp(ixtrack,itime)
          local_cfr    = l2cfr(ixtrack,itime)
          local_sza    = cos(d2r*REAL(sza(ixtrack,itime), KIND = r8))  ! JED fix
          local_vza    = cos(d2r*REAL(vza(ixtrack,itime), KIND = r8))  ! JED fix
          local_thg    = REAL(terrain_height(ixtrack,itime), KIND = r8)
          local_chg(:) = cli_heights(ixtrack,itime,:)

          ! -------------------------------------------------------------------------
          ! If local_cfr is < 0.05 the pressure retrieved is effective cloud pressure 
          ! and we should consider the pixel as clear sky
          ! -------------------------------------------------------------------------
          IF (local_cfr .LT. 0.05) local_cfr = 0.0d0

          !Bringing it to the lowest available pressure if needed
          IF (local_thg .GT. 1013.0) local_thg = 1013.0_r8
          !Bringing clouds heights to lowest available pressure if needed. Weird yes, but just in case
          IF (local_cld .GT. 1013.0) local_cld = 1013.0_r8

          ! -----------------------------------------------
          ! Find the two closest cloud top levels from the
          ! scattering weights table suitable for local_thg
          ! and local_cld
          ! -----------------------------------------------
          index_thg = MINLOC(ABS(re_pre - local_thg))
          ! Be sure that we are below the level of the surface in the table
          IF ( (local_thg .GT. re_pre(index_thg(1))) .AND. (index_thg(1) .GT. 1) ) index_thg = index_thg-1
          index_cld = MINLOC(ABS(re_pre - local_cld))
          ! Be sure that we are below the level of the cloud in the table
          IF ( (local_cld .GT. re_pre(index_cld(1))) .AND. (index_cld(1) .GT. 1) ) index_cld = index_cld-1

          ! ----------------------------------------------------------
          ! First compute back from the parametrization the scattering
          ! weights for the given wavelength and albedo.
          ! ----------------------------------------------------------          
          DO isozo = 1, 1
             DO iswav = iwavs(1), iwavf(1) !vl_nwav

                ! ---------------------------
                ! Initialize scattwe to zeros
                ! ---------------------------
                scattwe      = 0.0_r8
                scattwe_cld  = 0.0_r8

                DO ispre = 1, vl_ncld 
                   DO issza = 1, vl_nsza
                      DO isvza = 1, vl_nvza

                         IF (ISNAN(vl_Ir(isozo,ispre,issza,isvza,iswav))) THEN
                            vl_Ir(isozo,ispre,issza,isvza,iswav) = 0.0
                         END IF

                         Intensity = &
                              REAL(vl_I0(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              REAL(vl_I1(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              REAL(vl_I2(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              ( (albedo(ixtrack, itime)                                  * &
                              REAL(vl_Ir(isozo,ispre,issza,isvza,iswav), KIND = r8) ) / &
                              (1.0_r8 - (albedo(ixtrack,itime)                        * &
                              REAL(vl_Sb(isozo,ispre,iswav), KIND = r8))) )
                         
                         Intensity_cld = &
                              REAL(vl_I0(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              REAL(vl_I1(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              REAL(vl_I2(isozo,ispre,issza,isvza,iswav), KIND = r8)      + &
                              ( (ctrvar%amf_alb_cld                                             * &
                              REAL(vl_Ir(isozo,ispre,issza,isvza,iswav), KIND = r8) ) / &
                              (1.0_r8 - (ctrvar%amf_alb_cld                                  * &
                              REAL(vl_Sb(isozo,ispre,iswav), KIND = r8))) )
                         
                         DO isalt = 1, vl_nalt


                            Temp       = REAL(vl_tem(isozo,ispre,isalt), KIND = r8) - 273.15_r8
                            TempSquare = Temp * Temp

                            Oz_xs = REAL(vl_OzC0(iswav), KIND = r8) + &
                                    REAL(vl_OzC1(iswav), KIND = r8) * Temp + &
                                    REAL(vl_OzC2(iswav), KIND = r8) * TempSquare

                            ozo_abs = Oz_xs * REAL(vl_ozo(isozo,ispre,isalt), KIND = r8)

                            Jacobian = ( &
                                 REAL(vl_dI0(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 REAL(vl_dI1(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 REAL(vl_dI2(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 ( (albedo(ixtrack,itime)                                      * &
                                 REAL(vl_dIr(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)) / &
                                 (1.0_r8 - (albedo(ixtrack,itime)                              * &
                                 REAL(vl_Sb(isozo,ispre,iswav), KIND = r8))) )  )                &
                                 / REAL(vl_Factor, KIND = r8)

                            Jacobian_cld = ( &
                                 REAL(vl_dI0(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 REAL(vl_dI1(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 REAL(vl_dI2(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)  + &
                                 ( (ctrvar%amf_alb_cld                                                * &
                                 REAL(vl_dIr(isozo,ispre,issza,isvza,iswav,isalt), KIND = r8)) / &
                                 (1.0_r8 - (ctrvar%amf_alb_cld                                 * &
                                 REAL(vl_Sb(isozo,ispre,iswav), KIND = r8))) )  )                &
                                 / REAL(vl_Factor, KIND = r8)

                            ! -----------------------------------------
                            ! vl_ncld+1-ispre & vl_nalt+1-isalt to have
                            ! ascending orden for the interpolation
                            ! ----------------------------------------------------
                            ! Intensities for the calculation of the cloudy pixels
                            ! ----------------------------------------------------
                            Inte_clear(ispre, vl_nsza+1-issza, vl_nvza+1-isvza) = Intensity
                            Inte_cloud(ispre, vl_nsza+1-issza, vl_nvza+1-isvza) = Intensity_cld

                            IF (Jacobian .NE. 0.0_r8 .AND. Intensity .NE. 0.0_r8 &
                                 .AND. ozo_abs .NE. 0.0_r8) THEN
                               scattwe(ispre, vl_nsza+1-issza, vl_nvza+1-isvza, isalt) = &
                                    -Jacobian / Intensity / ozo_abs
                            END IF
                            IF (Jacobian_cld .NE. 0.0_r8 .AND. Intensity_cld .NE. 0.0_r8 &
                                 .AND. ozo_abs .NE. 0.0_r8) THEN
                               scattwe_cld(ispre, vl_nsza+1-issza, vl_nvza+1-isvza, isalt) = &
                                    -Jacobian_cld / Intensity_cld / ozo_abs
                            END IF
                            
                         END DO ! End altitudes (scattering look up tables)
                      END DO ! End vza loop
                   END DO ! End sza loop
                END DO ! End pressure loop

                ! ---------------------------------------
                ! Working out the cloud radiance fraction
                ! See note below, Boersma et al. 2011 and
                ! Martin et al. 2003 (crf used below)
                ! ---------------------------------------
                ezlocal_sza = local_sza
                ezlocal_vza = local_vza
                 
                Icr = 0.0_r8
                CALL ezspline_2d_interpolation (vl_nsza,vl_nvza,re_sza,re_vza,Inte_clear(index_thg(1),:,:), &
                     one,one,ezlocal_sza(one),ezlocal_vza(one),Icr(one,one), &
                     status)
                Icl = 0.0_r8
                CALL ezspline_2d_interpolation (vl_nsza,vl_nvza,re_sza,re_vza,Inte_cloud(index_cld(1),:,:), &
                     one,one,ezlocal_sza(one),ezlocal_vza(one),Icl(one,one), &
                     status)
                crf = 0.0_r8
                crf = local_cfr * Icl(one,one) / &
                     (local_cfr * Icl(one,one) + (1 - local_cfr) * Icr(one,one) )
                
                ! ----------------------------------
                ! Interpolate for each altitude
                ! to the given sza, vza and pressure
                ! ----------------------------------
                DO isalt = 1, CmETA ! Loop over altitudes, climatology
                   
                   cloud_scattw(one,one,one) = 0.0_r8
                   clear_scattw(one,one,one) = 0.0_r8
                   ! --------------------------------------------------
                   ! If we are below the level of the clouds, the cloud
                   ! scattering weights are not needed. Cloud fraction
                   ! must be GT than 0.0. Only interpolate for values
                   ! above the cloud top.
                   ! --------------------------------------------------
                   IF ( local_chg(isalt) .LE. local_cld .AND. local_cfr .GT. 0.0 ) THEN
                      cloud_scattw(one,one,one) =  linInterpol (              &
                           vl_nsza, vl_nvza, vl_nalt,                         &
                           re_sza,  re_vza,  re_alt,                          &
                           scattwe_cld(index_cld(1),:,:,:),                   &
                           local_sza, local_vza, local_chg(isalt),            &
                           status=status)
                   END IF
                   
                   ! --------------------------------------------------
                   ! If we are below the level of the land, no need to
                   ! work out those scattering weights.
                   ! --------------------------------------------------
                   IF ( local_chg(isalt) .LE. local_thg ) THEN
                      clear_scattw(one,one,one) = linInterpol (           &
                           vl_nsza, vl_nvza, vl_nalt,                     &
                            re_sza,  re_vza,  re_alt,                     &
                           scattwe(index_thg(1),:,:,:),                   &
                           local_sza, local_vza, local_chg(isalt),        &
                           status=status)
                   END IF

                   ! ---------------------------------------------------------------------------------
                   ! Boersma et al. 2011 AMT, 4, 2011
                   ! Cloud radiance fraction: Crf= Cfr * Icl / Ir
                   !  We define Icl = From the Vlidort calculation, see above
                   !            Icr = From the Vlidort calculation, see above
                   !             Ir = Cfr * Icl + (1 - Cfr) * Icr (Total pixel radiance)
                   ! 
                   ! Now the scattering weights become w = crf * scatt_cloud + (1 - crf) * scatt_clear
                   !  We add the scattweights calculated in the previous wavelengths.
                   ! ---------------------------------------------------------------------------------
                   scattw(ixtrack,itime,isalt) = (crf * cloud_scattw(one,one,one) + (1.0_r8 - crf) * clear_scattw(one,one,one)) &
                                                 + scattw(ixtrack,itime,isalt)

                END DO ! End looop over altitudes
                ! -------------------------------------------------------------------------------------------------
                ! If there is a negative value on the lower layer of the scattw it has to do with the interpolation
                ! Quick fix, work out a gradient from two layers above and apply it to this first layer assuming
                ! linear behaviour. Next version of the lookup tables should have at least two levels below surface
                ! level to prevent the need for this. Even worst the 0.7 scaling for extreme cases
                ! -------------------------------------------------------------------------------------------------
                IF (scattw(ixtrack, itime, 1) .LT. 0.0) THEN
                   grad = (scattw(ixtrack, itime,3) - scattw(ixtrack, itime, 2)) / log(local_chg(3)/local_chg(2))
                   scattw(ixtrack, itime, 1) = scattw(ixtrack, itime, 2) - grad * log(local_chg(2)/local_chg(1))
                   IF (scattw(ixtrack,itime,1) .LT. 0.0) scattw(ixtrack,itime,1) = scattw(ixtrack,itime,2) * 0.7
                END IF

             END DO ! End wavelength loop

          END DO ! End ozone profile loop
          scattw(ixtrack,itime,1:CmETA) = scattw(ixtrack,itime,1:CmETA) / nwavs          

          !  Set non-physical entries to zero.
          WHERE ( scattw(ixtrack,itime,1:CmETA) < 0.0_r8 )
             scattw(ixtrack,itime,1:CmETA) = 0.0_r8
          END WHERE

       END DO ! End loop xtrack
       IF ( pcfvar%verb_thresh_lev .GE. vb_lev_screen ) WRITE(*,*) 'Scattering weights line', itime

    END DO ! End loop lines
    
  END SUBROUTINE COMPUTE_SCATT

  SUBROUTINE compute_amf ( nt, nx, CmETA, climatology, &
                           scattw, saoamf, amfdiag )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT(IN) :: nt, nx, CmETA
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT(IN) :: climatology, scattw
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT(IN) :: amfdiag

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saoamf

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: ixtrack, itimes

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    DO itimes  = 0, nt-1 ! Swath lines loop
       DO ixtrack = 1, nx ! Xtrack pixel loop

          ! ------------------------------------------------
          ! Set all amf with amfdiag < 0 to r8_missval
          ! ------------------------------------------------
          IF ( amfdiag(ixtrack,itimes) .LT. 0 ) THEN
             saoamf(ixtrack,itimes) = r8_missval
             CYCLE
          ENDIF

          ! -------------------------
          ! Finally work out the AMFs
          ! -------------------------
          saoamf(ixtrack,itimes) = SUM(scattw(ixtrack, itimes, 1:CmETA) * &
                                        climatology(ixtrack,itimes,1:CmETA))     / &
                                   SUM(climatology(ixtrack,itimes,1:CmETA))
 
       END DO ! Finish xtrack pixel loop
    END DO ! Finish 
    
  END SUBROUTINE compute_amf

  SUBROUTINE write_albedo_he5(albedo, nt, nx, errstat)

    ! ==================================================================
    ! This routines writes the albedos obtained from the OMLER climatolo
    ! gy to the output file.
    ! ==================================================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (IN) :: albedo

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

    locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 1 /)
    he5_edge_2d   = (/ nx, nt /)      


    colloc = albedo
    CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(albedo_field)),             &
                               he5_start_2d, he5_stride_2d, he5_edge_2d,&
                               colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE write_albedo_he5

  SUBROUTINE write_climatology_he5(climatology, cli_heights, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the Target Gas Profiles from the GEOS-Chem
    ! climatology to the output file.
    ! ===============================================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: climatology, cli_heights

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc

    locerrstat = pge_errstat_ok

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)      


    colloc = climatology
    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(gasprofile_field)),         &
                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
                               colloc(1:nx,0:nt-1,1:nl) )
    errstat = MAX ( errstat, locerrstat )

    colloc = cli_heights
    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(clialtgrid_field)),         &
                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
                               colloc(1:nx,0:nt-1,1:nl) )
    errstat = MAX ( errstat, locerrstat )
    
    
  END SUBROUTINE write_climatology_he5

  SUBROUTINE write_scatt_he5(scattw, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the scattering weigths to the output file.
    ! ===============================================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: scattw !, akernels

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                               :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc

    locerrstat = pge_errstat_ok

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)      


    colloc = scattw
    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(scaweights_field)),         &
                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
                               colloc(1:nx,0:nt-1,1:nl) )

    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE write_scatt_he5


SUBROUTINE he5_amf_write ( &
         pge_idx, nx, nt, saocol, saodco, amfmol, amfgeo, amfdiag, &
         amfcfr, amfctp, errstat )

  USE OMSAO_indices_module, ONLY: pge_hcho_idx, pge_gly_idx

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nt
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saocol, saodco
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfmol, amfgeo
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfcfr, amfctp
  INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat
  REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: amfloc
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

  
  ! -------------------------------------------------------------
  ! Air mass factor plus diagnostic.
  ! -------------------------------------------------------------
  ! As yet, only OMBRO and OMHCHO have true, non-geometric AMFs.
  ! But we try to have symmetric data fields as much as possible,
  ! hence the presence of the "molecule specific" AMF and its
  ! diagnostic for all PGEs. Non-OMBRO and -OMHCHO PGEs carry a
  ! geometric AMF here.
  ! 
  ! For completeness, the geometric AMF is added.
  ! -------------------------------------------------------------

  locerrstat = pge_errstat_ok

  he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

  ! ----------------------------------------
  ! (1) AMF diagnostic. No rounding required
  ! ----------------------------------------
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfdiag_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfdiag(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! -----------------
  ! (2) Geometric AMF
  ! -----------------
  amfloc = REAL ( amfgeo, KIND=r4 )
  CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfgeo_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! -----------------
  ! (3) Molecular AMF
  ! -----------------
  amfloc = REAL ( amfmol, KIND=r4 )
  CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! ----------------------------------------------------------
  ! (4) OMHCHO, OMCHOCHO only: AMF cloud fraction and pressure
  ! ----------------------------------------------------------
  IF ( pge_idx == pge_hcho_idx .OR. pge_idx == pge_gly_idx ) THEN
     amfloc = REAL ( amfcfr, KIND=r4 )
     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfcfr_field)), &
          he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
     errstat = MAX ( errstat, locerrstat )
  
     amfloc = REAL ( amfctp, KIND=r4 )
     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfctp_field)), &
          he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
     errstat = MAX ( errstat, locerrstat )
  END IF

  ! -----------------------------------------------------------------------
  ! (5) All PGEs: Output of columns and column uncertainties. For some PGEs
  !     (e.g., OMBRO, OMHCHO, OMCHOCHO) those have been adjusted by the AMF,
  !     but we have as yet to perform the rounding for any of them.
  ! -----------------------------------------------------------------------
  colloc = saocol
  CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )
  
  colloc = saodco
  CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE he5_amf_write

END MODULE OMSAO_wfamf_module

