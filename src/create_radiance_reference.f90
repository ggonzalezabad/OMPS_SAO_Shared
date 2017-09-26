SUBROUTINE create_radiance_reference (nt, nx, nw, locerrstat)

  USE OMSAO_precision_module
  USE OMSAO_errstat_module, ONLY : pge_errstat_ok
  USE OMSAO_radiance_ref_module
  USE OMSAO_omidata_module, ONLY: &
       omi_ccdpix_selection, omi_nwav_radref, omi_radref_spec, omi_radref_wavl,     &
       omi_radref_qflg, omi_radref_sza, omi_radref_vza, omi_radref_wght,            &
       omi_ccdpix_exclusion, omi_sol_wav_avg 
  USE OMSAO_variables_module, ONLY : ctrvar, pcfvar
  USE OMSAO_parameters_module
  USE OMSAO_omps_reader

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4) :: nt, nx, nw

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4) :: locerrstat

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE(omps_nmev_type) :: OMPS_data_radiance_reference
  INTEGER (KIND=i2)                          :: omps_reader_status_reference
  INTEGER (KIND=i2), DIMENSION(nw)           :: qflg_mask
  REAL    (KIND=r4)                          :: lat_midpt
  REAL    (KIND=r4), DIMENSION(nx,0:nt-1)    :: latr4, tmp_szenith, tmp_vzenith
  REAL    (KIND=r4), DIMENSION(nx)           :: szacount
  REAL    (KIND=r8)                          :: specsum
  REAL    (KIND=r8), DIMENSION(nw)           :: cntr8, radref_wavl_ix
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: radref_spec, radref_wavl
  REAL    (KIND=r8), DIMENSION(nx,nw)        :: allcount, dumcount
  REAL    (KIND=r8), DIMENSION(nw,nx,0:nt-1) :: tmp_radiance_spec, tmp_radiance_wavl
  INTEGER (KIND=i4), DIMENSION(0:nt-1,2)     :: xtrange
  INTEGER (KIND=i4)                          :: fpix, lpix, midpt_line, iline, ix, &
                                                icnt, imin, imax, j1
  LOGICAL                                    :: yn_have_scanline
  LOGICAL,           DIMENSION (2)           :: yn_have_limits

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = pge_errstat_ok

  ! ------------------------------------
  ! Read OMPS radiance reference granule
  ! ------------------------------------
  omps_reader_status_reference = OMPS_NMEV_READER(OMPS_data_radiance_reference, &
                                                    pcfvar%l1b_radref_fname)

  latr4(1:nx,0:nt-1)                  = OMPS_data_radiance_reference%Latitude(1:nx,1:nt)
  tmp_szenith(1:nx,0:nt-1)            = OMPS_data_radiance_reference%SolarZenithAngle(1:nx,1:nt)
  tmp_vzenith(1:nx,0:nt-1)            = OMPS_data_radiance_reference%SatelliteZenithAngle(1:nx,1:nt)
  tmp_radiance_spec(1:nw,1:nx,0:nt-1) = OMPS_data_radiance_reference%Radiance(1:nw,1:nx,1:nt)
  tmp_radiance_wavl(1:nw,1:nx,0:nt-1) = OMPS_data_radiance_reference%BandCenterWavelengths(1:nw,1:nx,1:nt)
  
  ! ------------------------------
  ! Initialize some some variables
  ! ------------------------------
  radiance_reference_lnums = -1  ! This will be written to file, hence needs a value
  lat_midpt = SUM ( ctrvar%radref_latrange ) / 2.0_r4

  CALL omi_set_xtrpix_range ( &
       nt, nx, ctrvar%pixnum_lim(3:4), &
       xtrange(0:nt-1,1:2), fpix, lpix, locerrstat )

  ! ----------------------------------------------------------------------
  ! Locate the swath line numbers corresponding the center of the latitude
  ! range to average into radiance reference spectrum.
  ! ----------------------------------------------------------------------
  CALL find_swathline_by_latitude ( &
       nx, 0, nt-1, latr4(1:nx,0:nt-1), lat_midpt, &
       xtrange(0:nt-1,1:2), midpt_line, yn_have_scanline )

  ! --------------------------------------------------------------------
  ! If lower and upper bounds of the radiance reference block to average
  ! are identical, then we keep the midpoint line number as the only
  ! reference. Else locate the corresponding swath line numbers.
  ! --------------------------------------------------------------------
  IF ( ctrvar%radref_latrange(1) == ctrvar%radref_latrange(2) ) THEN
     radiance_reference_lnums(1:2) = midpt_line
     yn_have_limits(1:2)           = .TRUE.
  ELSE
     CALL find_swathline_by_latitude ( &
          nx, 0, midpt_line, latr4(1:nx,0:midpt_line), ctrvar%radref_latrange(1), &
          xtrange(0:midpt_line,1:2), radiance_reference_lnums(1), yn_have_limits(1)   )
     CALL find_swathline_by_latitude ( &
          nx, midpt_line, nt-1, latr4(1:nx,midpt_line:nt-1), ctrvar%radref_latrange(2), &
          xtrange(midpt_line:nt-1,1:2), radiance_reference_lnums(2), yn_have_limits(2) )
  END IF

  ! --------------------------------------------------------------------
  ! Now we can average the spectra and the wavelength arrays. Loop over
  ! the block of swath lines in multiples of NLINES_MAX (100 by default)
  ! --------------------------------------------------------------------
  allcount       = 0.0_r8  ; dumcount       = 0.0_r8 ;  szacount = 0.0_r4
  radref_wavl    = 0.0_r8  ; radref_spec    = 0.0_r8
  omi_radref_sza = 0.0_r4  ; omi_radref_vza = 0.0_r4

  DO iline = radiance_reference_lnums(1), radiance_reference_lnums(2)

     ! ------------------------------------------------------
     ! Skip this cross-track position if there isn't any data
     ! ------------------------------------------------------
     fpix = xtrange(iline,1)
     lpix = xtrange(iline,2)

     DO ix = 1, nx
        
        IF ( (ix < fpix) .OR. (ix > lpix) ) CYCLE

        ! ----------------------------------------------------------------
        ! For now with the quality flags only check for zero values (good)
        ! ----------------------------------------------------------------
        qflg_mask(1:nw) = OMPS_data_radiance_reference%PixelQualityFlags(1:nw,nx,nt)
        cntr8(1:nw) = 1.0_r8
        WHERE ( qflg_mask(1:nw) > 0_i2 )
           cntr8(1:nw) = 0.0_r8
        END WHERE

        ! ------------------------------------
        ! Only proceed if we have a good value
        ! ------------------------------------
        IF ( ANY ( cntr8(1:nw) > 0.0_r8 ) ) THEN

           tmp_radiance_spec(1:nw,ix,iline) = tmp_radiance_spec(1:nw,ix,iline) * cntr8(1:nw)
           
           specsum = SUM ( tmp_radiance_spec(1:nw,ix,iline) ) / SUM ( cntr8(1:nw) )
           IF ( specsum == 0.0_r8 ) specsum = 1.0_r8
           
           specsum = 1.0_r8
           
           radref_spec(ix,1:nw) = &
                radref_spec(ix,1:nw) + tmp_radiance_spec(1:nw,ix,iline)/specsum
           radref_wavl(ix,1:nw) = &
                radref_wavl(ix,1:nw) + tmp_radiance_wavl(1:nw,ix,iline)
           allcount(ix,1:nw) = allcount(ix,1:nw) + cntr8(1:nw)
           dumcount(ix,1:nw) = dumcount(ix,1:nw) + 1.0_r8
           
           IF ( tmp_szenith(ix,iline) /= r4_missval .AND. &
                tmp_vzenith(ix,iline) /= r4_missval         ) THEN
              omi_radref_sza(ix) = omi_radref_sza(ix) + tmp_szenith(ix,iline)
              omi_radref_vza(ix) = omi_radref_vza(ix) + tmp_vzenith(ix,iline)
              szacount      (ix) = szacount(ix) + 1.0_r4
           END IF
           
        END IF
        
     END DO ! xtrack loop
     
  END DO ! line loop
  
  ! -----------------------------------------------------------
  ! Now for the actual averaging and assignment of final arrays
  ! -----------------------------------------------------------
  DO ix = 1, nx
     ! -----------------------------------
     ! Average the wavelengths and spectra
     ! -----------------------------------
     WHERE ( allcount(ix,1:nw) /= 0.0_r8 )
        radref_spec(ix,1:nw) = radref_spec(ix,1:nw) / allcount(ix,1:nw)
     END WHERE
     WHERE ( dumcount(ix,1:nw) /= 0.0_r8 )
        radref_wavl(ix,1:nw) = radref_wavl(ix,1:nw) / dumcount(ix,1:nw)
     END WHERE

     ! -------------------------------------------
     ! Average the Solar and Viewing Zenith Angles
     ! -------------------------------------------
     IF ( szacount(ix) > 0.0_r4 ) THEN
        omi_radref_sza(ix) = omi_radref_sza(ix) / szacount(ix)
        omi_radref_vza(ix) = omi_radref_vza(ix) / szacount(ix)
     ELSE
        omi_radref_sza(ix) = r4_missval
        omi_radref_vza(ix) = r4_missval             
     END IF

     ! -------------------------------------------------------------------------------
     ! Determine the CCD pixel numbers based on the selected wavelength fitting window
     ! -------------------------------------------------------------------------------
     radref_wavl_ix = radref_wavl(ix,1:nw)
     DO j1 = 1, 3, 2
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winwav_lim(j1  ),KIND=r8), 'LE', &
             omi_ccdpix_selection(ix,j1  ) )
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winwav_lim(j1+1),KIND=r8), 'GE', &
             omi_ccdpix_selection(ix,j1+1) )
     END DO

     imin = omi_ccdpix_selection(ix,1)
     imax = omi_ccdpix_selection(ix,4)

     icnt = imax - imin + 1
     omi_nwav_radref(       ix) = icnt
     omi_radref_spec(1:icnt,ix) = radref_spec(ix,imin:imax)
     omi_radref_wavl(1:icnt,ix) = radref_wavl_ix(imin:imax)
     omi_radref_qflg(1:icnt,ix) = 0_i2
     omi_radref_wght(1:icnt,ix) = normweight

     ! -----------------------------------------------------------------
     ! Re-assign the average solar wavelength variable, sinfe from here
     ! on we are concerned with radiances.
     ! -----------------------------------------------------------------
     omi_sol_wav_avg(ix) =  SUM( omi_radref_wavl(1:icnt,ix) ) / REAL(icnt, KIND=r8)

     ! ------------------------------------------------------------------
     ! Set weights and quality flags to "bad" for missing spectral points
     ! ------------------------------------------------------------------
     allcount(ix,1:icnt) = allcount(ix,imin:imax)
     WHERE ( allcount(ix,1:icnt) == 0.0_r8 )
        omi_radref_qflg(1:icnt,ix) = 7_i2
        omi_radref_wght(1:icnt,ix) = downweight
     END WHERE

     ! ------------------------------------------------------------------------------
     ! If any window is excluded, find the corresponding indices. This has to be done
     ! after the array assignements above because we need to know which indices to
     ! exclude from the final arrays, not the complete ones read from the HE4 file.
     ! ------------------------------------------------------------------------------
     omi_ccdpix_exclusion(ix,1:2) = -1
     IF ( MINVAL(ctrvar%fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winexc_lim(1),KIND=r8), 'GE', &
             omi_ccdpix_exclusion(ix,1) )
        CALL array_locate_r8 ( &
             nw, radref_wavl_ix, REAL(ctrvar%fit_winexc_lim(2),KIND=r8), 'LE', &
             omi_ccdpix_exclusion(ix,2) )
     END IF
     
  END DO

  locerrstat = omps_reader_status_reference  

END SUBROUTINE create_radiance_reference
