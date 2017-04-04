!! MODULE to READ OMPS data
!! It may also contain the calculation
!! of the slit functions
!!
!! Created by gga June 2014
!!
!! All he5 internals are done by the belgium (BIRA)
!! module ReadH5dataset
!! with small modifications
!!
!! ------------------------------------------------
MODULE OMSAO_OMPS_SLIT_READER

  USE ReadH5dataset ! Contains the generic HDF5 reading routines

  IMPLICIT NONE

  CHARACTER(LEN=256) :: ompsSlitFN = &
       "/data/tempo2/ggonzale/OMPS_source/Kai/OMPS_L1B/tbl/OMPS_SLIT/BPS_TC_MACROPIXEL_V10_output_verif_final.hdf5"

  CHARACTER(LEN=256) :: ompsBcwlFN = &
       "/data/tempo2/ggonzale/OMPS_source/Kai/OMPS_L1B/tbl/OMPS_SLIT/CBC_TC_MACROPIXEL_V11_output_verif_final.hdf5"

  INTEGER, PUBLIC :: nbands, nXts, nSiltPts
  REAL(KIND=8), DIMENSION(1:196,1:36,1:41), PUBLIC :: slitMacroTab
  REAL(KIND=8), DIMENSION(1:196,1:36     ), PUBLIC :: wlMacroTab
  REAL(KIND=8), DIMENSION(           1:41), PUBLIC :: wlRelTab
  REAL(KIND=8),                             PUBLIC :: wlRelminTab, &
       wlRelmaxTab, &
       wlStepTab
  !    CHARACTER(LEN=10), PUBLIC :: dataOption = 'GROUND' !! 'FLIGHT'
  !! should use flight band pass
  CHARACTER(LEN=10), PUBLIC :: dataOption = 'FLIGHT'
  
  LOGICAL, PRIVATE :: SlitDataStaged = .FALSE. 
  PUBLIC :: OMPS_macroSlit_read
  PUBLIC :: OMPS_SlitF  !! single value
  PUBLIC :: OMPS_SlitA  !! array 
  PUBLIC :: OMPS_wl2band !! nearest band idx to wl
  PUBLIC :: OMPS_wlr2bands !! central band indices that cover
  !! input wavelength range.

CONTAINS

  SUBROUTINE OMPS_macroSlit_read( Slit_fn, wlBC_fn, & !input files
       nb, nx, npts,   & ! Number of bands, xtrack and points
       slitMacroData, wlMacro, wlRel,    & ! Slit function, wavelengths, wlrel
       wlRelmin, wlRelmax, wlStep,       & ! min, max, step
       LFAIL, msg) 

    CHARACTER(LEN=*),               INTENT(IN) :: Slit_fn, wlBC_fn
    INTEGER,                        INTENT(OUT) :: nb, nx, npts
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: slitMacroData
    REAL(KIND=8), DIMENSION(:,:  ), INTENT(OUT) :: wlMacro
    REAL(KIND=8), DIMENSION(:    ), INTENT(OUT) :: wlRel
    REAL(KIND=8),                   INTENT(OUT) :: wlRelmin, wlRelmax, &
         wlStep
    LOGICAL, INTENT(OUT) :: LFAIL 
    CHARACTER(LEN=*), INTENT(OUT) :: msg
    REAL(KIND=8), DIMENSION(:),   POINTER :: foo1d => NULL()
    REAL(KIND=8), DIMENSION(:,:), POINTER :: foo2d => NULL()
    
    !! local
!!$    INTEGER(HID_T) :: fileId
!!$    INTEGER(HID_T) :: dataGroupId
!!$    
!!$    INTEGER (HSIZE_T), DIMENSION(7) :: start, count
    
    real(KIND=8), dimension(1) :: foo
    integer,  dimension(1) :: ifoo
    integer :: NpixelMacro, NwaveRelMacro
    integer :: ic, ii,jj,k
    integer, DIMENSION(1:36) :: idx
    real(KIND=8), dimension(:,:), allocatable :: bandpassMacroFlight
    real(KIND=8), dimension(:  ), allocatable :: bandCenterMacro
    CHARACTER(LEN=256):: filename
    CHARACTER(LEN=20) :: modulename = 'OMPS_macroSlit_read'
    
    LFAIL = .FALSE.
    msg   = 'OK'
    !
    ! open HDF5-FORTRAN interface and file:
    !
  
    filename = Slit_fn
    ! Read NPIXELMACRO & NWAVERELMACRO
    CALL H5ReadDataset(filename, &
         "/DATA/NPIXELMACRO", ifoo(1))
    NpixelMacro = ifoo(1)
    CALL H5ReadDataset(filename, &
         "/DATA/NWAVERELMACRO", ifoo(1))
    NwaveRelMacro = ifoo(1)
    
    nx   = 36
    npts = NwaveRelMacro
    IF( MOD (NpixelMacro,38) /= 0 ) THEN
       LFAIL = .TRUE.
       msg = TRIM(modulename) //': expect NpixelMacro is a multiple of 38.'
       RETURN
    ENDIF
    nb   = NpixelMacro/38

    IF( SIZE( wlRel ) < npts ) THEN
       LFAIL = .TRUE.
       msg = TRIM(modulename) //': input szie wlRel too small.'
       RETURN
    ENDIF
    
    IF( SIZE( wlMacro, 1) < nb .OR. SIZE( wlMacro, 2) < nx ) THEN
       LFAIL = .TRUE.
       msg = TRIM(modulename) //': input szie wlMacro too small.'
       RETURN
    ENDIF
    
    IF( SIZE(slitMacroData,1)<nb .OR. SIZE(slitMacroData,2)<nx .OR. &
         SIZE(slitMacroData,3)<npts ) THEN
       LFAIL = .TRUE.
       msg = TRIM(modulename) //': input szie slitMacroData too small.'
       RETURN
    ENDIF
    
    ! Read WaveDeltaMacro, WaveMaxMacro and WaveMinMacro
    CALL H5ReadDataset(filename, &
         "/DATA/WAVEDELTAMACRO", foo(1))
    wlStep = foo(1)
    CALL H5ReadDataset(filename, &
         "/DATA/WAVEMAXMACRO", foo(1))
    wlRelMax = foo(1)
    CALL H5ReadDataset(filename, &
         "/DATA/WAVEMINMACRO", foo(1))
    wlRelMin = foo(1)

    IF( SIZE( wlRel ) < NwaveRelMacro ) THEN
       LFAIL = .TRUE.
       msg = TRIM(modulename) //': wlRel array size too small.'
       RETURN
    ENDIF

    ! ReadWaveRelMacro
    allocate( foo1d(1:SIZE(wlRel)))
    CALL H5ReadDataset(filename, &
         "/DATA/WAVERELMACRO", foo1d)
    wlRel = foo1d
    IF(associated( foo1d )) deallocate( foo1d )

    ! Read BandPassMacro_Flight
    allocate( bandpassMacroFlight(1:NpixelMacro,1:NwaveRelMacro))
    allocate( foo2d(1:NpixelMacro,1:NwaveRelMacro))
    CALL H5ReadDataset(filename, &
         "/DATA/BANDPASSMACRO_FLIGHT", &
         foo2d)
    bandpassMacroFlight = foo2d
    IF(associated( foo2d )) deallocate( foo2d )

    ! Copy data to slitMacroData from BandPassMacro_Flight
    DO ic = 1, nb
       ii = (ic - 1)*38+1; jj = ic*38
       idx(1 :18) = (/(k,k=ii,ii+17)/)
       idx(19:36) = (/(k,k=ii+20,jj)/)
       slitMacroData( ic,1:36,1:41) = &
            bandpassMacroFlight(idx(1:36),1:NwaveRelMacro)
    ENDDO
    ! Deallocate memory
    IF(allocated( bandpassMacroFlight)) deallocate( bandpassMacroFlight )

    ! Read band center Wavelength
    filename = wlBC_fn
    allocate( foo1d(1:NpixelMacro) )
    allocate( bandCenterMacro(1:NpixelMacro) )
    CALL H5ReadDataset(filename, &
         "/DATA/WAVEMACRO", foo1d)
    bandCenterMacro = foo1d
    ! Copy data to wlMacro from bandCenterMacro
    DO ic = 1, nb
       ii = (ic - 1)*38+1; jj = ic*38
       idx(1 :18) = (/(k,k=ii,ii+17)/)
       idx(19:36) = (/(k,k=ii+20,jj)/)
       wlMacro( ic,1:36) =  bandCenterMacro(idx(1:36))
    ENDDO
    !Deallocate memory
    IF(associated( foo1d )) deallocate( foo1d )
    IF(allocated( bandCenterMacro)) deallocate( bandCenterMacro )
       
    slitMacroTab = LOG( slitMacroTab )
    SlitDataStaged = .TRUE. 

  END SUBROUTINE OMPS_macroSlit_read

       !! The OMPS slit is determined by the band index (ib) and 
       !! the cross track position (iXt). An optional parameter
       !! slitAbscissaScaling expands or shrinks the spacing of
       !! tabulated wavelength spacing. Given an input
       !! wavelength relative to the central wavlength of the band,
       !! this subroutine returns the slit value at the input
       !! wavelength x. Log-linear interpolation is done on the
       !! tabulated when input x is on the scaled grid.
       FUNCTION OMPS_SlitF( x, ib, iXt, &
                            slitAbscissaScaling ) RESULT( slitF )
         REAL(KIND=8), INTENT(IN) :: x 
         !!scaling of x of slit table to shrink or expand the slit width
         REAL(KIND=8), INTENT(IN), OPTIONAL :: slitAbscissaScaling 
         INTEGER, INTENT(IN) :: ib, iXt

         !! result
         REAL(KIND=8) :: slitF
         REAL(KIND=8) :: xScaled, frc
         INTEGER :: il, ih
     
         CHARACTER(LEN=10) :: modulename = 'OMPS_SlitF'
         LOGICAL :: LFAIL 
         CHARACTER(LEN=255) :: msg

         IF( .NOT.  SlitDataStaged ) THEN
            CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
                          nbands, nXts, nSiltPts, slitMacroTab,&
                          wlMacroTab, wlRelTab, wlRelminTab,   &
                          wlRelmaxTab, wlStepTab, LFAIL, msg   )
            IF( LFAIL ) THEN
               WRITE(*,*) TRIM(modulename) // ': read OMPS slit file failed.'
               SlitDataStaged = .FALSE. 
               slitF = -1.d0
               RETURN
            ENDIF
         ENDIF
         IF( ib < 1 .OR. ib >  nbands ) THEN
            WRITE(*,*) TRIM(modulename) // ': input band out of range:', ib
            slitF = -1.d0
            RETURN
         ENDIF
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            WRITE(*,*) TRIM(modulename) // ': input iXt out of range:', iXt
            slitF = -1.d0
            RETURN
         ENDIF

         IF( PRESENT( slitAbscissaScaling ) ) THEN
            xScaled = x/slitAbscissaScaling
            IF( xScaled < wlRelminTab .OR. xScaled > wlRelmaxTab ) THEN
               slitF = 0.d0
            ELSE
               frc = (xScaled - wlRelminTab )/ wlStepTab 
               il  = INT(frc)
               frc = frc - DBLE(il)
               il = il + 1; ih = il + 1
               slitF = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                     + slitMacroTab(ib,iXt,ih)*      frc
               slitF = EXP( slitF )
            ENDIF
         ELSE
            IF( x < wlRelminTab .OR. x > wlRelmaxTab ) THEN
               slitF = 0
            ELSE
               frc = (x - wlRelminTab )/ wlStepTab 
               il  = INT(frc)
               frc = frc - DBLE(il)
               il = il + 1; ih = il + 1
               slitF = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                     + slitMacroTab(ib,iXt,ih)*      frc
               slitF = EXP( slitF )
            ENDIF
         ENDIF
         RETURN 
       END FUNCTION OMPS_SlitF

       !! The OMPS slit is determined by the band index (ib) and 
       !! the cross track position (iXt). An optional parameter
       !! slitAbscissaScaling expands or shrinks the spacing of
       !! tabulated wavelength spacing. Given an array of input
       !! wavelength relative to the central wavlength of the band,
       !! this subroutine returns the slit values at the input
       !! wavelengths xa. Log-linear interpolation is done on the
       !! tabulated when input xa is on the scaled grid.
       FUNCTION OMPS_SlitA( xa, ib, iXt, &
                            slitAbscissaScaling ) RESULT( slitA )
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xa 
          !!scaling of x of slit table to shrink or expand the slit width
         REAL(KIND=8), INTENT(IN), OPTIONAL :: slitAbscissaScaling
         INTEGER, INTENT(IN) :: ib, iXt

         !! result
         REAL(KIND=8), DIMENSION(SIZE(xa)) :: slitA
         REAL(KIND=8)                      :: xScaled, frc, slitFoo
         INTEGER :: ixa, nxa, il, ih
     
         CHARACTER(LEN=10) :: modulename = 'OMPS_SlitA'
         LOGICAL :: LFAIL 
         CHARACTER(LEN=255) :: msg

         IF( .NOT.  SlitDataStaged ) THEN
            CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
                          nbands, nXts, nSiltPts, slitMacroTab,&
                          wlMacroTab, wlRelTab, wlRelminTab,   &
                          wlRelmaxTab, wlStepTab, LFAIL, msg)
            IF( LFAIL ) THEN
               WRITE(*,*) TRIM(modulename) // ': read OMPS slit file failed.'
               SlitDataStaged = .FALSE. 
               slitA(:) = -1.d0
               RETURN
            ENDIF
         ENDIF

         IF( ib < 1 .OR. ib >  nbands ) THEN
            WRITE(*,*) TRIM(modulename) // ': input band out of range:', ib
            slitA(:) = -1.d0
            RETURN
         ENDIF
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            WRITE(*,*) TRIM(modulename) // ': input iXt out of range:', iXt
            slitA(:) = -1.d0
            RETURN
         ENDIF

         nxa = SIZE(xa(:))
         IF( PRESENT( slitAbscissaScaling ) ) THEN
            DO ixa = 1, nxa
              xScaled = xa(ixa)/slitAbscissaScaling
              IF( xScaled < wlRelminTab .OR. xScaled > wlRelmaxTab ) THEN
                 slitA(ixa) = 0.d0
              ELSE
                 frc = (xScaled - wlRelminTab )/ wlStepTab 
                 il  = INT(frc)
                 frc = frc - DBLE(il)
                 il = il + 1; ih = il + 1
                 slitFoo = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                         + slitMacroTab(ib,iXt,ih)*      frc
                 slitA(ixa) = EXP( slitFoo )
              ENDIF
            ENDDO
         ELSE
            DO ixa = 1, nxa
              IF( xa(ixa) < wlRelminTab .OR. xa(ixa) > wlRelmaxTab ) THEN
                 slitA(ixa) = 0.d0
              ELSE
                 frc = (xa(ixa) - wlRelminTab )/ wlStepTab 
                 il  = INT(frc)
                 frc = frc - DBLE(il)
                 il = il + 1; ih = il + 1
                 slitFoo = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                         + slitMacroTab(ib,iXt,ih)*      frc
                 slitA(ixa) = EXP( slitFoo )
              ENDIF
            ENDDO
         ENDIF
         RETURN 
       END FUNCTION OMPS_SlitA

       !! Given a central wavelength wlc, and the cross track
       !! position iXt, this function returns the closest band index
       !! (therefore the tabulated slit data)
       FUNCTION OMPS_wl2band( wlc, iXt, LFAIL, msg ) RESULT(ib)
         REAL(KIND=8), INTENT(IN) :: wlc
         INTEGER, INTENT(IN) :: iXt
         LOGICAL, INTENT(OUT) :: LFAIL
         CHARACTER(LEN=*), INTENT(OUT) :: msg

         !! local
         INTEGER :: ib
         CHARACTER(LEN=12) :: modulename = 'OMPS_wl2band'
 
         LFAIL = .FALSE.; msg = 'OK'
         IF( .NOT.  SlitDataStaged ) THEN
            CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
                          nbands, nXts, nSiltPts, slitMacroTab,&
                          wlMacroTab, wlRelTab, wlRelminTab,   &
                          wlRelmaxTab, wlStepTab, LFAIL, msg)
            IF( LFAIL ) THEN
               msg = TRIM(modulename) //':'// TRIM(msg)
               SlitDataStaged = .FALSE. 
               ib = -1
               RETURN
            ENDIF
         ENDIF
         
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            LFAIL = .TRUE.
            WRITE(msg,*) TRIM(modulename) // ': iXt out of range:', iXt
            ib = -1
            RETURN
         ENDIF 
           
         IF( wlc <= wlMacroTab(1,iXt ) ) THEN
            ib = 1
         ELSE IF( wlc >= wlMacroTab(nbands,iXt ) ) THEN
            ib = nbands
         ELSE
            ib = MINVAL(MAXLOC(wlMacroTab(1:nbands,iXt ), &
                        MASK= (wlMacroTab(1:nbands,iXt ) <= wlc) ))

            IF( ABS(wlc-wlMacroTab(ib,iXt))  > &
                ABS(wlMacroTab(ib+1,iXt)-wlc) ) THEN 
               ib = ib+1
            ENDIF
         ENDIF
         RETURN 
       END FUNCTION OMPS_wl2band

       !! Given a wavlength range wlr, and the cross track position iXt
       !! this routine returns the number of bands (nwls) and the band center
       !! wavelengths (wls) as well as the band index (bandIdx) that covers
       !! (i.e., just bigger than ) this input range wlr.
       SUBROUTINE OMPS_wlr2bands( wlr, iXt, wls, bandIdx, nwls, LFAIL, msg ) 
         REAL(KIND=8), DIMENSION(2), INTENT(IN) :: wlr
         INTEGER, INTENT(IN) :: iXt
         INTEGER, INTENT(OUT) :: nwls
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: wls
         INTEGER,      DIMENSION(:), INTENT(OUT) :: bandIdx
         LOGICAL, INTENT(OUT) :: LFAIL
         CHARACTER(LEN=*), INTENT(OUT) :: msg
  
         !! local
         INTEGER :: ibs, ibe, ii
         CHARACTER(LEN=14) :: modulename = 'OMPS_wlr2bands'
 
         LFAIL = .FALSE.; msg = 'OK'
         IF( .NOT.  SlitDataStaged ) THEN
            CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
                          nbands, nXts, nSiltPts, slitMacroTab,&
                          wlMacroTab, wlRelTab, wlRelminTab,   &
                          wlRelmaxTab, wlStepTab, LFAIL, msg)
            IF( LFAIL ) THEN
               msg = (modulename) // ':'// TRIM(msg)
               SlitDataStaged = .FALSE. 
               RETURN
            ENDIF
         ENDIF
         
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            LFAIL = .TRUE.
            WRITE( msg, *) TRIM(modulename) // ': iXt out of range:', iXt
            nwls = -1
            RETURN
         ENDIF 
           
         IF( wlr(1) <= wlMacroTab(1,iXt ) ) THEN
            ibs = 1
         ELSE IF( wlr(1) >= wlMacroTab(nbands,iXt ) ) THEN
            ibs = nbands
         ELSE
            ibs = MAXVAL(MAXLOC(wlMacroTab(1:nbands,iXt ), &
                         MASK= (wlMacroTab(1:nbands,iXt ) <= wlr(1)) ))
         ENDIF

         IF( wlr(2) <= wlMacroTab(1,iXt ) ) THEN
            ibe = 1
         ELSE IF( wlr(2) >= wlMacroTab(nbands,iXt ) ) THEN
            ibe = nbands
         ELSE
            ibe = MINVAL(MINLOC(wlMacroTab(1:nbands,iXt ), &
                         MASK= (wlMacroTab(1:nbands,iXt ) >= wlr(2)) ))
         ENDIF

         nwls = ibe - ibs + 1

         IF( SIZE( wls ) < nwls .oR. SIZE( bandIdx) < nwls )  THEN
            LFAIL = .FALSE.
            msg = TRIM(modulename) //': input size wls or bandIdx too small' 
            wls(:) = -1.d0
            bandIdx(:) = -1
            RETURN
         ENDIF

         wls(1:nwls) = wlMacroTab(ibs:ibe,iXt )
         bandIdx(1:nwls) = (/ (ii, ii = ibs, ibe)/)
         RETURN 
       END SUBROUTINE OMPS_wlr2bands
  
END MODULE OMSAO_OMPS_SLIT_READER
