MODULE OMPS_macroSlitNOAA_m
     USE HDF5
     USE H5Util_class
     IMPLICIT NONE
     CHARACTER(LEN=80) :: ompsSlitFN = &
       "../tbl/BPS_TC_MACROPIXEL_V10_output_verif_final.hdf5"

     CHARACTER(LEN=80) :: ompsBcwlFN = &
       "../tbl/CBC_TC_MACROPIXEL_V11_output_verif_final.hdf5"

     type(H5SDS_T), PRIVATE :: sds_bandpassMacroFlight = &
            H5SDS_T("BANDPASSMACRO_FLIGHT",     &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_bandpassMacroGround = &
            H5SDS_T("BANDPASSMACRO_GROUND",     &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_NpixelMacro =         &
            H5SDS_T("NPIXELMACRO",              &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_NwaveRelMacro =       &
            H5SDS_T("NWAVERELMACRO",            &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_WaveDeltaMacro =      &
            H5SDS_T("WAVEDELTAMACRO",           &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_WaveMaxMacro =        &
            H5SDS_T("WAVEMAXMACRO",             &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_WaveMinMacro =        &
            H5SDS_T("WAVEMINMACRO",             &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)
     type(H5SDS_T), PRIVATE :: sds_WaveRelMacro =        &
            H5SDS_T("WAVERELMACRO",             &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)

     type(H5SDS_T), PRIVATE :: sds_WaveMacro =        &
            H5SDS_T("WAVEMACRO",             &
                    -1,-1,-1,(/0,0,0,0,0,0,0/),"N/A","N/A",(/0,0/),0)

     
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

       SUBROUTINE OMPS_macroSlit_read( Slit_fn, wlBC_fn, nb, nx, npts,   &
                                     slitMacroData, wlMacro, wlRel,    &
                                     wlRelmin, wlRelmax, wlStep,       &
                                     LFAIL, msg, option ) 
!        USE DataTypeDef
!      FUNCTION OMPS_macroSlit_read( Slit_fn, wlBC_fn ) RESULT( status )
         CHARACTER(LEN=*), INTENT(IN) :: Slit_fn, wlBC_fn
         INTEGER,                    INTENT(OUT) :: nb, nx, npts
         REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: slitMacroData
         REAL(KIND=8), DIMENSION(:,:  ), INTENT(OUT) :: wlMacro
         REAL(KIND=8), DIMENSION(:    ), INTENT(OUT) :: wlRel
         REAL(KIND=8),                   INTENT(OUT) :: wlRelmin, wlRelmax, &
                                         wlStep
         LOGICAL, INTENT(OUT) :: LFAIL 
         CHARACTER(LEN=*), INTENT(OUT) :: msg
         CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: option
     
         !! local
         INTEGER(HID_T) :: fileId
         INTEGER(HID_T) :: dataGroupId

         INTEGER :: hdferr, errorStatus

         real(KIND=8), dimension(1) :: foo
         integer,  dimension(1) :: ifoo
         integer :: NpixelMacro, NwaveRelMacro
         integer :: ic, ii,jj,k
         integer, DIMENSION(1:36) :: idx
         real(KIND=8), dimension(:,:), allocatable :: bandpassMacroFlight
         real(KIND=8), dimension(:,:), allocatable :: bandpassMacroGround
         real(KIND=8), dimension(:  ), allocatable :: bandCenterMacro
         CHARACTER(LEN=20) :: modulename = 'OMPS_macroSlit_read'
        
         LFAIL = .FALSE.
         msg   = 'OK'
         !
         ! open HDF5-FORTRAN interface and file:
         !

         call h5open_f(hdferr)
         call h5fopen_f(Slit_fn, H5F_ACC_RDONLY_F, fileId, hdferr)

         IF( hdferr /= 0 ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': fail open failed : '//TRIM(Slit_fn) 
            RETURN
         ENDIF
         !
         ! acquire group ID: there are four data groups in a TC SDR file: 
         !

         call h5gopen_f(fileId, "DATA", dataGroupId, hdferr)
         IF( hdferr /= 0 ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': open DATA group failed.'
            RETURN
         ENDIF


         errorStatus   = H5Util_selectDataset(dataGroupId, sds_NpixelMacro)
         errorStatus   = H5Util_readDataset(sds_NpixelMacro, ifoo )
         NpixelMacro   = ifoo(1)
         errorstatus   = h5Util_disposeDataset(sds_NpixelMacro)

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_NwaveRelMacro)
         errorStatus   = H5Util_readDataset(sds_NwaveRelMacro, ifoo )
         NwaveRelMacro = ifoo(1)
         errorstatus   = h5Util_disposeDataset(sds_NwaveRelMacro)

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

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_WaveDeltaMacro)
         errorStatus   = H5Util_readDataset(sds_WaveDeltaMacro, foo )
         wlStep        = foo(1)
         errorstatus   = h5Util_disposeDataset(sds_WaveDeltaMacro)

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_WaveMaxMacro)
         errorStatus   = H5Util_readDataset(sds_WaveMaxMacro, foo )
         wlRelMax      = foo(1)
         errorstatus   = h5Util_disposeDataset(sds_WaveMaxMacro)

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_WaveMinMacro)
         errorStatus   = H5Util_readDataset(sds_WaveMinMacro, foo )
         wlRelMin      = foo(1)
         errorstatus   = h5Util_disposeDataset(sds_WaveMinMacro)

         IF( SIZE( wlRel ) < NwaveRelMacro ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': wlRel array size too small.'
            RETURN
         ENDIF

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_WaveRelMacro)
         errorStatus   = H5Util_readDataset(sds_WaveRelMacro, wlRel )
         errorstatus   = h5Util_disposeDataset(sds_WaveRelMacro)
   
         IF( PRESENT( option ) ) THEN
            IF( option == 'FLIGHT' .OR. option == 'flight' ) THEN
               errorStatus   = H5Util_selectDataset( dataGroupId, &
                                              sds_bandpassMacroFlight )
        
               allocate( bandpassMacroFlight(1:NpixelMacro,1:NwaveRelMacro))
               errorStatus   = H5Util_readDataset(sds_bandpassMacroFlight, &
                                                  bandpassMacroFlight )
               errorstatus   = h5Util_disposeDataset(sds_bandpassMacroFlight)
               DO ic = 1, nb
                 ii = (ic - 1)*38+1; jj = ic*38
                 idx(1 :18) = (/(k,k=ii,ii+17)/)
                 idx(19:36) = (/(k,k=ii+20,jj)/)
                 slitMacroData( ic,1:36,1:41) =  &
                       bandpassMacroFlight(idx(1:36),1:NwaveRelMacro)
               ENDDO

            ELSE IF( option == 'GROUND' .OR. option == 'ground' ) THEN
               errorStatus   = H5Util_selectDataset( dataGroupId, &
                                              sds_bandpassMacroGround )

               allocate( bandpassMacroGround(1:NpixelMacro,1:NwaveRelMacro))
               errorStatus   = H5Util_readDataset(sds_bandpassMacroGround, &
                                                  bandpassMacroGround )
               errorstatus   = h5Util_disposeDataset(sds_bandpassMacroGround)

               !! slitMacroData( 1:196,1:36,1:41) 
               DO ic = 1, nb
                 ii = (ic - 1)*38+1; jj = ic*38
                 idx(1 :18) = (/(k,k=ii,ii+17)/)
                 idx(19:36) = (/(k,k=ii+20,jj)/)
                 slitMacroData( ic,1:36,1:41) =  &
                       bandpassMacroGround(idx(1:36),1:NwaveRelMacro)
               ENDDO
            ELSE
               LFAIL = .TRUE.
               msg = TRIM(modulename) // ': unknown option:'//TRIM(option)
               RETURN
            ENDIF
         ELSE
            errorStatus   = H5Util_selectDataset( dataGroupId, &
                                              sds_bandpassMacroFlight )

            allocate( bandpassMacroFlight(1:NpixelMacro,1:NwaveRelMacro))
            errorStatus   = H5Util_readDataset(sds_bandpassMacroFlight, &
                                                  bandpassMacroFlight )
            errorstatus   = h5Util_disposeDataset(sds_bandpassMacroFlight)
            DO ic = 1, nb
              ii = (ic - 1)*38+1; jj = ic*38
              idx(1 :18) = (/(k,k=ii,ii+17)/)
              idx(19:36) = (/(k,k=ii+20,jj)/)
              slitMacroData( ic,1:36,1:41) = &
                   bandpassMacroFlight(idx(1:36),1:NwaveRelMacro)
            ENDDO

         ENDIF

         ! shutdown: 

         ! free up memory: 
         IF(allocated( bandpassMacroFlight)) deallocate( bandpassMacroFlight )
         IF(allocated( bandpassMacroGround)) deallocate( bandpassMacroGround )

         ! dispose group ID:

         call h5gclose_f(dataGroupId,     hdferr)

         ! close file and HDF5-FORTRAN interface: 

         call h5fclose_f(fileId, hdferr)
         call h5close_f(hdferr)

         !! Read Band Center Wavelength
         call h5open_f(hdferr)
         call h5fopen_f(wlBC_fn, H5F_ACC_RDONLY_F, fileId, hdferr)

         IF( hdferr /= 0 ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': fail open failed : '//TRIM(wlBC_fn) 
            RETURN
         ENDIF
         call h5gopen_f(fileId, "DATA", dataGroupId, hdferr)
         IF( hdferr /= 0 ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': open DATA group failed.'
            RETURN
         ENDIF

         errorStatus   = H5Util_selectDataset(dataGroupId, sds_WaveMacro)
         IF( sds_WaveMacro%dims(1) /= NpixelMacro ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) //': sds_WaveMacro%dim(1) /= NpixelMacro'
            RETURN
         ENDIF
         allocate( bandCenterMacro(1:NpixelMacro) )
         errorStatus   = H5Util_readDataset(sds_WaveMacro, bandCenterMacro )
         errorstatus   = h5Util_disposeDataset(sds_WaveMacro)

         DO ic = 1, nb
           ii = (ic - 1)*38+1; jj = ic*38
           idx(1 :18) = (/(k,k=ii,ii+17)/)
           idx(19:36) = (/(k,k=ii+20,jj)/)
           wlMacro( ic,1:36) =  bandCenterMacro(idx(1:36))
         ENDDO

         IF(allocated( bandCenterMacro)) deallocate( bandCenterMacro )

         call h5gclose_f(dataGroupId,     hdferr)
         call h5fclose_f(fileId, hdferr)
         call h5close_f(hdferr)

!        slitMacroTab = LOG( slitMacroTab )
         SlitDataStaged = .TRUE. 
       END SUBROUTINE OMPS_macroSlit_read

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
                          wlRelmaxTab, wlStepTab, LFAIL, msg,  &
                          dataOption )
            IF( LFAIL ) THEN
               WRITE(*,*) TRIM(modulename) // ': read OMPS slit file failed.'
               SlitDataStaged = .FALSE. 
               slitF = -1.d0
               RETURN
            ENDIF
         ENDIF
         IF( ib < 1 .OR. ib >  nbands ) THEN
            WRITE(*,*) TRIM(modulename) // ': input band out of range.'
            slitF = -1.d0
            RETURN
         ENDIF
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            WRITE(*,*) TRIM(modulename) // ': input iXt out of range.'
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
               IF( il < 1 ) THEN
                  slitF = slitMacroTab(ib,iXt,1)
               ELSE IF( ih > 41 ) THEN
                  slitF = slitMacroTab(ib,iXt,41)
               ELSE  
                   slitF = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                         + slitMacroTab(ib,iXt,ih)*      frc
               ENDIF
            ENDIF
!            print * , x, wlsteptab
         ELSE
            IF( x < wlRelminTab .OR. x > wlRelmaxTab ) THEN
               slitF = 0
            ELSE
               IF( x <= wlRelminTab ) THEN
                  slitF = slitMacroTab(ib,iXt,1)
               ELSE IF( x >= wlRelmaxTab ) THEN
                  slitF = slitMacroTab(ib,iXt,41)
               ELSE 
                  frc = (x - wlRelminTab )/ wlStepTab 
                  il  = INT(frc)
                  frc = frc - DBLE(il)
                  il = il + 1; ih = il + 1
                  IF( il < 1 ) THEN
                     slitF = slitMacroTab(ib,iXt,1)
                  ELSE IF( ih > 41 ) THEN
                     slitF = slitMacroTab(ib,iXt,41)
                  ELSE  
                     slitF = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                           + slitMacroTab(ib,iXt,ih)*      frc
                  ENDIF
   
              !!  WRITE(*,*) 'il, ih =', il, ih, frc
              !!  slitF = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
              !!        + slitMacroTab(ib,iXt,ih)*      frc
               ENDIF
            ENDIF
         
         ENDIF
          
         RETURN 

       END FUNCTION OMPS_SlitF

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
                          wlRelmaxTab, wlStepTab, LFAIL, msg,  &
                          dataOption )
            IF( LFAIL ) THEN
               WRITE(*,*) TRIM(modulename) // ': read OMPS slit file failed.'
               SlitDataStaged = .FALSE. 
               slitA(:) = -1.d0
               RETURN
            ENDIF
         ENDIF

         IF( ib < 1 .OR. ib >  nbands ) THEN
            WRITE(*,*) TRIM(modulename) // ': input band out of range.'
            slitA(:) = -1.d0
            RETURN
         ENDIF
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            WRITE(*,*) TRIM(modulename) // ': input iXt out of range.'
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
                 IF( il < 1 ) THEN
                    slitFoo = slitMacroTab(ib,iXt,1)
                 ELSE IF( ih > 41 ) THEN
                    slitFoo = slitMacroTab(ib,iXt,41)
                 ELSE  
                    slitFoo = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                            + slitMacroTab(ib,iXt,ih)*      frc
                 ENDIF
                 slitA(ixa) = slitFoo 
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
                 IF( il < 1 ) THEN
                    slitFoo = slitMacroTab(ib,iXt,1)
                 ELSE IF( ih > 41 ) THEN
                    slitFoo = slitMacroTab(ib,iXt,41)
                 ELSE
                    slitFoo = slitMacroTab(ib,iXt,il)*(1.d0-frc) &
                            + slitMacroTab(ib,iXt,ih)*      frc
                 ENDIF
                 slitA(ixa) = slitFoo 
              ENDIF
            ENDDO
         ENDIF
         RETURN 
       END FUNCTION OMPS_SlitA

       FUNCTION OMPS_wl2band( wlc, iXt, LFAIL, msg ) RESULT(ib)
         REAL(KIND=8), INTENT(IN) :: wlc
         INTEGER, INTENT(IN) :: iXt
         LOGICAL, INTENT(OUT) :: LFAIL
         CHARACTER(LEN=*), INTENT(OUT) :: msg

         !! local
         INTEGER :: ib
         CHARACTER(LEN=12) :: modulename = 'OMPS_wl2band'

         ib = 0
 
         LFAIL = .FALSE.; msg = 'OK'
         IF( .NOT.  SlitDataStaged ) THEN
            CALL OMPS_macroSlit_read( ompsSlitFN, ompsBcwlFN,  &
                          nbands, nXts, nSiltPts, slitMacroTab,&
                          wlMacroTab, wlRelTab, wlRelminTab,   &
                          wlRelmaxTab, wlStepTab, LFAIL, msg,  &
                          dataOption )
            IF( LFAIL ) THEN
               msg = TRIM(modulename) //':'// TRIM(msg)
               SlitDataStaged = .FALSE. 
               RETURN
            ENDIF
         ENDIF
         
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) // ': iXt out of range.'
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
                          wlRelmaxTab, wlStepTab, LFAIL, msg,  &
                          dataOption )
            IF( LFAIL ) THEN
               msg = (modulename) // ':'// TRIM(msg)
               SlitDataStaged = .FALSE. 
               RETURN
            ENDIF
         ENDIF
         
         IF( iXt < 1 .OR. iXt > nXts ) THEN
            LFAIL = .TRUE.
            msg = TRIM(modulename) // ': iXt out of range.'
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
    
END MODULE OMPS_macroSlitNOAA_m
