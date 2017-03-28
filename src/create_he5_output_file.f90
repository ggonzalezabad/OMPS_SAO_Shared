PROGRAM create_he5_output_file !(file_name, nt, nx, nl)

!!$  USE OMSAO_precision_module
!!$  USE OMSAO_errstat_module
  USE HDF5

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER(LEN=256), PARAMETER :: file_name='test.he5'
  INTEGER  (KIND=4) , PARAMETER :: nt=400, nx=36, nl=45

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=4) :: he5error

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=4), PARAMETER :: maxwvl = 500
  INTEGER (KIND=4)            :: errstat, rank, fillval_int
  INTEGER (HID_T) :: file_id, group_id, dspace_id, dset_id, list_id, type_id
  REAL    (KIND=8) :: fillval

  ! ------
  ! Groups
  ! ------
  CHARACTER(LEN=12), PARAMETER :: geolocation = '/Geolocation'
  CHARACTER(LEN=5),  PARAMETER :: datafields  = '/Data'

  ! ----------------
  ! Data field names
  ! ----------------
  ! Geolocation fields
  ! ------------------
  CHARACTER(LEN=8),  PARAMETER :: lat  = 'Latitude'
  CHARACTER(LEN=9),  PARAMETER :: lon  = 'Longitude'
  CHARACTER(LEN=17), PARAMETER :: saa  = 'SolarAzimuthAngle'
  CHARACTER(LEN=16), PARAMETER :: sza  = 'SolarZenithAngle'
  CHARACTER(LEN=18), PARAMETER :: sca  = 'SpacecraftAltitude'
  CHARACTER(LEN=15), PARAMETER :: trp  = 'TerrainPressure'
  CHARACTER(LEN=19), PARAMETER :: vaa  = 'ViewingAzimuthAngle'
  CHARACTER(LEN=18), PARAMETER :: vza  = 'ViewingZenithAngle'
  CHARACTER(LEN=23), PARAMETER :: gpqf = 'GroundPixelQualityFlags'
  CHARACTER(LEN=22), PARAMETER :: iqf  = 'InstrumentQualityFlags'
  CHARACTER(LEN=16), PARAMETER :: sed  = 'SunEarthDistance'

  ! -----------------
  ! Define dimensions
  ! -----------------
  INTEGER(HSIZE_T), DIMENSION(1) :: he5_1D
  INTEGER(HSIZE_T), DIMENSION(2) :: he5_2D
  INTEGER(HSIZE_T), DIMENSION(3) :: he5_3D

  he5error = pge_errstat_ok; errstat = pge_errstat_ok
!!$  file_name = 'test_'//file_name

  ! Initialize HDF5 interface
  CALL h5open_f(errstat)
  ! Create output file
  CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, errstat)
  ! Create geolocation group
  CALL h5gcreate_f(file_id, geolocation, group_id, errstat)

  ! Create datasets in geolocation group
  rank = 1
  he5_1d(1) = 1
  fillval = -9999
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, list_id, errstat)
  CALL h5pset_fill_value_f(list_id, H5T_NATIVE_DOUBLE, fillval, errstat)
  CALL h5screate_simple_f(rank, he5_1d, dspace_id, errstat)
  CALL h5dcreate_f(group_id, sed, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)
  CALL h5sclose_f(dspace_id, errstat)

  rank = 1
  he5_1d(1) = nt
  CALL h5screate_simple_f(rank, he5_1d, dspace_id, errstat)
  CALL h5dcreate_f(group_id, sca, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)
  CALL h5sclose_f(dspace_id, errstat)

  CALL h5screate_simple_f(rank, he5_1d, dspace_id, errstat)
  CALL h5dcreate_f(group_id, iqf, H5T_NATIVE_INTEGER, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)
  CALL h5sclose_f(dspace_id, errstat)

  rank = 2
  he5_2d(1) = nx; he5_2d(2) = nt 
  CALL h5screate_simple_f(rank, he5_2d, dspace_id, errstat)
  CALL h5dcreate_f(group_id, lat, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, lon, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id )
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, saa, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, sza, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, vaa, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, vza, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)

  CALL h5dcreate_f(group_id, trp, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)
  CALL h5sclose_f(dspace_id, errstat)
 
  CALL h5screate_simple_f(rank, he5_2d, dspace_id, errstat)
  CALL h5dcreate_f(group_id, gpqf, H5T_NATIVE_INTEGER, dspace_id, dset_id, errstat, list_id)
  CALL h5dclose_f(dset_id,errstat)
  CALL h5sclose_f(dspace_id, errstat)
  CALL h5gclose_f(group_id, errstat)


  ! Create datasets in Data group
  CALL h5gcreate_f(file_id, datafields, group_id, errstat)
  write(*,*) errstat

  CALL h5gclose_f(group_id, errstat)
  write(*,*) errstat
  
  ! Close the file
  CALL h5fclose_f(file_id, errstat)
  write(*,*) errstat

END PROGRAM create_he5_output_file
