Module lk_hdfOpenField_mod
  use hdf5
  implicit none

  interface lk_hdfOpenField
     module procedure lk_openScalar
     module procedure lk_openVector
     module procedure lk_open3Darray
  end interface lk_hdfOpenField

  private
  public :: lk_hdfOpenField

contains

  integer function lk_openScalar(fileName,fieldName,myScalar)
    character*200, intent(in) :: fileName, fieldName
    real(kind=8), intent(out) :: myScalar

    character*200 :: cleanFileName, cleanFieldName
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: error

    cleanFileName = trim(adjustl(fileName))
    cleanFieldName = trim(adjustl(fieldName))

    !INITIALIZE THE HDF LIBRARY
    call h5open_f(error)
    call h5fopen_f(cleanFileName, H5F_ACC_RDONLY_F, file_id, error)
    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Could not read the file with the EOS table"
       lk_openScalar=-1
       return
    endif
    call h5dopen_f(file_id, cleanFieldName, dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)

    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Problem opening the scalar dataset(#Err1)", cleanFieldName
       lk_openScalar=-1
       return
    endif

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, myScalar, dims, error)

    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Problem opening the vector dataset(#Err2)", cleanFieldName
       lk_openScalar=-1
       return
    endif
    !Close dataset and file
    call h5dclose_f(dset_id,error)
    call h5fclose_f(file_id,error)

    lk_openScalar = error

  end function lk_openScalar


  integer function lk_openVector(fileName,fieldName,myVector)
    character*200, intent(in) :: fileName, fieldName
    real(kind=8), allocatable, dimension(:), intent(out) :: myVector

    character*200 :: cleanFileName, cleanFieldName
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: dims1(1)
    integer(HSIZE_T) :: nElements
    integer :: error

    !Initialize the number of elements
    nElements = 0

    cleanFileName = trim(adjustl(fileName))
    cleanFieldName = trim(adjustl(fieldName))

    !INITIALIZE THE HDF LIBRARY
    call h5open_f(error)
    call h5fopen_f(cleanFileName, H5F_ACC_RDONLY_F, file_id, error)
    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Could not read the file with the EOS table"
       lk_openVector=-1
       return
    endif
    call h5dopen_f(file_id, cleanFieldName, dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_npoints_f(dspace_id, nElements, error)
    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Problem opening the vector dataset(#Err1)", cleanFieldName
       lk_openVector=-1
       return
    endif
    allocate(myVector(nElements))
    dims1(1) = nElements
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, myVector, dims1, error)

    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Problem opening the vector dataset(#Err2)", cleanFieldName
       lk_openVector=-1
       return
    endif
    !Close dataset and file
    call h5dclose_f(dset_id,error)
    call h5fclose_f(file_id,error)

    lk_openVector = error

  end function lk_openVector

  integer function lk_open3Darray(fileName,fieldName,myArray)
    character*200, intent(in) :: fileName, fieldName
    real(kind=8), dimension(:,:,:), intent(out) :: myArray

    character*200 :: cleanFileName, cleanFieldName
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T):: dims(3), maxdims(3)
    integer :: error

    cleanFileName = trim(adjustl(fileName))
    cleanFieldName = trim(adjustl(fieldName))

    !INITIALIZE THE HDF LIBRARY
    call h5open_f(error)
    call h5fopen_f(cleanFileName, H5F_ACC_RDONLY_F, file_id, error)
    if(error.eq.-1) then
       write(*,*) "hdfOpenField::Could not read the file with the EOS table"
       lk_open3Darray=-1
       return
    endif
    call h5dopen_f(file_id, cleanFieldName, dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
    if(error.eq.-1) then
        write(*,*) "hdfOpenField::Problem opening the array dataset(#Err1)", cleanFieldName
       lk_open3Darray=-1
       return
    endif
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, myArray, dims, error)

   if(error.eq.-1) then
        write(*,*) "hdfOpenField::Problem opening the array dataset(#Err2)", cleanFieldName
       lk_open3Darray=-1
       return
    endif
    !Close dataset and file
    call h5dclose_f(dset_id,error)
    call h5fclose_f(file_id,error)

    lk_open3Darray = error
  end function lk_open3Darray

end Module lk_hdfOpenField_mod
