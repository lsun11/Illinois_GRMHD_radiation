subroutine output_data_to_file(filename,num_cols,col_names,output_data,header_flag)
  character                                :: filename*50
  integer                                  :: num_cols,header_flag
  real*8, dimension(num_cols)              :: output_data
  character, dimension(num_cols)           :: col_names*20
  integer i,strlen,u
  u = 1
!  write (*,100) output_data
!100 format (999E8.3E3)
!  i = len(filename)
!  do while (filename(i:i) .eq.  )
!     i = i - 1
!  enddo
!  strlen = i
!  do i=0,3
!     filename(1:strlen+i-1) = file_suffix(i+1)
!  end do
!  filename_final = filename // file_suffix //  
!  open (u, FILE=filename, STATUS=OLD, ACCESS=APPEND)
  open (u, FILE=filename, POSITION='APPEND')
  if(header_flag .eq. 1) then
     write(u,'(XXX999A)') col_names
  end if
  write(u,'(999ES20.10E3)') output_data
! USE FOR DEBUGGING, OUTPUTS 15 DIGITS:
!  write(u,(999E25.15E3)) output_data
  close(u)
end subroutine output_data_to_file
subroutine output_multiline_data_to_file(filename,num_cols,col_names,num_rows,output_data,header_flag)
  character                                :: filename*50
  integer                                  :: num_cols,num_rows,header_flag
  real*8, dimension(num_rows,num_cols)              :: output_data
  character, dimension(num_cols)           :: col_names*20
  integer i,strlen,u
  u = 1
  open (u, FILE=filename, POSITION='APPEND')
  if(header_flag .eq. 1) then
     write(u,'(XXX999A)') col_names
  end if
  do i=1,num_rows
     write(u,'(999ES20.10E3)') output_data(i,1:num_cols)
  end do
! USE FOR DEBUGGING, OUTPUTS 15 DIGITS:
!  write(u,(999E25.15E3)) output_data
  close(u)
end subroutine output_multiline_data_to_file
