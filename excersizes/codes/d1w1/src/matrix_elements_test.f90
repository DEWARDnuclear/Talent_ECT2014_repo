PROGRAM matrix_elements_test

  USE types
  USE symmat_mod

  IMPLICIT NONE


  ! Calling calculate_eigs_.. with an allocatable matrix causes segfault when the matrix is to be printed again. Using a normal matrix works. Error was with using incorrect work array in second call to dsyev, this must have messed with the stack.



  REAL(kind=r_kind), allocatable :: A(:,:), A_eigs(:)
  INTEGER :: A_size
  INTEGER :: II
  !REAL(kind=r_kind) :: A(5,5), A_eigs(5)
  !INTEGER :: A_size
  A_size = 5

  ALLOCATE(A(A_size,A_size),A_eigs(A_size))

  CALL generate_symmetrix_matrix(A_size,A)
  WRITE(*,*) 'The matrix'
  CALL print_matrix(A_size,A)

  CALL calculate_eigs_symmetric_matrix(A_size,A,A_eigs)
  
  WRITE(*,*) 'Eigenvalues'
  DO II = 1,A_size
     WRITE(*,'(F16.10)') A_eigs(II)
  END DO

  WRITE(*,*) 'The matrix'
  CALL print_matrix(A_size,A)




END PROGRAM matrix_elements_test


