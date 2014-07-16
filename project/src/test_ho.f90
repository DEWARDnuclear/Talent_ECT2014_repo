PROGRAM test_ho

USE types
USE Ho_basis
USE numerical_integration

IMPLICIT NONE

INTEGER :: gpts
!ho stuff
REAL(kind=r_kind) :: bosc, ol
INTEGER :: n, npr, l, lpr
REAL(kind=r_kind), ALLOCATABLE :: wf1(:), wf2(:)



gpts = 70 !larger than 70 causes numerical errors or floating poit exeptions

CALL init_grid_GH(gpts)

ALLOCATE(wf1(grid_size_GH),wf2(grid_size_GH))

bosc = 1.0_r_kind

n = 50
npr = 50
l = 1
lpr = 1

CALL RadHO(n,l,bosc,grid_points_GH,wf1,grid_size_GH)

CALL RadHO(npr,lpr,bosc,grid_points_GH,wf2,grid_size_GH)

ol  = overlap_ho(n,l,npr,lpr)

WRITE(*,*) 'ol = ', ol


END PROGRAM test_ho




