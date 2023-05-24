! --- Exercise 3 - Matrix-Matrix multiplication ---------------------- !
! The aim of this program is to test three different matrix-matrix 
! multiplication methods and check their performances also with 
! different compiler flags.
! This exercise is meant also for practicing with checkpoints,
! debugging and error handling.
!
! AUTHOR: Michele Guadagnini - Mt. 1230663
! -------------------------------------------------------------------- !

MODULE LoopMultiply
    IMPLICIT NONE
    CONTAINS
    
!   Matrix multiplication with 1st loop order 
    FUNCTION LoopMul1(AA, BB) RESULT(CC)
        INTEGER*2 ii, jj, kk
        DOUBLE PRECISION, DIMENSION(:,:) :: AA
        DOUBLE PRECISION, DIMENSION(:,:) :: BB
        DOUBLE PRECISION, DIMENSION(size(AA,dim=1),size(BB,dim=2)) :: CC
        
        CC = 0.0
        DO ii = 1,size(AA,dim=1)
            DO jj = 1,size(BB,dim=2)
                DO kk = 1,size(AA,dim=2)
                    CC(ii,jj) = CC(ii,jj) + AA(ii,kk)*BB(kk,jj)
                ENDDO
            ENDDO
        ENDDO            
        RETURN
    END FUNCTION
    
!   Matrix multiplication with 2nd loop order 
    FUNCTION LoopMul2(AA, BB) RESULT(CC)
        INTEGER*2 ii, jj, kk
        DOUBLE PRECISION, DIMENSION(:,:) :: AA
        DOUBLE PRECISION, DIMENSION(:,:) :: BB
        DOUBLE PRECISION, DIMENSION(size(AA,dim=1),size(BB,dim=2)) :: CC
        
        CC = 0.0
        DO jj = 1,size(BB,dim=2)
            DO kk = 1,size(AA,dim=2)
                DO ii = 1,size(AA,dim=1)
                    CC(ii,jj) = CC(ii,jj) + AA(ii,kk)*BB(kk,jj)
                ENDDO
            ENDDO
        ENDDO            
        RETURN
    END FUNCTION
    
!   Subroutine for printing the matrices
    SUBROUTINE PrintMatrix(DD, string)
        INTEGER*2 :: ll
        CHARACTER(len=*) string
        DOUBLE PRECISION, DIMENSION(:,:) :: DD
        
        !For too large matrices (more than 36 elements) it prints the shape instead of all elements
        IF ((size(DD,dim=1)*size(DD,dim=2)).ge.(6*6)) THEN
            WRITE(*,fmt="(A,A,'(',I6,',',I6,')')") string, " with shape: ", shape(DD)
        ELSE 
            WRITE(*,"(A)") string
            DO ll = 1,size(DD,dim=2)
                WRITE(*,fmt="(*(f6.3))") DD(:,ll)
            ENDDO
        ENDIF
        
        RETURN
    END SUBROUTINE
    
!   Default initialization subroutine
    SUBROUTINE DefaultInit(MatSizes)
        INTEGER*2, DIMENSION(4) :: MatSizes, default
        
        !Default initialization
        default = 100      !vector of length 4
        WRITE(*,*) "The possible arguments are: # rows A, # cols A, # rows B, # cols B, (logical debug flag)"
        WRITE(*,fmt="(A,I6,A)") "Setting default values: size ", default(1), " square matrices with debug inactive"
        
        MatSizes = default
        
        RETURN
    END SUBROUTINE
    
!   Check shapes of matrices
    SUBROUTINE CheckShapes(MatSizes)
        INTEGER*2, DIMENSION(4) :: MatSizes
        INTEGER ii
        
        !Check for non-positive values in shapes
        DO ii = 1,size(MatSizes)
            IF (MatSizes(ii) .lt. 1) THEN
                WRITE(*,"(A)") "Error: negative or null matrix size provided"
                CALL DefaultInit(MatSizes)
                RETURN
            ENDIF
        ENDDO
        
        !Check for shapes inconsistence
        IF (MatSizes(2) /= MatSizes(3)) THEN
            WRITE(*,"(A)") "Error: inconsistent matrix shapes for matrix-matrix multiplication"
            CALL DefaultInit(MatSizes)
            RETURN
        ENDIF
        
        RETURN
    END SUBROUTINE
    
!   Check and print execution time
    SUBROUTINE PrintExecTime(start,fin)
        DOUBLE PRECISION start, fin
        
        !Check for delta time to be positive
        IF (start .gt. fin) THEN
            WRITE(*,*) "Error: execution time must be positive!"
        ENDIF
        
        WRITE(*,"(A,f10.6)") "Exec. time [s]: ", fin-start
        
        RETURN
    END SUBROUTINE

END MODULE

PROGRAM MatrixTest
    USE LoopMultiply
    USE Checkpoints
    IMPLICIT NONE
    
    INTEGER*2, DIMENSION(4) :: MatSizes
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AA, BB, CC
    DOUBLE PRECISION start, fin
    CHARACTER(8), DIMENSION(5) :: args
    INTEGER idx, totargs, expargs, ios
    LOGICAL :: Debug = .FALSE.
    
!   Checking # of input values and their consistence with variables type
    expargs = 4
    totargs = command_argument_count()
    IF (totargs == expargs) THEN
        DO idx = 1,totargs
            CALL get_command_argument(idx, args(idx))
        ENDDO
        READ(args,*,IOSTAT=ios) MatSizes
        IF (ios.ne.0) THEN  !if conversion of input fails, default initialization is called
            WRITE(*,"(A)") "Error: invalid input value provided"
            CALL DefaultInit(MatSizes) 
        ENDIF
    ELSEIF (totargs == expargs+1) THEN
        DO idx = 1,totargs
            CALL get_command_argument(idx, args(idx))
        ENDDO
        READ(args,*,IOSTAT=ios) MatSizes,Debug
        IF (ios.ne.0) THEN
            WRITE(*,"(A)") "Error: invalid input value provided"
            CALL DefaultInit(MatSizes)
        ELSE
            CALL Checkpoint(Debug, "Input arguments accepted!", 157)
        ENDIF
    ELSEIF (totargs == 0) THEN
        CALL DefaultInit(Matsizes)
    ELSE    
        WRITE(*,fmt="(A,I2,A,I2,A)") "Error: ", expargs," arguments expected, but ", totargs, " provided"
        CALL DefaultInit(MatSizes)
    ENDIF    
    
!   Checking matrix dimensions to be positive and consistent with matrix-matrix multiplication
    CALL CheckShapes(MatSizes)
    CALL Checkpoint(Debug, "Matrices has consistent shapes.", 168)
    
!   Allocate the memory for all matrices
    ALLOCATE(AA(MatSizes(1),MatSizes(2)))
    ALLOCATE(BB(MatSizes(3),MatSizes(4)))
    ALLOCATE(CC(MatSizes(1),MatSizes(4)))
    CALL Checkpoint(Debug, "Memory correctly allocated.", 174)
        
!   Random initialization of the matrices
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(AA)
    CALL RANDOM_NUMBER(BB)
    CALL Checkpoint(Debug, "Matrices initialization complete.", 180)
    
!   Printing initial matrices
    CALL PrintMatrix(AA, "Matrix A")
    CALL PrintMatrix(BB, "Matrix B")
    
!   First Loop order
    CALL CPU_TIME(start)
    CC = LoopMul1(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "First loop order multiplication successful.", 190)
    CALL PrintMatrix(CC, "Matrix C (1st loop order)")
    CALL PrintExecTime(start,fin)

!   Second Loop order
    CALL CPU_TIME(start)
    CC = LoopMul2(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "Second loop order multiplication successful.", 198)
    CALL PrintMatrix(CC, "Matrix C (2nd loop order)")
    CALL PrintExecTime(start,fin)
    
!   Fortran intrinsic function
    CALL CPU_TIME(start)
    CC = MATMUL(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "MATMUL multiplication successful.", 206)
    CALL PrintMatrix(CC, "Matrix C (MATMUL function)")
    CALL PrintExecTime(start,fin)
    
!   Deallocate memory and stop program
    DEALLOCATE(AA)
    DEALLOCATE(BB)
    DEALLOCATE(CC) 
    CALL Checkpoint(Debug, "Execution is complete!")
    STOP
END PROGRAM MatrixTest
    
