!This module contains the matrix-matrix multiplication functions, the function to read the size from file,
! a subroutine to print the matrices and a subroutine to print the resulting exec. timings in a file
MODULE LoopMultiply
    USE Debugger
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
    
!   Read matrix size from file, check that it is a strictly positive integer and returns it. 
!   If the reading is not successful it calls CatchError subroutine
    FUNCTION ReadSizeFromFile(FileToRead, Debug) RESULT(MatSize) 
        CHARACTER(len=*) FileToRead
        INTEGER*2 MatSize, iosOpen, iosRead
        LOGICAL Debug
        
        FileToRead = TRIM(FileToRead)
        OPEN(unit=10, file=FileToRead, status='old', IOSTAT=iosOpen)
        IF (iosOpen.ne.0) THEN
            CALL CatchError("cannot open the file to read the matrix size!", Fatal=.TRUE., LNumber=54)
        ENDIF
        READ(10, *, iostat=iosRead) MatSize  
        IF (iosRead.ne.0) THEN
            CALL CatchError("cannot read the size value!", Fatal=.TRUE., LNumber=58)
        ENDIF        
        CLOSE(10)
        
        IF (MatSize.lt.1) THEN 
            CALL CatchError("matrix size must be a strictly positive integer number!", Fatal=.TRUE., LNumber=63)
            RETURN
        ELSE
            CALL Checkpoint(Debug, "Size of matrices read successfully.", 66)
        ENDIF
        
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
    
!   Check and print execution time
    SUBROUTINE WriteExecTimeToFile(FileToWriteID, MatSize, start, fin)
        DOUBLE PRECISION start, fin
        INTEGER FileToWriteID
        INTEGER*2 MatSize
        
        !Check for delta time to be positive
        IF (start .gt. fin) THEN
            CALL CatchError("Error: execution time must be positive!", Fatal=.FALSE., LNumber=98)
            !The output is set to -1000 to emphasize it is wrong
            WRITE(FileToWriteID,"(I6,A,f10.6)") MatSize, " ", -1000.0
        ELSE
            WRITE(FileToWriteID,"(I6,A,f10.6)") MatSize, " ", fin-start
        ENDIF
               
        RETURN
    END SUBROUTINE

END MODULE
