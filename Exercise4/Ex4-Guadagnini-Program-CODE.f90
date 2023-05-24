! --- Exercise 4 - Automated plots and fits -------------------------- !
! The aim of this program is to use FORTRAN along with Python and 
! gnuplot to do automated plots and fits of the scaling timings of 
! different matrix-matrix multiplication methods.
! FORTRAN is used to run the multiplications and save the data.
!
! AUTHOR: Michele Guadagnini - Mt. 1230663
! -------------------------------------------------------------------- !

PROGRAM MatrixTest
    USE Debugger
    USE LoopMultiply
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: expargs = 4   !# of input command arguments
    CHARACTER(40), DIMENSION(expargs) :: args, filenames    !raw and parsed input arguments
    INTEGER*2 MatSize                   !Matrix size
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AA, BB, CC     !Matrices
    DOUBLE PRECISION start, fin         !Variables to compute exec. time
    INTEGER idx, jj, totargs, ios, iosOut, OpenErrors
    CHARACTER(10) debugstr
    LOGICAL :: Debug = .FALSE.
    
!   (Activate debug) Reading the filenames to read size and print results
    totargs = command_argument_count()
    IF (totargs == expargs+1) THEN
        CALL get_command_argument(totargs, debugstr)
        debugstr = TRIM(debugstr)
        IF (debugstr == "TRUE") THEN
            Debug = .TRUE.
        ENDIF
        totargs = totargs -1
    ENDIF    
    IF (totargs == expargs) THEN
        DO idx = 1,totargs
            CALL get_command_argument(idx, args(idx))
        ENDDO
        READ(args,"(A)",IOSTAT=ios) filenames
    ENDIF
    IF ((ios.ne.0) .or. (totargs /= expargs)) THEN 
        CALL CatchError("cannot read the filenames needed to read the size and store the results!", Fatal=.TRUE., LNumber=41)
    ENDIF
    CALL Checkpoint(Debug, "Filenames read successfully.", 43)
    
!   Reading the size of the matrices from file
    MatSize = ReadSizeFromFile(filenames(1), Debug)
    CALL Checkpoint(Debug, "Matrix size read successfully.", 47)
    
!   Opening the files for the outputs
    OpenErrors = 0
    DO jj=2,4
        OPEN(unit=jj*10, file=filenames(jj), position='append', status='unknown', IOSTAT=iosOut)
        OpenErrors = OpenErrors + ABS(iosOut)
    ENDDO
    IF (OpenErrors.ne.0) THEN
        CALL CatchError("cannot open one or more output files.", Fatal=.TRUE., LNumber=56)
    ENDIF
    CALL Checkpoint(Debug, "Output files successfully opened.", 58)
    
!   Allocate the memory for all matrices
    ALLOCATE(AA(MatSize,MatSize))
    ALLOCATE(BB(MatSize,MatSize))
    ALLOCATE(CC(MatSize,MatSize))
    CALL Checkpoint(Debug, "Memory correctly allocated.", 64)
        
!   Random initialization of the matrices
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(AA)
    CALL RANDOM_NUMBER(BB)
    CALL Checkpoint(Debug, "Matrices initialization complete.", 70)
    
!   First Loop order
    CALL CPU_TIME(start)
    CC = LoopMul1(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "First loop order multiplication successful.", 76)
    CALL WriteExecTimeToFile(20, MatSize, start, fin)

!   Second Loop order
    CALL CPU_TIME(start)
    CC = LoopMul2(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "Second loop order multiplication successful.", 83)
    CALL WriteExecTimeToFile(30, MatSize, start, fin)
    
!   Fortran intrinsic function
    CALL CPU_TIME(start)
    CC = MATMUL(AA, BB)
    CALL CPU_TIME(fin)
    CALL Checkpoint(Debug, "MATMUL multiplication successful.", 90)
    CALL WriteExecTimeToFile(40, MatSize, start, fin)
    
!   Closing the output files
    DO jj = 2,4
        CLOSE(jj*10)
    ENDDO
    CALL Checkpoint(Debug, "Output files closed successfully.", 97)
    
!   Deallocate memory and stop program
    DEALLOCATE(AA)
    DEALLOCATE(BB)
    DEALLOCATE(CC) 
    CALL Checkpoint(Debug, "Execution is complete!")
    STOP
END PROGRAM MatrixTest    
