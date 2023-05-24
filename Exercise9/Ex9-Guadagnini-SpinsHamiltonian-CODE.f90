! This module contains the subroutines to compute and diagonalize
! the Hamiltonian of N spin particles with nearest-neighbor
! interaction. 
! It is part of the Week 9-10 assignment regarding 1D Ising model.
!
! AUTHOR: Michele Guadagnini - ID 1230663


MODULE SpinsHamiltonian
    USE Debugger 
    IMPLICIT NONE
    
    TYPE Parameters
        INTEGER :: NN, kk
        DOUBLE PRECISION :: Lambda
        LOGICAL :: Debug
    END TYPE
    
    CONTAINS
    
    !it reads the parameters from a file
    SUBROUTINE InitParamsFromFile(f_config, Pars, argNames)
        CHARACTER(len=40) :: Line, f_config
        CHARACTER(len=2) ni
        INTEGER :: ii, jj, Npars, ios
        TYPE(Parameters) :: Pars
        CHARACTER(len=8), DIMENSION(:) :: argNames
        CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: argVals
        
        Npars = size(argNames)
        ALLOCATE(argVals(Npars))
        
        WRITE(*,"(A)") "Reading the input parameters from the file: "//f_config
        
        OPEN(unit=43, file=f_config, status='old',IOSTAT=ios)
        IF (ios .ne. 0) THEN
            CALL CatchError("Cannot open config file.", Fatal=.TRUE.)
        ELSE        
            DO ii=1,Npars !assuming number of lines is equal to the number of arguments to read
                READ(43,"(A)",IOSTAT=ios) Line
                IF (ios .ne. 0) THEN
                    WRITE(ni,"(I2)") ii
                    CALL CatchError("Cannot read line # "//ni//" from config file.", Fatal=.TRUE.)
                ENDIF
                
                DO jj=1,Npars !search for all the parameters in the line, until one is found
                    IF (INDEX(Line, argNames(jj)) .ne. 0) THEN
                        argVals(jj) = TRIM( Line(8:) )
                        EXIT
                    ENDIF
                ENDDO
            ENDDO
            
            READ(argVals, *, IOSTAT=ios) Pars   !look a the definition of type of Pars
            IF (ios .ne. 0) THEN
                CALL CatchError("Invalid value for one or more parameters.", FATAL=.TRUE.)
            ENDIF
        ENDIF
        CLOSE(43)
        
        DEALLOCATE(argVals)
        RETURN
    END SUBROUTINE
    
    !it computes the field term of the hamiltonian
    FUNCTION Ham_Zsum(Pars) RESULT(Mat_Zsum)
        TYPE(Parameters) Pars
        INTEGER aa, ii
        DOUBLE PRECISION, DIMENSION(2**Pars%NN, 2**Pars%NN) :: Mat_Zsum
        DOUBLE PRECISION temp
        
        ! Z contribution is diagonal, so it can be computed in this way
        Mat_Zsum = 0d0
        DO aa = 1, 2**Pars%NN
            Mat_Zsum(aa,aa) = DFLOAT(Pars%NN)
            DO ii = 1,Pars%NN
                temp = 2d0*DFLOAT( MOD((aa-1)/(2**(ii-1)), 2) )
                Mat_Zsum(aa,aa) = Mat_Zsum(aa,aa) - temp
            ENDDO
        ENDDO
        
        RETURN
    END FUNCTION
    
    !it computes the interacting term of the hamiltonian
    FUNCTION Ham_Xsum(Pars) RESULT(Mat_Xsum)
        TYPE(Parameters) Pars
        INTEGER aa, bb, ii
        DOUBLE PRECISION, DIMENSION(2**Pars%NN, 2**Pars%NN) :: Mat_Xsum
        
        Mat_Xsum = 0d0
        DO aa = 1, 2**Pars%NN
            DO bb = 1, 2**Pars%NN
                DO ii = 1,Pars%NN-1
                    IF ( IEOR(aa-1,bb-1) .eq. (2**(ii-1) + 2**ii) ) THEN
                        Mat_Xsum(aa,bb) = Mat_Xsum(aa,bb) + 1d0
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        
        RETURN
    END FUNCTION
    
    !It computes eigenvalues and eigenvectors of a real symmetric matrix
    SUBROUTINE SymmetricEigenpairs(NN, Ham, EigVals)
        INTEGER NN, INFO, LWORK
        DOUBLE PRECISION, DIMENSION(:,:) :: Ham
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
        DOUBLE PRECISION, DIMENSION(:) :: EigVals
        
        WRITE(*,"(A)") "Starting matrix diagonalization..."
        
        !Query the optimal workspace
        ALLOCATE(WORK(2))
        LWORK = -1
        CALL DSYEV("V","U", NN, Ham, NN, EigVals, WORK, LWORK, INFO)
        LWORK = INT(WORK(1), kind=4)
        DEALLOCATE(WORK)
        
        IF (INFO /= 0) THEN
            WRITE(*,"(A,I2)") " Unable to set optimal workspace. INFO: ", INFO
            STOP
        ENDIF
        
        !Solving the eigenproblem
        ALLOCATE(WORK(LWORK))
        CALL DSYEV("V","U", NN, Ham, NN, EigVals, WORK, LWORK, INFO)
        WRITE(*,"(A,I2)") " Eigenproblem solved with exit status: ", INFO
        
        DEALLOCATE(WORK)
        
        RETURN
    END SUBROUTINE

    !It prints the first k eigenfunctions in parallel columns
    SUBROUTINE PrintEigVecs(filename, eigvs, kk)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:,:) :: eigvs
        INTEGER ii, kk, jj
        
        OPEN(unit=10, file=filename, status='unknown')
        !WRITE(10,"(A)") "# k-esim eigenfunction "
        DO ii = 1,size(eigvs, dim=1)
            WRITE(10,'(*(ES16.8))') ( eigvs(ii,jj), jj=1,kk )
        ENDDO
        CLOSE(10)
        
        RETURN
    END SUBROUTINE
    
    !It prints the first k eigenvalues in column
    SUBROUTINE PrintEigVals(Pars, filename, eigs)
        TYPE(Parameters) Pars
        CHARACTER(len=*) filename
        CHARACTER(len=2) Nchar
        DOUBLE PRECISION, DIMENSION(:) :: eigs
        INTEGER ii, kset
        LOGICAL f_exist
                
        WRITE(Nchar,"(I0.2)") Pars%NN
        filename = "Data/"//TRIM(filename)//"_N"//TRIM(Nchar)//".dat"
        
        INQUIRE(file=filename, exist=f_exist)
        
        IF (Pars%Debug) THEN
            PRINT *, f_exist
            PRINT *, filename
        ENDIF
        
        OPEN(unit=20 , file=filename, status='unknown', action='write', position='append')        
        IF (.NOT.f_exist) THEN
            WRITE(20, "(A)") "#     Lambda      eigenvalues ->"
        ENDIF
        
        kset = MIN(SIZE(eigs), Pars%kk) ! to avoid errors of exceeded array boundaries in the line below
        WRITE(20,"(*(F15.8))") Pars%Lambda, ( eigs(ii), ii=1,kset)
        CLOSE(20)
        
        RETURN
    END SUBROUTINE   
    
END MODULE
