! This module contains the subroutines used to apply the Real Space 
! Renormalization Group (RSRG) method.
! It is part of the Week 11-12 assignment regarding Renormalization Group.
!
! AUTHOR: Michele Guadagnini - ID 1230663


MODULE RSRG
    USE Debugger 
    IMPLICIT NONE
    
    TYPE Parameters
        INTEGER :: NN, TruncN
        DOUBLE PRECISION :: Lambda, ExitTh
        INTEGER :: MaxItr, kk
        LOGICAL :: Debug
    END TYPE
    
    INTERFACE OPERATOR (.kron.)
        MODULE PROCEDURE kronProduct
    END INTERFACE
    
    CONTAINS
    
!-- PARAMETERS CONFIGURATION -------------------------------------------

    !it reads the parameters from a file
    SUBROUTINE InitParamsFromFile(f_config, Pars, argNames)
        CHARACTER(len=40) :: Line, f_config
        INTEGER :: ii, jj, Npars, ios
        TYPE(Parameters) :: Pars
        CHARACTER(len=8), DIMENSION(:) :: argNames
        CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: argVals
        
        Npars = size(argNames)
        ALLOCATE(argVals(Npars))
        
        WRITE(*,"(A)") "Reading the input parameters from the file: "//f_config
        
        OPEN(unit=43, file=f_config, status='old', IOSTAT=ios)
        CALL CatchError((ios .ne. 0), "Cannot open config file.", Fatal=.TRUE.)
    
        DO ii=1,Npars !assuming number of lines is equal to the number of arguments to read
            ! READING LINE ii
            READ(43,"(A)",IOSTAT=ios) Line
            CALL CatchError((ios .ne. 0), "Cannot read a line from config file. Line:", Fatal=.TRUE., LNumber=ii)
            
            ! SEARCH FOR ALL THE PARAMETERS IN THE LINE, UNTIL ONE IS FOUND
            DO jj=1,Npars 
                IF (INDEX(Line, argNames(jj)) .ne. 0) THEN
                    argVals(jj) = TRIM( Line(8:) )  !cutting out the value from the line
                    EXIT    !exit the loop over jj
                ENDIF
            ENDDO
        ENDDO
        
        READ(argVals, *, IOSTAT=ios) Pars   !save parameters values, checking type
        CALL CatchError((ios .ne. 0), "Invalid value for one or more parameters.", FATAL=.TRUE.)
        CLOSE(43)
        
        CALL InitDebugger(Pars%Debug)
        WRITE(*,"(A)") "Parameters configured to the following values: "
        DO ii=1,Npars
            WRITE(*,"(A,A,A)") "    ", argNames(ii), argVals(ii)
        ENDDO
        WRITE(*,*)
        
        DEALLOCATE(argVals)
        RETURN
    END SUBROUTINE
    
!-- GENERAL MATHEMATICAL TOOLS -----------------------------------------

    !Absolute difference between 2 DOUBLE PRECISION numbers
    FUNCTION AbsDiff(AA, BB) RESULT(Res)
        DOUBLE PRECISION :: AA, BB, Res

        Res = ABS( AA - BB ) 
        
        RETURN
    END FUNCTION
    
    !kronecker delta function
    FUNCTION d_kron(xx, yy) RESULT(dd)
        INTEGER xx,yy
        DOUBLE PRECISION dd

        IF (xx.eq.yy) THEN
            dd = 1d0
        ELSE
            dd = 0d0
        ENDIF
    
        RETURN
    END FUNCTION
    
    !Pauli matrices calculation based on kronecker delta
    FUNCTION PauliMat(xx) RESULT(sigma)
        INTEGER xx      !xx=1 -> sigma_x, xx=2 -> sigma_y ...
        DOUBLE COMPLEX, DIMENSION(2,2) :: sigma
        
        sigma(1,1) = DCMPLX( d_kron(xx,3),           0d0)
        sigma(1,2) = DCMPLX( d_kron(xx,1), -d_kron(xx,2))
        sigma(2,1) = DCMPLX( d_kron(xx,1), +d_kron(xx,2))
        sigma(2,2) = DCMPLX(-d_kron(xx,3),           0d0)
        
        RETURN
    END FUNCTION
    
    !Double Precision Identity Matrix of size ms
    FUNCTION DIdentity(ms) RESULT(IDms)
        INTEGER ms, ll
        DOUBLE PRECISION, DIMENSION(ms,ms) :: IDms
        
        IDms = 0d0
        DO ll = 1,ms
            IDms(ll,ll) = 1d0
        ENDDO
        
        RETURN
    END FUNCTION
    
    !kronecker (direct) product
    FUNCTION kronProduct(Mat1, Mat2) RESULT(Res)
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Mat1
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Mat2
        DOUBLE PRECISION, DIMENSION(SIZE(Mat1,dim=1)*SIZE(Mat2,dim=1),SIZE(Mat1,dim=2)*SIZE(Mat2,dim=2)) :: Res
        INTEGER Nrows1, Ncols1, Nrows2, Ncols2
        INTEGER ii, jj, aa,bb,cc,dd
        
        Nrows1 = SIZE(Mat1, dim=1)
        Ncols1 = SIZE(Mat1, dim=2) 
        Nrows2 = SIZE(Mat2, dim=1)
        Ncols2 = SIZE(Mat2, dim=2)
        
        Res = 0d0
        DO ii = 1,Nrows1
            DO jj = 1,Ncols1
                aa = (ii-1)*Nrows2 +1
                bb = ii*Nrows2
                cc = (jj-1)*Ncols2 +1
                dd = jj*Ncols2
                Res(aa:bb,cc:dd) = Mat1(ii,jj) * Mat2
            ENDDO
        ENDDO
    
        RETURN
    END FUNCTION  
    
    !matrix projection 
    FUNCTION Projection(Proj, Mat) RESULT(Res)
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Proj
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Mat
        DOUBLE PRECISION, DIMENSION(SIZE(Proj,dim=2),SIZE(Proj,dim=2)) :: Res 
        
        Res = MATMUL(TRANSPOSE(Proj), MATMUL(Mat, Proj))        
        
        RETURN
    END FUNCTION
    
    !It computes eigenvalues and eigenvectors of a real symmetric matrix
    SUBROUTINE SymmetricEigenpairs(NN, HH, EigVals)
        INTEGER NN, INFO, LWORK
        DOUBLE PRECISION, DIMENSION(:,:) :: HH
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
        DOUBLE PRECISION, DIMENSION(:) :: EigVals
        
        !Query the optimal workspace
        ALLOCATE(WORK(2))
        LWORK = -1
        CALL DSYEV("V","U", NN, HH, NN, EigVals, WORK, LWORK, INFO)
        CALL CatchError((INFO.ne.0), "    Unable to set optimal workspace. INFO:", .TRUE., INFO)
        
        LWORK = INT(WORK(1), kind=4)
        DEALLOCATE(WORK)
        
        !Solving the eigenproblem
        ALLOCATE(WORK(LWORK))
        CALL DSYEV("V","U", NN, HH, NN, EigVals, WORK, LWORK, INFO)
        CALL CatchError((INFO.ne.0), "    Unsuccessful diagonalization. INFO:", .TRUE., INFO)
        WRITE(*,"(A,I2)") "    Eigenproblem solved with exit status: ", INFO
        
        DEALLOCATE(WORK)
        
        RETURN
    END SUBROUTINE
    
!-- ISING 1D w TRANSVERSE FIELD SPECIFIC SUBROUTINES/FUNCTIONS ---------
    
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
    
    !It prints the first up to 10x10 matrix of eigenvectors
    SUBROUTINE PrintEigVecs(filename, eigvs)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:,:) :: eigvs
        INTEGER ii, kk, jj, maxrows
        
        kk = MIN(size(eigvs, dim=1),10)
        maxrows = MIN(size(eigvs, dim=1),10)
        
        OPEN(unit=10, file=filename, status='unknown')
        DO ii = 1,maxrows
            WRITE(10,'(*(ES16.8))') ( eigvs(ii,jj), jj=1,kk )
        ENDDO
        CLOSE(10)
        
        RETURN
    END SUBROUTINE
    
    !It prints the first kk eigenvalues in column
    SUBROUTINE PrintEigVals(Pars, filename, eigs, DescrSz, lastIT)
        TYPE(Parameters) Pars
        CHARACTER(len=*) filename   !to this is passed a constant, so it cannot be modified
        CHARACTER(len=40) fullpath
        CHARACTER(len=2) Nchar, Tchar
        DOUBLE PRECISION, DIMENSION(:) :: eigs
        INTEGER ii, kset, lastIT
        DOUBLE PRECISION DescrSz
        LOGICAL f_exist
                
        WRITE(Nchar,"(I0.2)") Pars%NN
        WRITE(Tchar,"(I0.2)") Pars%TruncN
        fullpath = "Data/"//TRIM(filename)//"_N"//TRIM(Nchar)//"_T"//TRIM(Tchar)//".dat"
        INQUIRE(file=fullpath, exist=f_exist)   !check if file already exist

        IF (Pars%Debug) THEN
            PRINT *, "Fullname: ", fullpath
            PRINT *, "Exists? : ", f_exist
        ENDIF
        
        OPEN(unit=20 , file=TRIM(fullpath), status='unknown', action='write', position='append')        
        IF (.NOT.f_exist) THEN
            WRITE(20, "(A)")    "# last itr    DescrSz     Lambda      eigenvalues/(# particles - 1) ->"
        ENDIF
        kset = MIN(SIZE(eigs), Pars%kk) ! to avoid errors of exceeded array boundaries in the line below
        WRITE(20,"(I3, ES16.8, *(F15.8))") lastIT, DescrSz, Pars%Lambda, (eigs(ii)/(DescrSz-1d0), ii=1,kset)
        CLOSE(20)
        
        RETURN
    END SUBROUTINE
    
END MODULE
