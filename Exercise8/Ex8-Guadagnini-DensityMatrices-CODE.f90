! This module contains the subroutines used for allocate the pure 
! wavefunction and to compute the density matrices.
! It is part of the Week 8 assignment regarding Density Matrices
!
! AUTHOR: Michele Guadagnini - ID 1230663


MODULE DensityMatrices
    USE Debugger 
    IMPLICIT NONE
    
    TYPE Parameters
        INTEGER :: Hdim, Nsys
        LOGICAL :: Sep, Debug
    END TYPE
    
    TYPE PureState
        INTEGER :: DD, NN
        LOGICAL :: separable
        DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: WF      !wavefunction  
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
            CALL CatchError("Cannot open config file. Using Defaults", Fatal=.TRUE.)
        ELSE        
            DO ii=1,Npars !assuming number of lines is equal to the number of arguments to read
                READ(43,"(A)",IOSTAT=ios) Line
                IF (ios .ne. 0) THEN
                    WRITE(ni,"(I2)") ii
                    CALL CatchError("Cannot read line "//ni//" from config file.", Fatal=.TRUE.)
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
                CALL CatchError("Invalid value for one or more parameters. Stopping...", FATAL=.TRUE.)
            ENDIF
        ENDIF
        CLOSE(43)
        
        DEALLOCATE(argVals)
        RETURN
    END SUBROUTINE
    
    !It allocates the memory depending on the type of state (separable or not)
    SUBROUTINE AllocPurePsi(Pars, psi)
        TYPE(Parameters) Pars
        TYPE(PureState) psi
        
        psi%separable = Pars%Sep
        psi%DD = Pars%Hdim
        psi%NN = Pars%Nsys
        
        IF (Pars%Sep) THEN
            IF (.NOT.(ALLOCATED(psi%WF))) THEN
                ALLOCATE(psi%WF(psi%DD*psi%NN))  !dim = DD*NN
            ENDIF
        ELSE
            IF (.NOT.(ALLOCATED(psi%WF))) THEN
                ALLOCATE(psi%WF(psi%DD**psi%NN)) !dim = DD^NN
            ENDIF
        ENDIF
        
        RETURN
    END SUBROUTINE

    !It initialize the psi at random
    SUBROUTINE RandomPsi(psi)
        TYPE(PureState) psi
        DOUBLE PRECISION, DIMENSION(SIZE(psi%WF)) :: Re, Im
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(Re)
        CALL RANDOM_NUMBER(Im)
        psi%WF = DCMPLX(Re,Im)        
        
        !normalization if it is the total wavefunction
        IF (.NOT.psi%separable) THEN
            psi%WF = psi%WF / L2Norm(psi%WF)
        ENDIF
        
        RETURN
    END SUBROUTINE
    
    !It computes the density matrix from psi
    SUBROUTINE PureDensityMat(psi, Dmat)
        TYPE(PureState) psi
        DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: Dmat
        DOUBLE COMPLEX, DIMENSION(psi%DD**psi%NN) :: temp
        INTEGER :: MatSize, kk, iidx, jidx
        
        MatSize = psi%DD ** psi%NN
        
        IF (.NOT.(ALLOCATED(Dmat))) THEN
            ALLOCATE(Dmat(MatSize,MatSize))
        ENDIF
        
        !Distinguish between SEPARABLE and NOT SEPARABLE cases 
        IF (psi%separable) THEN            
            DO kk=1,MatSize
                !computing the indeces i and j from k 
                jidx = psi%DD +1 + MOD((kk-1),psi%DD)
                iidx = (kk-1)/psi%DD +1
                temp(kk) = psi%WF(iidx)*psi%WF(jidx)
            ENDDO
            temp = temp / L2Norm(temp) !normalization
            DO kk=1,MatSize
                Dmat(kk,:) = temp(kk)*DCONJG(temp(:))
            ENDDO
        ELSE    !Not separable case
            DO kk=1,MatSize
                Dmat(kk,:) = psi%WF(kk)*DCONJG(psi%WF(:))
            ENDDO            
        ENDIF
        
        !Check normalization
        IF (ABS(DREAL(CTrace(Dmat)) - 1d0) .gt. 1d-14)  THEN   !Tr(rho) -1 should be 0
            WRITE(*,*) CTrace(Dmat)
            CALL CatchError("The resulting Density Matrix does NOT respect normalization condition.", Fatal=.FALSE.)
        ENDIF

        RETURN
    END SUBROUTINE
    
    !It computes the reduced density matrix left or right
    SUBROUTINE RedDensityMat(Pars, DMat, RedDMat, left)
        TYPE(Parameters) Pars
        DOUBLE COMPLEX, DIMENSION(:,:) :: DMat
        DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: RedDMat
        DOUBLE COMPLEX temp
        LOGICAL left     !If it is true, it computes rho_A from rho_AB, H_AB = H_A x H_B
        INTEGER RedMatSize, ii, jj, ki, kj, tt
        
        RedMatSize = Pars%Hdim ** (Pars%Nsys-1)
        IF (.NOT.ALLOCATED(RedDMat)) THEN
            ALLOCATE(RedDMat(RedMatSize,RedMatSize))
        ENDIF
        
        !WARNING: supposing to have N=2
        IF (left) THEN      
        !Tracing out the right system (i.e. computing rho_A from rho_AB)
            DO ii=1,Pars%Hdim
                ki = 1 + Pars%Hdim*(ii-1)
                DO jj=1,Pars%Hdim
                    kj = 1 + Pars%Hdim*(jj-1)
                    temp = DCMPLX(0d0,0d0)
                    DO tt=0,Pars%Hdim-1
                        temp = temp + DMat(ki+tt,kj+tt)
                    ENDDO
                    RedDMat(ii,jj) = temp
                ENDDO
            ENDDO
        ELSE
        !Tracing out the left system (i.e. computing rho_B from rho_AB)
            DO ii=1,Pars%Hdim
                DO jj=1,Pars%Hdim
                    temp = DCMPLX(0d0,0d0)
                    DO tt=0,Pars%Hdim-1
                        temp = temp + DMat(ii+tt*Pars%Hdim,jj+tt*Pars%Hdim)
                    ENDDO
                    RedDMat(ii,jj) = temp
                ENDDO
            ENDDO                    
        ENDIF    
        
        RETURN
    END SUBROUTINE
    
    !It computes the norm of the wavefunction as complex variable (to avoid type conversions from CMPLX to REAL)
    FUNCTION L2Norm(WF) RESULT(norm)
        DOUBLE COMPLEX, DIMENSION(:) :: WF
        DOUBLE COMPLEX norm
        
        norm = ZSQRT(SUM(WF*DCONJG(WF)))
        
        RETURN
    END FUNCTION
    
    !It computes the trace of a complex matrix
    FUNCTION CTrace(Mat) RESULT(trac)
        DOUBLE COMPLEX, DIMENSION(:,:) :: Mat
        DOUBLE COMPLEX trac
        INTEGER jj
        
        !TODO check that matrix is square
        
        trac = DCMPLX(0d0,0d0)
        DO jj=1,SIZE(Mat, dim=1)
            trac = trac + Mat(jj,jj)
        ENDDO
        
        RETURN
    END FUNCTION
    
    SUBROUTINE WriteOnFile(CMat, Outfile)
        DOUBLE COMPLEX, DIMENSION(:,:) :: CMat
        CHARACTER(len=*) Outfile
        INTEGER ll
        
        PRINT *, "Writing Density Matrix on file..."
        
        OPEN(unit=23, file=Outfile, status='unknown')
        WRITE(23,"(A,('('sf6.3xspf6.3x'i)'))") "MATRIX TRACE: ", CTrace(CMat)
        WRITE(23,*) NEW_LINE('A'), "MATRIX ELEMENTS: "
        DO ll = 1,SIZE(CMat,dim=2)
            WRITE(23, "(*('('sf6.3xspf6.3x'i)':x))") CMat(:,ll)
        ENDDO
        CLOSE(23)
        
        RETURN
    END SUBROUTINE
    
END MODULE
