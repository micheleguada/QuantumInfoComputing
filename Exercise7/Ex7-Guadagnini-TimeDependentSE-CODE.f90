! This module discretizes the space and time, computes the kinetic and
! potential terms of the hamiltonian, computes the time evolution.
! It uses the Split Operators Method exploiting Fourier Transform. 
! It is part of the Week 7 assignment regarding time-dependent
! Schroedinger Equation.
!
! AUTHOR: Michele Guadagnini - ID 1230663


MODULE TimeDependentSE
    USE Debugger
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    
    TYPE Parameters
        INTEGER NdivX, NdivT
        DOUBLE PRECISION LL, TT        
        LOGICAL Debug
    END TYPE
    
    DOUBLE PRECISION, PARAMETER :: Pi = 2d0*DASIN(1d0)
    INTEGER, PARAMETER :: FFTW_ESTIMATE = 64
    
    CONTAINS
    
    
!Subroutines used for input and initialization ---------------------------------------
    
    !Default parameters initialization
    SUBROUTINE InitDefaults(Pars)
        TYPE(Parameters) Pars
        
        Pars%LL = 10d0
        Pars%NdivX = 101
        Pars%NdivT = 101
        Pars%TT = 5d0
        Pars%Debug = .FALSE.
        
        RETURN
    END SUBROUTINE
    
    !it reads the parameters from a file whose name is passed as command line argument
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
            CALL CatchError("Cannot open config file. Using Defaults", Fatal=.FALSE.)
            CALL InitDefaults(Pars)
        ELSE        
            DO ii=1,Npars !assuming number of lines is equal to the number of arguments to read
                READ(43,"(A)",IOSTAT=ios) Line
                IF (ios .ne. 0) THEN
                    WRITE(ni,"(I2)") ii
                    CALL CatchError("Cannot read line "//ni//" from config file.", Fatal=.FALSE.)
                ENDIF
                
                DO jj=1,Npars !search for all the parameters in the line, until one is found
                    IF (INDEX(Line, argNames(jj)) .ne. 0) THEN
                        argVals(jj) = TRIM( Line(8:) )
                        EXIT
                    ENDIF
                ENDDO
            ENDDO
            
            READ(argVals, *, IOSTAT=ios) Pars   !NdivX, NdivT, LL, TT, Debug
            IF (ios .ne. 0) THEN
                CALL CatchError("Invalid value for one or more parameters. Setting them to default.", FATAL=.FALSE.)
                CALL InitDefaults(Pars)
            ENDIF
        ENDIF
        CLOSE(43)
        
        DEALLOCATE(argVals)
        RETURN
    END SUBROUTINE
    
    !Discretization of the space and time
    SUBROUTINE Discretize(Pars, X_vec, T_vec, dx, dt) 
        TYPE(Parameters) Pars
        DOUBLE PRECISION xmin, tmin, dx, dt
        DOUBLE PRECISION, DIMENSION(:) :: X_vec, T_vec
        INTEGER ii
        
        WRITE(*,"(A)") "Discretizing the space..."
        xmin = -Pars%LL/2d0
        dx = Pars%LL/DFLOAT(Pars%NdivX-1) 
        DO ii=1,Pars%NdivX
            X_vec(ii) = xmin + dx*DFLOAT(ii-1)
        ENDDO
        
        WRITE(*,"(A)") "Discretizing the time..."
        tmin = 0d0
        dt = Pars%TT/DFLOAT(Pars%NdivT-1)
        DO ii=1,Pars%NdivT
            T_vec(ii) = tmin + dt*DFLOAT(ii-1)
        ENDDO
        
        RETURN
    END SUBROUTINE
    
    
!Subroutines used mainly to compute the ground state wavefunction --------------------
       
    !It computes the harmonic oscillator potential from the basis of discretization
    SUBROUTINE HarmonicOscillatorPot(Pars, X_vec, v_pot, tt)
        TYPE(Parameters) Pars
        DOUBLE PRECISION, DIMENSION(:) :: X_vec, v_pot
        DOUBLE PRECISION, OPTIONAL :: tt
           
        IF (PRESENT(tt)) THEN
            v_pot = 0.5d0*((X_vec - tt/Pars%TT)**2)
        ELSE
            WRITE(*,"(A)") "Computing the harmonic oscillator potential..."
            v_pot = 0.5d0*(X_vec**2)
        ENDIF
        
        RETURN
    END SUBROUTINE
    
    !Initialization of the matrix
    SUBROUTINE InitHamiltonian(Pars, v_pot, Ham, dx)
        TYPE(Parameters) Pars
        DOUBLE PRECISION dx
        INTEGER ii
        DOUBLE PRECISION, DIMENSION(:,:) :: Ham
        DOUBLE PRECISION, DIMENSION(:) :: v_pot
        
        WRITE(*,"(A)") "Initializing the hamiltonian matrix..."
        
        Ham = 0d0
        DO ii=1,Pars%NdivX
            !main diagonal
            Ham(ii,ii) = 1d0/(dx**2) + v_pot(ii)
            
            !upper diagonal (lower diagonal is not necessary due to LAPACK's subroutine implementation)
            IF (ii .ne. Pars%NdivX) THEN
                Ham(ii,ii+1) = - (0.5d0)/(dx**2)
            ENDIF
        ENDDO   
        
        RETURN
    END SUBROUTINE
    
    !It computes eigenvalues and eigenvectors of a real symmetric matrix
    SUBROUTINE SymmetricEigenpairs(Ndiv, Ham, EigVals)
        INTEGER Ndiv, INFO, LWORK
        DOUBLE PRECISION, DIMENSION(:,:) :: Ham
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
        DOUBLE PRECISION, DIMENSION(:) :: EigVals
        
        WRITE(*,"(A)") "Starting matrix diagonalization..."
        
        !Query the optimal workspace
        ALLOCATE(WORK(2))
        LWORK = -1
        CALL DSYEV("V","U", Ndiv, Ham, Ndiv, EigVals, WORK, LWORK, INFO)
        LWORK = INT(WORK(1), kind=4)
        DEALLOCATE(WORK)
        
        IF (INFO /= 0) THEN
            WRITE(*,"(A,I2)") " Unable to set optimal workspace. INFO: ", INFO
            STOP
        ENDIF
        
        !Solving the eigenproblem
        ALLOCATE(WORK(LWORK))
        CALL DSYEV("V","U", Ndiv, Ham, Ndiv, EigVals, WORK, LWORK, INFO)
        WRITE(*,"(A,I2)") " Eigenproblem solved with exit status: ", INFO
        
        DEALLOCATE(WORK)
        
        RETURN
    END SUBROUTINE
    
    
!Subroutines used to compute the time evolution --------------------------------------

    !It computes momentum eigenvalues to use as momentum space grid on fourier transformed wavefunction
    SUBROUTINE MomentumSpace(Pars, P_vec, dp)
        TYPE(Parameters) Pars
        DOUBLE PRECISION dp
        DOUBLE PRECISION, DIMENSION(:) :: P_vec
        INTEGER ii    
        WRITE(*,"(A)") "Computing the momentum space..."
        dp = 2d0*Pi/Pars%LL
        !The order of the P_vec elements must match the output of FFT
        DO ii=1,Pars%NdivX
            IF (ii .le. (Pars%NdivX/2+1)) THEN
                P_vec(ii) = dp*DFLOAT(ii-1)
            ELSE
                P_vec(ii) = dp*DFLOAT(ii -1 - Pars%NdivX)             
            ENDIF
        ENDDO    
        RETURN
    END SUBROUTINE
    
    !It computes the forward Fourier transform of a complex vector
    SUBROUTINE ComplexFFT(Pars, Cvec, flag)
        TYPE(Parameters) Pars
        DOUBLE COMPLEX, DIMENSION(:) :: Cvec
        DOUBLE COMPLEX, DIMENSION(size(Cvec)) :: temp
        TYPE(C_PTR) plan
        INTEGER flag    !-1 for forward transf., +1 for backward
        CALL dfftw_plan_dft_1d(plan, Pars%NdivX, Cvec, temp, flag, FFTW_ESTIMATE)
        CALL dfftw_execute_dft(plan, Cvec, temp)
        CALL dfftw_destroy_plan(plan) 
        Cvec = temp    
        RETURN
    END SUBROUTINE
    
    !it computes the L2 norm of the wavefunction
    FUNCTION L2Norm(Psi, ds) RESULT(Norm)
        DOUBLE COMPLEX, DIMENSION(:) :: Psi
        DOUBLE PRECISION :: Norm, ds
        
        Norm = DSQRT( SUM(DREAL(Psi)**2 + DIMAG(Psi)**2)*ds )
        
        RETURN
    END FUNCTION
    
    
!Additional subroutines for printings ------------------------------------------------
    
    !Store eigenfunctions for each timestep
    SUBROUTINE PrintTimeEvol(output, tstep, Psi, X_vec, TT_str)
        CHARACTER(len=40) :: output, fullname, TT_str
        CHARACTER(len=5) :: tstep_str
        DOUBLE COMPLEX, DIMENSION(:) :: Psi
        DOUBLE PRECISION, DIMENSION(:) :: X_vec
        INTEGER :: tstep, ii
        
        CALL system("mkdir -p TimeEvol"//TRIM(TT_str))
        
        WRITE(tstep_str, "(I0.4)") tstep
        fullname = "TimeEvol"//TRIM(TT_str)//"/"//TRIM(output)//TRIM(tstep_str)//".dat"
        
        OPEN(unit=30, file=fullname, status='unknown')
        WRITE(30,"(A)") "# basis            wavefunction"
        DO ii = 1,size(X_vec)
            WRITE(30,'(*(ES16.8))') X_vec(ii), DSQRT(DREAL(Psi(ii))**2 + DIMAG(Psi(ii))**2)
        ENDDO
        CLOSE(30)
    
        RETURN
    END SUBROUTINE  
    
    SUBROUTINE PrintPotEvol(output, tstep, pot, X_vec)
        CHARACTER(len=40) :: output, fullname
        CHARACTER(len=5) :: tstep_str
        DOUBLE PRECISION, DIMENSION(:) :: pot
        DOUBLE PRECISION, DIMENSION(:) :: X_vec
        INTEGER :: tstep, ii
        
        WRITE(tstep_str, "(I0.4)") tstep
        fullname = "PotEvol/"//TRIM(output)//TRIM(tstep_str)//".dat"
        
        OPEN(unit=30, file=fullname, status='unknown')
        WRITE(30,"(A)") "# basis            Pot"
        DO ii = 1,size(X_vec)
            WRITE(30,'(*(ES16.8))') X_vec(ii), pot(ii)
        ENDDO
        CLOSE(30)
    
        RETURN
    END SUBROUTINE 
    
    !It prints the first k eigenfunctions in parallel columns
    SUBROUTINE PrintEigVecs(filename, v_bas, eigvs, kk, dx)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:,:) :: eigvs
        DOUBLE PRECISION, DIMENSION(:) :: v_bas
        INTEGER ii, kk, jj
        DOUBLE PRECISION dx
        
        OPEN(unit=10, file=filename, status='unknown')
        WRITE(10,"(A)") "# basis           -> k-esim eigenfunctions"
        DO ii = 1,size(eigvs, dim=1)
            !eigenfunctions are normalized to unit norm by dividing it by 1/sqrt(dx)
            WRITE(10,'(*(ES16.8))') v_bas(ii), ( eigvs(ii,jj)/DSQRT(dx), jj=1,kk )
        ENDDO
        CLOSE(10)
        
        RETURN
    END SUBROUTINE
    
    !It prints the first k eigenvalues in column
    SUBROUTINE PrintEigVals(filename, eigs, kk)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:) :: eigs
        INTEGER ii, kk
        
        OPEN(unit=20 , file=filename, status='unknown')
        WRITE(20, "(A)") "# Expected | Computed eigenvalues"
        DO ii = 1,kk
            WRITE(20,"(F5.1, F15.8)") (DFLOAT(ii)-0.5d0), eigs(ii)
        ENDDO
        CLOSE(20)
        
        RETURN
    END SUBROUTINE
    
END MODULE
