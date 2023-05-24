! --- Exercise 7 - Time Dependent Schroedinger Equation -------------- !
! The aim of this program is to compute the time evolution of the ground 
! state of a time-dependent harmonic oscillator system. 
! In order to solve the time-dependent Schroedinger Equation with the 
! Hamiltonian: H = ( p^2 + (q - q(t))^2 )/2 the Split Operator Method 
! is applied.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !

PROGRAM SchroedingerEq
    USE Debugger
    USE TimeDependentSE
    IMPLICIT NONE
    
    TYPE(Parameters) Pars
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X_vec, T_vec, P_vec, v_pot
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Hamiltonian
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EigVals          !eigenvalues
    DOUBLE PRECISION :: dx, dt, dp
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Psi_0
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: Psi_t, FFTtest
    INTEGER tstep
    
    CHARACTER(len=40) :: EigValsFile, EigVecsFile
    CHARACTER(len=40) :: f_config, output, TT_str
    CHARACTER(len=8), DIMENSION(5) :: argNames
    
    INTEGER ii
    
    !   Results file names (only for debugging)
    EigValsFile = "Eigenvalues.dat"
    EigVecsFile = "Eigenfunctions.dat"
    
!-- READING THE ARGUMENTS FROM FILE --------------------------------------------------
    
!   Names of the input parameters (in the same order as defined in the TYPE(Parameters) definition)
    argNames = (/"NdivX =", "NdivT =", "LL    =", "TT    =", "Debug ="/)
    
!   Parameters file name
    f_config = "Config_Pars.txt"
    
!   Partial name of the file containing results of time evolution
    output = "Tstep"

    CALL InitParamsFromFile(f_config, Pars, argNames)
    
    CALL InitDebugger(Pars%Debug)
    CALL Checkpoint(Pars%Debug, "Input parameters set.")
    
    WRITE(*,"(A)") "Parameters configured to the following values: "
    WRITE(*,"(A,F8.3)") "    LL    = ", Pars%LL
    WRITE(*,"(A,F8.3)") "    TT    = ", Pars%TT
    WRITE(*,"(A,I8)")   "    NdivX = ", Pars%NdivX
    WRITE(*,"(A,I8)")   "    NdivT = ", Pars%NdivT
    WRITE(*,*)
    
    WRITE(TT_str,"(F0.2)") Pars%TT
    
!   Memory allocation
    ALLOCATE( X_vec(Pars%NdivX), v_pot(Pars%NdivX), EigVals(Pars%NdivX) )
    ALLOCATE( T_vec(Pars%NdivT) )
    ALLOCATE( Hamiltonian(Pars%NdivX,Pars%NdivX) )
    CALL Checkpoint(Pars%Debug, "Memory for ground state computation correctly allocated.")
    
!-- COMPUTING THE GROUND STATE AT TIME 0 ---------------------------------------------

    WRITE(*,"(A)") "Starting to compute the ground state eigenfunction at t=0 ..."

!   Discretization
    CALL Discretize(Pars, X_vec, T_vec, dx, dt)
    CALL Checkpoint(Pars%Debug, "Discretization of the space and time successful.")
    
!   Compute potential
    CALL HarmonicOscillatorPot(Pars, X_vec, v_pot)
    CALL Checkpoint(Pars%Debug, "Hamiltonian Potential computed successfully.")
    
    !it prints the discretization basis and potential if debug is active
    IF (Pars%Debug) THEN        
        OPEN(unit=80 , file="test.txt", status='unknown')
        WRITE(80, "(A)") "# basis           potential"
        DO ii = 1,Pars%NdivX
            WRITE(80,"(2ES16.8)") X_vec(ii), v_pot(ii)
        ENDDO
        CLOSE(80)
    ENDIF
    
!   Compute the hamiltonian matrix
    CALL InitHamiltonian(Pars, v_pot, Hamiltonian, dx)
    CALL Checkpoint(Pars%Debug, "Hamiltonian matrix computed successfully.")
    
!   Diagonalize the matrix
    CALL SymmetricEigenpairs(Pars%NdivX, Hamiltonian, EigVals)
    CALL Checkpoint(Pars%Debug, "Matrix diagonalization completed.")
    
!   Storing some results on file for debugging
    IF (Pars%Debug) THEN
        CALL PrintEigVals(EigValsFile, EigVals, 4)
        CALL PrintEigVecs(EigVecsFile, X_vec, Hamiltonian, 4, dx)
    ENDIF
    
!   Storing on memory the ground state eigenfunction normalized
    ALLOCATE(Psi_0(size(X_vec)), Psi_t(size(X_vec)))
    Psi_0 = Hamiltonian(:,1)/DSQRT(dx)
    CALL Checkpoint(Pars%Debug, "Ground state computed and stored.")
    
!   Test the fourier transform subroutines
    IF (Pars%Debug) THEN
        ALLOCATE(FFTtest(Pars%NdivX))
        OPEN(unit=67, file="FFTtest.txt", status="unknown")
        FFTtest = DCMPLX(Psi_0, 0d0)
        WRITE(67,*) "# FOURIER TRANSFORM SUBROUTINES TEST #"
        WRITE(67,*) "##################### Momentum representation #######################"
        WRITE(67,*) "# discrete space  |  Psi  |  F(Psi) "
        CALL ComplexFFT(Pars, FFTtest, -1)
        DO ii=1,size(X_vec)
            WRITE(67,'(*(ES16.8))') X_vec(ii), Psi_0(ii), DSQRT(DREAL(FFTtest(ii))**2+DIMAG(FFTtest(ii))**2)
        ENDDO
        WRITE(67,*) "##################### Space representation    #######################"
        WRITE(67,*) "# discrete space  |  Psi  |  F^(-1)[F(Psi)] "
        CALL ComplexFFT(Pars, FFTtest, +1)
        FFTtest = FFTtest/DCMPLX(Pars%NdivX, 0d0)
        DO ii=1,size(X_vec)
            WRITE(67,'(*(ES16.8))') X_vec(ii), Psi_0(ii), DSQRT(DREAL(FFTtest(ii))**2+DIMAG(FFTtest(ii))**2)
        ENDDO
        WRITE(67,*) ""
        WRITE(67,*) "##################### Norms computation       #######################"
        WRITE(67,*) "# SQRT(|Psi|^2)  |  SQRT( |F^(-1)[F(Psi)]|^2 ) "
        WRITE(67,'(*(ES16.8))') L2Norm(DCMPLX(Psi_0, 0d0), dx), L2Norm(FFTtest, dx)
        CLOSE(67)
        DEALLOCATE(FFTtest)
    ENDIF
    
    DEALLOCATE(EigVals)
    DEALLOCATE(Hamiltonian)
    
!-- COMPUTING THE TIME EVOLUTION ---------------------------------------

    WRITE(*,"(A)") "Starting to compute the time evolution of the ground state ..."
    
    ALLOCATE(P_vec(Pars%NdivX))
    CALL MomentumSpace(Pars, P_vec, dp)
    CALL Checkpoint(Pars%Debug, "Momentum space computed successfully.")

    Psi_t = DCMPLX(Psi_0, 0d0)
    CALL PrintTimeEvol(output, 0, Psi_t, X_vec, TT_str)
    
    WRITE(*,"(A)") "Starting the Split Operator Method computations ..."
    
    DO tstep=2, size(T_vec)
        !Update the potential array to the new time step
        CALL HarmonicOscillatorPot(Pars, X_vec, v_pot, T_vec(tstep) )
    
        ! 1) half potential operation
        Psi_t = Psi_t * CDEXP( DCMPLX(0d0, -v_pot*(dt*0.5d0)) )
        
        ! 2) Fourier transform of Psi_t
        CALL ComplexFFT(Pars, Psi_t, -1)
        
        ! 3) kinetic term operation in momentum space
        Psi_t = Psi_t * CDEXP( DCMPLX(0d0, -0.5d0*(P_vec**2)*dt) )
        
        ! 4) Inverse Fourier transform of Psi_t
        CALL ComplexFFT(Pars, Psi_t, +1)        
        Psi_t = Psi_t/DCMPLX(Pars%NdivX, 0d0) !normalization
        
        ! 5) second half potential operation
        Psi_t = Psi_t * CDEXP( DCMPLX(0d0, -v_pot*(dt*0.5d0)) )
        
        !Saving the computed time evolution in a file for each time step
        IF ( (MOD(tstep,100) .eq. 0) .or. (tstep .eq. size(T_vec)) ) THEN
            CALL PrintTimeEvol(output, tstep, Psi_t, X_vec, TT_str)
            IF (Pars%Debug) THEN
                CALL PrintPotEvol(output, tstep, v_pot, X_vec)
            ENDIF
        ENDIF
    ENDDO
    CALL Checkpoint(Pars%Debug, "Time evolution computation complete.")
    
!   Deallocating the memory
    DEALLOCATE(X_vec, v_pot, P_vec)
    DEALLOCATE(Psi_0, Psi_t)
    DEALLOCATE(T_vec)
    CALL Checkpoint(Pars%Debug, "Memory corretly deallocated. Execution complete!")
    
END PROGRAM    
