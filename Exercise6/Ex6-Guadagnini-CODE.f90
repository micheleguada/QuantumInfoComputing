! --- Exercise 6 - Continuous time-ind. Schroedinger Equation -------- !
! The aim of this program is to compute the first k eigenvalues and 
! eigenvectors of the continuous time-independent Schroedinger Equation
! for the harmonic oscillator with the Hamiltonian: H = p^2 + w^2 q^2.
! The chosen method is the Finite Difference Method.
! h_bar and 2m are set to be 1.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !

PROGRAM SchroedingerEq
    USE Debugger
    USE ContinuousTimeIndSE
    IMPLICIT NONE
    
    DOUBLE PRECISION ww, LL, hh             !frequency (omega), Length (L), unit step (h)
    DOUBLE PRECISION width                  !dimension of the system in unit of: 1/ww
    INTEGER kk, Ndiv                        !# of eigenpairs to print (k), # of division of the space (Ndiv)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: v_bas, v_pot     !basis vector, potential vector
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Hamiltonian
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EigVals          !eigenvalues
    CHARACTER(len=40) :: EigValsFile, EigVecsFile
    
    LOGICAL Debug
    INTEGER totargs, idx, ios, ii, expargs
    CHARACTER(len=10) argDebug
    CHARACTER(len=10), DIMENSION(4) :: args
    expargs = size(args)+1
    ios = 0
    Debug = .FALSE. 
    
!   Results file names
    EigValsFile = "Eigenvalues.dat"
    EigVecsFile = "Eigenfunctions.dat" 
    
!-- READING THE ARGUMENTS FROM COMMAND LINE --------------------------------------------------

    totargs = command_argument_count()
    IF ( (totargs .gt. expargs) .or. (totargs .lt. (expargs-1)) ) THEN
        CALL CatchError("Wrong number of arguments. Using defaults.", Fatal=.FALSE.)
        CALL InitDefaults(ww, width, Ndiv, kk)
    ENDIF
    
    IF (totargs == expargs) THEN
        CALL get_command_argument(expargs, argDebug)
        argDebug = TRIM(argDebug)
        IF (argDebug == "TRUE") THEN
            Debug = .TRUE.
            !Debugger informations printing
            CALL InitDebugger(Debug)
        ENDIF
        totargs = expargs-1
    ENDIF
        
    IF (totargs == expargs-1) THEN
        DO idx = 1,expargs-1
            CALL get_command_argument(idx, args(idx))
        ENDDO        
        READ(args,*,IOSTAT=ios) ww, width, Ndiv, kk
    ENDIF
    IF (ios .ne. 0) THEN 
        CALL CatchError("Invalid value for one or more parameters. Using Defaults", Fatal=.FALSE.)
        CALL InitDefaults(ww, width, Ndiv, kk)
    ENDIF
    CALL Checkpoint(Debug, "Parameters values set.")
    
!   Memory allocation
    ALLOCATE(v_bas(Ndiv))
    ALLOCATE(v_pot(Ndiv))
    ALLOCATE(Hamiltonian(Ndiv,Ndiv))
    ALLOCATE(EigVals(Ndiv))
    CALL Checkpoint(Debug, "Memory correctly allocated.")
    
!-- SOLVING THE SCHROEDINGER EQUATION --------------------------------------------------------

    WRITE(*,"(A)") "Starting to solve the Harmonic Oscillator Schroedinger equation..."
    WRITE(*,"(A)") "Values of the parameters are: "
    WRITE(*,"(A,F8.5)") "    w     = ", ww
    WRITE(*,"(A,F8.5)") "    width = ", width
    WRITE(*,"(A,I8)")   "    Ndiv  = ", Ndiv
    WRITE(*,"(A,I8)")   "    k     = ", kk
    WRITE(*,*)

!   Discretization
    CALL Discretize(ww, Ndiv, LL, hh, v_bas, width)
    CALL Checkpoint(Debug, "Discretization of the space successful.")
    
!   Compute potential
    CALL HarmonicOscillatorPot(ww, v_bas, v_pot)
    CALL Checkpoint(Debug, "Hamiltonian Potential computed successfully.")
    
    !it prints the discretization basis and potential if debug is active
    IF (Debug) THEN
        OPEN(unit=80 , file="test.txt", status='unknown')
        WRITE(80, "(A)") "# basis           potential"
        DO ii = 1,Ndiv
            WRITE(80,"(2ES16.8)") v_bas(ii), v_pot(ii)
        ENDDO
        CLOSE(80)
    ENDIF
    
!   Compute the hamiltonian matrix
    CALL InitHamiltonian(Ndiv, v_pot, Hamiltonian, hh)
    CALL Checkpoint(Debug, "Hamiltonian matrix computed successfully.")
    
    !it prints the hamiltonian matrix if debug is active and the matrix is small
    IF ((Debug) .and. (Ndiv .le. 10)) THEN
        OPEN(unit=80 , file="testMat.txt", status='unknown')
        WRITE(80, "(A)") "# hamiltonian matrix"
        DO ii = 1,Ndiv
            WRITE(80,"(*(F10.6))") Hamiltonian(ii,:)
        ENDDO
        CLOSE(80)
    ENDIF
    
!   Diagonalize the matrix
    CALL SymmetricEigenpairs(Ndiv, Hamiltonian, EigVals)
    CALL Checkpoint(Debug, "Matrix diagonalization completed.")
    
!   Storing the results on file
    WRITE(*,"(A)") "Printing the results on file..."
    CALL PrintEigVals(EigValsFile, EigVals, kk)
    CALL PrintEigVecs(EigVecsFile, v_bas, Hamiltonian, kk, hh)
    CALL Checkpoint(Debug, "Resulting eigenpairs printed on corresponding files.")
    
!   Deallocating the memory
    DEALLOCATE(v_bas)
    DEALLOCATE(v_pot)
    DEALLOCATE(Hamiltonian)
    DEALLOCATE(EigVals)
    CALL Checkpoint(Debug, "Memory corretly deallocated. Execution complete!")
    
END PROGRAM    
