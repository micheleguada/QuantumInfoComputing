! --- Exercise 9 - 1D Ising model ------------------------------------ !
! The aim of this program is to print the first k  eigenvalues and 
! eigenvectors by diagonalizing the hamiltonian:
!   H = lambda * SUM(P_z^i) + SUM(P_x^i P_x^(i+1))
! where P_z and P_x denotes the Pauli matrices.
! It is part of Week 9-10 assignment.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !


PROGRAM IsingModel1D
    USE Debugger
    USE SpinsHamiltonian
    IMPLICIT NONE
    
    TYPE(Parameters) Pars    
    CHARACTER(len=40) :: f_config
    CHARACTER(len=8), DIMENSION(4) :: argNames
    
    CHARACTER(len=40) :: vals_file, vecs_file
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EigVals
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ham
    
    vals_file = "EigVals"   ! + .dat"
    vecs_file = "EigVecs.dat"
    
!-- READING THE ARGUMENTS FROM FILE ------------------------------------
    
!   Names of the input parameters (in the same order as defined in the TYPE(Parameters) definition)
    argNames = (/"NN    = ", "kk    = ", "Lambda= ", "Debug = "/)
    
!   Parameters file name
    f_config = "Config_Pars.txt"

    CALL InitParamsFromFile(f_config, Pars, argNames)
    
    CALL Checkpoint(Pars%Debug, "Input parameters set.")
    CALL InitDebugger(Pars%Debug)
    
    WRITE(*,"(A)") "Parameters configured to the following values: "
    WRITE(*,"(A,I2)")    "    NN     = ", Pars%NN
    WRITE(*,"(A,I2)")    "    kk     = ", Pars%kk
    WRITE(*,"(A,F5.3)" ) "    Lambda = ", Pars%Lambda
    WRITE(*,*)
    
!-- Memory allocation
    IF (.not.(ALLOCATED(Ham))) THEN
        ALLOCATE(Ham(2**Pars%NN,2**Pars%NN))
    ENDIF
    IF (.not.(ALLOCATED(EigVals))) THEN
        ALLOCATE(EigVals(2**Pars%NN))
    ENDIF    
    CALL Checkpoint(Pars%Debug, "Memory corretly allocated.")
    
!-- COMPUTE THE MATRIX REPRESENTATION OF THE HAMILTONIAN ---------------

    WRITE(*,"(A)") "Computing the Hamiltonian..."

    Ham = Pars%Lambda * Ham_Zsum(Pars)
    CALL Checkpoint(Pars%Debug, "Z contribution computed.")
    
    Ham = Ham + Ham_Xsum(Pars)
    CALL Checkpoint(Pars%Debug, "X contribution computed.")
    
    ! check the computed hamiltonian
    IF (Pars%Debug) THEN
        CALL PrintEigVecs("Ham_Z_debug.dat", Ham_Zsum(Pars), 2**Pars%NN) !check the Z contribution
        CALL PrintEigVecs("Ham_X_debug.dat", Ham_Xsum(Pars), 2**Pars%NN) !check the X contribution
        CALL PrintEigVecs("Ham_debug.dat"  , Ham           , 2**Pars%NN)
    ENDIF
    
!-- DIAGONALIZE THE HAMILTONIAN ----------------------------------------

    WRITE(*,"(A)") "Diagonalizing the Hamiltonian..."

    CALL SymmetricEigenpairs(2**Pars%NN, Ham, EigVals)
    CALL Checkpoint(Pars%Debug, "Hamiltonian diagonalization completed.")
    
!-- PRINT THE FIRST k EIGENPAIRS ---------------------------------------

    WRITE(*,"(A)") "Printing the results..."

    CALL PrintEigVals(Pars, vals_file, EigVals)
    CALL Checkpoint(Pars%Debug, "First k eigenvalues stored on file: "//vals_file)
    
    IF (Pars%Debug) THEN 
        ! check diagonalized hamiltonian (eigenvectors)
        CALL PrintEigVecs(vecs_file, Ham, Pars%kk)
    ENDIF
    
!-- Memory deallocation
    IF (ALLOCATED(Ham)) THEN
        DEALLOCATE(Ham)
    ENDIF
    IF (ALLOCATED(EigVals)) THEN
        DEALLOCATE(EigVals)
    ENDIF       
    CALL Checkpoint(Pars%Debug, "Memory corretly deallocated. Execution complete.")   
    
END PROGRAM
