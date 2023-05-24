! --- Exercise 10 - Renormalization Group with 1D Ising model -------- !
! The aim of this program is to apply Renormalization Group methods to
! the 1D Ising model with hamiltonian:
!   H = lambda * SUM(P_z^i) + SUM(P_x^i P_x^(i+1))
! where P_z and P_x denotes the Pauli matrices.
! It is part of Week 10-11 assignment.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !


PROGRAM RGIsingModel1D
    USE Debugger
    USE RSRG
    IMPLICIT NONE
    
    TYPE(Parameters) Pars    
    CHARACTER(len=40), PARAMETER :: f_config  = "Config_Pars.txt"     !Parameters file name
    CHARACTER(len=40), PARAMETER :: vals_file = "EigVals"
    CHARACTER(len=8), DIMENSION(7) :: argNames
    
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Init_EigVals, EigVals
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Init_Ham, Init_EigVectors, Init_HamL, Init_HamR
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ham, EigVectors, HamL, HamR
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TruncHam, TruncHamL, TruncHamR
    
    INTEGER idx, tempSz
    DOUBLE PRECISION :: OldGS, OldSz, DescrSz
    
    OldGS = 0d0
    OldSz = 0d0
    
!-- READING THE ARGUMENTS FROM FILE ------------------------------------
    
!   Names of the input parameters (in the same order as defined in the TYPE(Parameters) definition)
    argNames = (/"NN    = ", "TruncN= ", "Lambda= ", "ExitTh= ", "MaxItr= ", "kk    = ", "Debug = "/)

    CALL InitParamsFromFile(f_config, Pars, argNames)
    CALL Checkpoint(Pars%Debug, "Input parameters set.")
    
    CALL CatchError((Pars%TruncN .gt. Pars%NN), "The size of the truncated hamiltonian must "&
                    //"be smaller or equal than the initial hamiltonian size.", Fatal=.TRUE.)
    
!-- Memory allocation
    IF (.not.(ALLOCATED(Ham))) THEN
        ALLOCATE(Init_Ham(2**Pars%NN    , 2**Pars%NN    ))
        ALLOCATE( Ham(2**(2*Pars%TruncN), 2**(2*Pars%TruncN)))
    ENDIF
    IF (.not.(ALLOCATED(EigVectors))) THEN
        ALLOCATE(Init_EigVectors(2**Pars%NN    , 2**Pars%NN    ))
        ALLOCATE( EigVectors(2**(2*Pars%TruncN), 2**(2*Pars%TruncN)))
    ENDIF
    IF (.not.(ALLOCATED(EigVals))) THEN
        ALLOCATE(Init_EigVals(2**Pars%NN        ))
        ALLOCATE(     EigVals(2**(2*Pars%TruncN)))
    ENDIF
    IF (.not.(ALLOCATED(HamL)) .and. .not.ALLOCATED(HamR)) THEN
        ALLOCATE(Init_HamL(2**Pars%NN    , 2**Pars%NN    ))
        ALLOCATE(Init_HamR(2**Pars%NN    , 2**Pars%NN    ))
        ALLOCATE( HamL(2**(2*Pars%TruncN), 2**(2*Pars%TruncN)))
        ALLOCATE( HamR(2**(2*Pars%TruncN), 2**(2*Pars%TruncN)))
    ENDIF
    IF (.not.(ALLOCATED(TruncHam))) THEN
        ALLOCATE( TruncHam(2**Pars%TruncN, 2**Pars%TruncN))
        ALLOCATE(TruncHamL(2**Pars%TruncN, 2**Pars%TruncN))
        ALLOCATE(TruncHamR(2**Pars%TruncN, 2**Pars%TruncN))
    ENDIF
    CALL Checkpoint(Pars%Debug, "Memory corretly allocated.")
    
!-- COMPUTE THE INITIAL HAMILTONIAN MATRIX -----------------------------

    WRITE(*,"(A)") "Computing the initial Hamiltonian..."

    Init_Ham = Pars%Lambda * Ham_Zsum(Pars)
    CALL Checkpoint(Pars%Debug, "Z contribution computed.")
    
    Init_Ham = Init_Ham + Ham_Xsum(Pars)
    CALL Checkpoint(Pars%Debug, "X contribution computed.")
    
    ! prints to check the computed hamiltonian
    IF (Pars%Debug) THEN
        CALL PrintEigVecs("Ham_Z_debug.dat", Ham_Zsum(Pars)) !check the Z contribution
        CALL PrintEigVecs("Ham_X_debug.dat", Ham_Xsum(Pars)) !check the X contribution
        CALL PrintEigVecs("Ham_debug.dat"  , Init_Ham      )
    ENDIF
    
    WRITE(*,"(A)") "Diagonalizing the initial Hamiltonian..."
    Init_EigVectors = Init_Ham    !creating a copy in order to not overwrite it during diagonalization
    CALL SymmetricEigenpairs(2**Pars%NN, Init_EigVectors, Init_EigVals)
    CALL Checkpoint(Pars%Debug, "Hamiltonian diagonalization completed.")

    ! compute the initial left and right interaction terms    
    TruncHamL = DIdentity(2**(Pars%TruncN-1)) .kron. DREAL(PauliMat(1))
    TruncHamR = DREAL(PauliMat(1)) .kron. DIdentity(2**(Pars%TruncN-1))
    
!-- LOOP OF Real Space Renorm. Group ALGORITHM -------------------------
    WRITE(*,*) " "  !simple newline 

    ! number of particles in the system (to be updated at each step)
    DescrSz = DFLOAT(Pars%NN)
    
    DO idx = 1,Pars%MaxItr
        WRITE(*,"(A,I4)") "STARTING ITERATION #: ", idx
        
        !## 1) Projection ##!
        WRITE(*,"(A)") "  Projecting the Hamiltonian and interaction operators..."        
        IF (idx .eq. 1) THEN  !to have (NN) and iterated size (TruncN) independent
            ! padding the interaction terms
            tempSz = Pars%NN - Pars%TruncN      
            Init_HamL = TruncHamL .kron. DIdentity(2**tempSz)
            Init_HamR = TruncHamR .kron. DIdentity(2**tempSz)
            ! projections
            TruncHam  = Projection( Init_EigVectors(:,1:2**Pars%TruncN), Init_Ham )
            TruncHamL = Projection( Init_EigVectors(:,1:2**Pars%TruncN), Init_HamL)
            TruncHamR = Projection( Init_EigVectors(:,1:2**Pars%TruncN), Init_HamR)
        ELSE
            ! padding the interaction terms
            HamL = TruncHamL .kron. DIdentity(2**Pars%TruncN)
            HamR = TruncHamR .kron. DIdentity(2**Pars%TruncN)
            ! projections
            TruncHam  = Projection( EigVectors(:,1:2**Pars%TruncN), Ham )
            TruncHamL = Projection( EigVectors(:,1:2**Pars%TruncN), HamL)
            TruncHamR = Projection( EigVectors(:,1:2**Pars%TruncN), HamR)
        ENDIF
        CALL Checkpoint(Pars%Debug, "Projections completed.")
        
        !## 2) Double the system ##!
        WRITE(*,"(A)") "  Building the doubled Hamiltonian..."
        Ham = (TruncHam .kron. DIdentity(2**Pars%TruncN)) + & 
              (DIdentity(2**Pars%TruncN) .kron. TruncHam) + &
              (TruncHamL .kron. TruncHamR)
        CALL Checkpoint(Pars%Debug, "Doubled Hamiltonian computation completed.")     
        
        !## 3) Diagonalization ##!
        WRITE(*,"(A)") "  Diagonalizing the doubled Hamiltonian..."
        EigVectors = Ham  !creating a copy in order to not overwrite it
        CALL SymmetricEigenpairs(2**(2*Pars%TruncN), EigVectors, EigVals)
        CALL Checkpoint(Pars%Debug, "Hamiltonian diagonalization completed.")
        
        OldSz   = DescrSz     ! size of the described system in the previous step
        DescrSz = 2d0*DescrSz ! size of the present described system

        ! Check exit condition
        IF ( AbsDiff( OldGS, EigVals(1)/(DescrSz-1d0) ) .lt. Pars%ExitTh) THEN
            WRITE(*,"(A,ES16.8)") "  Improvement in the ground state eigenvalue: < ",Pars%ExitTh 
            WRITE(*,"(A,I4)")     "  Exiting at iteration: ", idx
            EXIT    !exit the RSRG loop
        ELSEIF (idx == Pars%MaxItr) THEN
            WRITE(*,"(A,I4)")     "  Reached maximum number of iterations: ",idx
            EXIT
        ELSE
            OldGS = EigVals(1)/(DescrSz-1d0)
        ENDIF
        WRITE(*,*) " "  !simple newline 
    ENDDO
    
!-- SAVING THE RESULTS -------------------------------------------------
    
    WRITE(*,*) " "  !simple newline
    WRITE(*,"(A,ES16.8)") "The final size of the described system is: ", DescrSz
    WRITE(*,"(A)")        "Saving the resulting eigenvalues..."
    CALL PrintEigVals(Pars, vals_file, EigVals, DescrSz, idx)
    CALL Checkpoint(Pars%Debug, "Resulting eigenvalues stored.")
    
!-- Memory deallocation
    IF (ALLOCATED(Ham)) THEN
        DEALLOCATE(Ham)
        DEALLOCATE(Init_Ham)
    ENDIF
    IF (ALLOCATED(EigVals)) THEN
        DEALLOCATE(EigVals)
        DEALLOCATE(Init_EigVals)
    ENDIF       
    IF (ALLOCATED(EigVectors)) THEN
        DEALLOCATE(EigVectors)
        DEALLOCATE(Init_EigVectors)
    ENDIF
    IF (ALLOCATED(HamL) .and. ALLOCATED(HamR)) THEN
        DEALLOCATE(HamL)
        DEALLOCATE(HamR)
        DEALLOCATE(Init_HamL)
        DEALLOCATE(Init_HamR)
    ENDIF
    IF (ALLOCATED(TruncHam)) THEN
        DEALLOCATE(TruncHam)
        DEALLOCATE(TruncHamL)
        DEALLOCATE(TruncHamR)
    ENDIF
    CALL Checkpoint(Pars%Debug, "Memory corretly deallocated. Execution complete.")   
    
END PROGRAM
