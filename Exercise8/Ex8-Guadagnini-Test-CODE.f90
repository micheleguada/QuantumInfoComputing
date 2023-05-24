! --- Exercise 8 - Density Matrices ---------------------------------- !
! The aim of this program is to test the subroutines for creating
! density matrices from a pure state |psi> that are contained in the 
! module "DensityMatrices".
! It is part of Week 8 assignment.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !


PROGRAM DensityMatOfPureStates
    USE Debugger
    USE DensityMatrices
    IMPLICIT NONE
    
    TYPE(Parameters) Pars
    TYPE(PureState) psi
    
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: DensityMat
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ReducedDensityMat
    
    CHARACTER(len=40) :: f_config
    CHARACTER(len=8), DIMENSION(4) :: argNames
    
!-- READING THE ARGUMENTS FROM FILE --------------------------------------------------
    
!   Names of the input parameters (in the same order as defined in the TYPE(Parameters) definition)
    argNames = (/"Hdim  = ", "Nsys  = ", "Sep   = ", "Debug = "/)
    
!   Parameters file name
    f_config = "Config_Pars.txt"

    CALL InitParamsFromFile(f_config, Pars, argNames)
    
    CALL InitDebugger(Pars%Debug)
    CALL Checkpoint(Pars%Debug, "Input parameters set.")
    
    WRITE(*,"(A)") "Parameters configured to the following values: "
    WRITE(*,"(A,I2)")   "    Hdim = ", Pars%Hdim
    WRITE(*,"(A,I2)")   "    Nsys = ", Pars%Nsys
    WRITE(*,"(A,L)" )   "    Sep  = ", Pars%Sep
    WRITE(*,*)
    
!-- TESTING THE GENERAL DENSITY MATRIX SUBROUTINES
    
!   Allocate memory to store the state
    CALL AllocPurePsi(Pars, psi)
    CALL Checkpoint(Pars%Debug, "Memory to store the state successfully allocated.")
    
!   Initialize it at random and normalize
    CALL RandomPsi(psi)
    CALL Checkpoint(Pars%Debug, "State randomly intialized.")
    
    IF ( (Pars%Debug).AND.(.NOT.psi%separable) ) THEN
        WRITE(*,*) "--> Norm of psi: ", L2Norm(psi%WF)     !it is of complex type but Imaginary part should be 0
    ENDIF
    
!   Compute density matrix from pure psi
    CALL PureDensityMat(psi, DensityMat)
    CALL Checkpoint(Pars%Debug, "Pure density matrix computed.")
    
!   Writing density matrix on file
    CALL WriteOnFile(DensityMat, "DensityMat.txt")
    CALL Checkpoint(Pars%Debug, "Density matrix saved on file.")
    
!   FROM NOW ON THE PROGRAM IS EXECUTED ONLY IF # OF SUBSYSTEMS: N=2
    IF (Pars%Nsys > 2) THEN
        WRITE(*,*) "For # of subsystems greater than 2 this program ends here."
        STOP
    ENDIF
    
!   Computing reduced density matrix of left system
    CALL RedDensityMat(Pars, DensityMat, ReducedDensityMat, .TRUE.)
    CALL Checkpoint(Pars%Debug, "Left reduced density matrix computed.")
    
    CALL WriteOnFile(ReducedDensityMat, "RedMat_Left.txt")
    CALL Checkpoint(Pars%Debug, "Left reduced density matrix saved on file.")
 
!   Computing reduced density matrix of right system
    CALL RedDensityMat(Pars, DensityMat, ReducedDensityMat, .FALSE.)
    CALL Checkpoint(Pars%Debug, "Right reduced density matrix computed.")
    
    CALL WriteOnFile(ReducedDensityMat, "RedMat_Right.txt")
    CALL Checkpoint(Pars%Debug, "Right reduced density matrix saved on file.")
    
!   Memory deallocation
    IF (ALLOCATED(psi%WF)) THEN 
        DEALLOCATE(psi%WF)
    ENDIF
    IF (ALLOCATED(DensityMat)) THEN 
        DEALLOCATE(DensityMat)
    ENDIF
    IF (ALLOCATED(ReducedDensityMat)) THEN 
        DEALLOCATE(ReducedDensityMat)
    ENDIF
    CALL Checkpoint(Pars%Debug, "Memory deallocated. Execution complete.")
    
END PROGRAM
