! --- Exercise 8 - Density Matrices ---------------------------------- !
! The aim of this program is to test the subroutines that are contained 
! in the module "DensityMatrices" on a system composed by 2 particles
! with Spin 1/2.
! It is part of Week 8 assignment.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !


PROGRAM DensityMat2SpinHalf
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
    f_config = "Config_2SpinHalf.txt"

    CALL InitParamsFromFile(f_config, Pars, argNames)
    
    CALL InitDebugger(Pars%Debug)
    CALL Checkpoint(Pars%Debug, "Input parameters set.")
    
    WRITE(*,"(A)") "Parameters configured to the following values: "
    WRITE(*,"(A,I2)")   "    Hdim = ", Pars%Hdim
    WRITE(*,"(A,I2)")   "    Nsys = ", Pars%Nsys
    WRITE(*,"(A,L)" )   "    Sep  = ", Pars%Sep
    WRITE(*,*)
    
    IF ( (Pars%Hdim.ne.2) .OR. (Pars%Nsys.ne.2) ) THEN
        WRITE(*,*) "This program require to have in input the parameters 'Hdim' and 'Nsys' both = 2."
        WRITE(*,*) "Please, edit the file: "//f_config//" accordingly and restart."
        STOP
    ENDIF
    
!-- TESTING THE DENSITY MATRIX SUBROUTINES ON A 2-SPIN 1/2 SYSTEM
    
!   Allocate memory to store the state
    CALL AllocPurePsi(Pars, psi)
    CALL Checkpoint(Pars%Debug, "Memory to store the state successfully allocated.")
    
!   Initialize it in a separable or non-separable state
    IF (Pars%Sep) THEN
        ! |Psi> = |0>_A (x) 1/sqrt(2)*(|0>_B + |1>_B)
        psi%WF = (/ DCMPLX(1d0,0d0), DCMPLX(0d0,0d0), DCMPLX(1d0/DSQRT(2d0),0d0), DCMPLX(1d0/DSQRT(2d0),0d0) /)
    ELSE
        ! |Psi> = 1/sqrt(2)*(|01> - |10>)
        psi%WF = (/ DCMPLX(0d0,0d0), DCMPLX(1d0/DSQRT(2d0),0d0), DCMPLX(-1d0/DSQRT(2d0),0d0), DCMPLX(0d0,0d0) /)
    ENDIF
    CALL Checkpoint(Pars%Debug, "State of the system set.")
    
    IF ( (Pars%Debug).AND.(.NOT.psi%separable) ) THEN
        WRITE(*,*) "--> Norm of psi: ", L2Norm(psi%WF)     !it is of complex type but Imaginary part should be 0
    ENDIF
    
!   Compute density matrix from pure psi
    CALL PureDensityMat(psi, DensityMat)
    CALL Checkpoint(Pars%Debug, "Density Matrix computed.")
    
!   Writing density matrix on file
    CALL WriteOnFile(DensityMat, "DensityMat_2SpinHalf.txt")
    CALL Checkpoint(Pars%Debug, "Density Matrix saved on file.")
    
!   Computing reduced density matrix of left system
    CALL RedDensityMat(Pars, DensityMat, ReducedDensityMat, .TRUE.)
    CALL Checkpoint(Pars%Debug, "Left density matrix computed.")
    
    CALL WriteOnFile(ReducedDensityMat, "RedMat_Left_2SpinHalf.txt")
    CALL Checkpoint(Pars%Debug, "Left density matrix saved on file.")
 
!   Computing reduced density matrix of right system
    CALL RedDensityMat(Pars, DensityMat, ReducedDensityMat, .FALSE.)
    CALL Checkpoint(Pars%Debug, "Right density matrix computed.")
    
    CALL WriteOnFile(ReducedDensityMat, "RedMat_Right_2SpinHalf.txt") 
    CALL Checkpoint(Pars%Debug, "Right density matrix saved on file.")
    
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
