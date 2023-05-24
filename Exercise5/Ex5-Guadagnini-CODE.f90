! --- Exercise 5 - Eigenproblem & Random Matrix Theory --------------- !
! The aim of this program is to study the distribution of normalized 
! spacings of eigenvalues of different type of random matrices. The 2 
! types considered are: HERMITIAN matrix and DIAGONAL matrix.
!
! AUTHOR: Michele Guadagnini - ID 1230663
! -------------------------------------------------------------------- !


PROGRAM RandMatTheory
    USE Debugger
    USE NormalizedSpacings
    IMPLICIT NONE
    
    INTEGER :: NN = 200                               !default size of the matrix
    INTEGER :: Nbin = 20                              !default number of bins
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: AA    !Matrix
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ss, ss_D    !vector of spacings
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Eigs, Eigs_D  !vector of eigenvalues
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: probs, xmids, probs_D, xmids_D
    DOUBLE PRECISION xmin, xmax
    
    CHARACTER(10) debugstr, Nstr, Nbinstr
    LOGICAL :: Debug = .FALSE.
    INTEGER totargs, ios, INFO
    
!-- READING THE ARGUMENTS FROM COMMAND LINE --------------------------------------------------

!   (Activate debug) Reading the size
    totargs = command_argument_count()
    IF ( (totargs .gt. 3) .or. (totargs == 1) ) THEN
        CALL CatchError("Wrong number of arguments. Using hard-coded defaults.", Fatal=.FALSE.)
        totargs = 0
    ENDIF
    
    IF (totargs == 3) THEN
        CALL get_command_argument(1, Nstr)
        CALL get_command_argument(2, Nbinstr)
        CALL get_command_argument(3, debugstr)
        debugstr = TRIM(debugstr)
        IF (debugstr == "TRUE") THEN
            Debug = .TRUE.
            CALL InitDebugger(Debug)
        ENDIF
        READ(Nstr,*,IOSTAT=ios) NN
        READ(Nbinstr,*,IOSTAT=ios) Nbin
    ELSEIF (totargs == 2) THEN
        CALL get_command_argument(1, Nstr)
        CALL get_command_argument(2, Nbinstr)
        READ(Nstr,*,IOSTAT=ios) NN
        READ(Nbinstr,*,IOSTAT=ios) Nbin
    ENDIF
    IF (ios .ne. 0) THEN 
        CALL CatchError("Invalid value for the size of the matrix. Using Defaults", Fatal=.FALSE.)
    ENDIF
    CALL Checkpoint(Debug, "Matrix size set successfully.")
    
!-- DIAGONAL MATRIX SPACINGS COMPUTATION ---------------------------------------------------
    WRITE(*,"(A)") "Starting Diagonal matrix spacings computation."

!   Allocate the memory
    ALLOCATE(Eigs_D(NN))
    ALLOCATE(ss_D(NN-1))
    ALLOCATE(probs_D(Nbin))
    ALLOCATE(xmids_D(Nbin))

!   Random initialization of the diagonal matrix
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(Eigs_D)
    Eigs_D = (Eigs_D*2-1)*10
    CALL DLASRT("I", NN, Eigs_D, INFO)       !sorting the diagonal elements in increasing order
    CALL Checkpoint(Debug, "Random Eigenvalues successfully generated.")

!   Computing the normalized spacings
    CALL NormSpacings(Eigs_D, ss_D)
    CALL Checkpoint(Debug, "Normalized spacings successfully computed.")
    
!   Computing the PDF
    xmin = MIN(0.0, MINVAL(ss_D))
    xmax = MIN(5.0, MAXVAL(ss_D))
    CALL ComputePDF(ss_D, Nbin, probs_D, xmids_D, xmin, xmax)
    CALL Checkpoint(Debug, "Diag. PDF successfully computed.")
    
!   Print on file the Diagonal matrix spacings
    CALL PrintColumnsOnFile("DiagSpacings.dat", xmids_D, probs_D)
    CALL Checkpoint(Debug, "Results printed on file.")    
    WRITE(*,"(A)") "Diagonal matrix spacings distribution computed and saved."
    
!   Deallocate memory
    DEALLOCATE(Eigs_D)
    DEALLOCATE(ss_D)
    DEALLOCATE(probs_D)
    DEALLOCATE(xmids_D)
    
!-- HERMITIAN MATRIX SPACINGS COMPUTATION --------------------------------------------------
    WRITE(*,"(A)") "Starting Hermitian matrix spacings computation."
    
!   Allocate the memory for the matrix
    ALLOCATE(AA(NN,NN))
    ALLOCATE(Eigs(NN))
    ALLOCATE(ss(NN-1))
    CALL Checkpoint(Debug, "Memory correctly allocated.")

!   Random initialization of the matrix
    CALL RandInitHermitianMat(NN, AA)
    CALL Checkpoint(Debug, "Matrix initialization complete.")
    
!   Computing the eigenvalues
    CALL HermEigenvalues(NN, AA, Eigs)
    CALL Checkpoint(Debug, "Eigenvalues successfully computed.")

!   Computing the normalized spacings
    CALL NormSpacings(Eigs, ss)
    CALL Checkpoint(Debug, "Normalized spacings successfully computed.")

!   Computing the PDF
    ALLOCATE(probs(Nbin))
    ALLOCATE(xmids(Nbin))
    xmin = MIN(0.0, MINVAL(ss))
    xmax = MIN(4.0, MAXVAL(ss))
    CALL ComputePDF(ss, Nbin, probs, xmids, xmin, xmax)
    CALL Checkpoint(Debug, "Herm. PDF successfully computed.")
    
!   Print on file the Hermitian matrix spacings
    CALL PrintColumnsOnFile("HermSpacings.dat", xmids, probs)
    CALL Checkpoint(Debug, "Results printed on file.")    
    WRITE(*,"(A)") "Hermitian matrix spacings distribution computed and saved."   
    
!   Deallocate memory and stop program
    DEALLOCATE(AA)
    DEALLOCATE(Eigs)
    DEALLOCATE(ss)
    DEALLOCATE(probs)
    DEALLOCATE(xmids)
    CALL Checkpoint(Debug, "Execution is complete!")
    STOP
END PROGRAM
