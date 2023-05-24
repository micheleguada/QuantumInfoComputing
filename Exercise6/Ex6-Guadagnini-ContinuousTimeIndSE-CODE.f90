! This module discretizes the space, computes potential of the hamiltonian, 
! diagonalize a real symmetric matrix, compute the eigenvalues and eigenvectors. 
! It is part of the Week 6 assignment regarding continuous time-independent
! Schroedinger Equation.
!
! AUTHOR: Michele Guadagnini - ID 1230663


MODULE ContinuousTimeIndSE
    USE Debugger
    IMPLICIT NONE
    CONTAINS
    
    !Default parameters initialization
    SUBROUTINE InitDefaults(ww, width, Ndiv, kk)
        DOUBLE PRECISION ww, width
        INTEGER Ndiv, kk
        
        ww = 1d0
        width = 10d0
        Ndiv = 101
        kk = 4
        
        RETURN
    END SUBROUTINE
    
    !Discretization of the space
    SUBROUTINE Discretize(ww, Ndiv, LL, hh, v_bas, width) 
        DOUBLE PRECISION ww, hh, xmin, LL, width
        INTEGER Ndiv, ii
        DOUBLE PRECISION, DIMENSION(:) :: v_bas
        
        WRITE(*,"(A)") "Discretizing the space..."   
        
        LL = width/ww           !width of the system: L ~ width * 1/w
        xmin = -LL/2d0
        hh = LL/DFLOAT(Ndiv-1)  !Ndiv-1 in order to have symmetric values of x with respect to 0
        DO ii=1,Ndiv
            v_bas(ii) = xmin + hh*DFLOAT(ii-1)
        ENDDO
        
        RETURN
    END SUBROUTINE
    
    !It computes the harmonic oscillator potential from the basis of discretization
    SUBROUTINE HarmonicOscillatorPot(ww, v_bas, v_pot)
        DOUBLE PRECISION ww
        DOUBLE PRECISION, DIMENSION(:) :: v_bas, v_pot
    
        WRITE(*,"(A)") "Computing the harmonic oscillator potential..."
    
        v_pot = ww*ww*(v_bas*v_bas)
        
        RETURN
    END SUBROUTINE
    
    !Initialization of the matrix
    SUBROUTINE InitHamiltonian(Ndiv, v_pot, Ham, hh)
        INTEGER Ndiv, ii
        DOUBLE PRECISION hh
        DOUBLE PRECISION, DIMENSION(:,:) :: Ham
        DOUBLE PRECISION, DIMENSION(:) :: v_pot
        
        WRITE(*,"(A)") "Initializing the hamiltonian matrix..."
        
        Ham = 0d0
        DO ii=1,Ndiv
            !main diagonal
            Ham(ii,ii) = 2d0/(hh**2) + v_pot(ii)
            
            !upper diagonal
            IF (ii .ne. Ndiv) THEN
                Ham(ii,ii+1) = - (1d0)/(hh**2)
            ENDIF
            
            !lower diagonal is not necessary due to LAPACK's subroutine implementation
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
    
    !It prints the first k eigenfunctions in parallel columns
    SUBROUTINE PrintEigVecs(filename, v_bas, eigvs, kk, hh)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:,:) :: eigvs
        DOUBLE PRECISION, DIMENSION(:) :: v_bas
        INTEGER ii, kk, jj
        DOUBLE PRECISION hh
        
        OPEN(unit=10 , file=filename, status='unknown')
        WRITE(10,"(A)") "# basis           -> k-esim eigenfunctions"
        DO ii = 1,size(eigvs, dim=1)
            !eigenfunctions are normalized to unit norm by dividing it by 1/sqrt(hh)
            WRITE(10,'(*(ES16.8))') v_bas(ii), ( eigvs(ii,jj)/DSQRT(hh), jj=1,kk )
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
            WRITE(20,"(I6, F15.8)") (2*ii-1), eigs(ii)
        ENDDO
        CLOSE(20)
        
        RETURN
    END SUBROUTINE
    
END MODULE
