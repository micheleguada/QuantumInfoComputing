! This module diagonalize a matrix, compute the 
! normalized spacings between an array of eigenvalues and
! create and normalize an historgram. It is part of the
! Week 5 assignment.
!
! AUTHOR: Michele Guadagnini - ID 1230663

MODULE NormalizedSpacings
    USE Debugger
    IMPLICIT NONE
    CONTAINS
    
    !Hermitian initialization of the matrix
    SUBROUTINE RandInitHermitianMat(NN, AA)
        INTEGER NN, ii, jj
        DOUBLE COMPLEX, DIMENSION(:,:) :: AA
        DOUBLE PRECISION Re, Im
        
        CALL RANDOM_SEED()
        DO ii=1,NN
            DO jj=ii,NN
                CALL RANDOM_NUMBER(Re)
                Re = (Re*2-1)*10
                IF (ii .eq. jj) THEN
                    AA(jj,ii) = DCMPLX(Re,0.)
                ELSE
                    CALL RANDOM_NUMBER(Im)
                    Im = (Im*2-1)*10
                    AA(ii,jj) = DCMPLX(Re,Im)
                ENDIF
            ENDDO
        ENDDO     
                
        RETURN
    END SUBROUTINE
    
    !It computes eigenvalues
    SUBROUTINE HermEigenvalues(NN, AA, Eigs)
        INTEGER NN, INFO, LWORK
        DOUBLE COMPLEX, DIMENSION(:,:) :: AA
        DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: WORK
        DOUBLE PRECISION, DIMENSION(:) :: Eigs
        DOUBLE PRECISION, DIMENSION(3*NN-2) :: RWORK
        
        !Query the optimal workspace
        ALLOCATE(WORK(2))
        LWORK = -1
        CALL ZHEEV("N","U", NN, AA, NN, Eigs, WORK, LWORK, RWORK, INFO)
        LWORK = INT(WORK(1), kind=4)
        DEALLOCATE(WORK)
        
        IF (INFO /= 0) THEN
            WRITE(*,"(A,I2)") " Unable to set optimal workspace. INFO: ", INFO
            STOP
        ENDIF
        
        !Solving the eigenproblem
        ALLOCATE(WORK(LWORK))
        CALL ZHEEV("N","U", NN, AA, NN, Eigs, WORK, LWORK, RWORK, INFO)
        WRITE(*,"(A,I2)") " Eigenproblem solved with exit status: ",INFO
        
        DEALLOCATE(WORK)
        
        RETURN
    END SUBROUTINE  
    
    !It computes the normalized spacings
    SUBROUTINE NormSpacings(Eigs, ss)
        DOUBLE PRECISION, DIMENSION(:) :: Eigs
        DOUBLE PRECISION, DIMENSION(size(Eigs)-1) :: ss
        DOUBLE PRECISION S_mean
        INTEGER ii
        
        DO ii = 1, size(Eigs)-1
            ss(ii) = Eigs(ii+1) - Eigs(ii)
        ENDDO
        S_mean = SUM(ss)/(MAX(1,size(Eigs)-1))
        ss = ss / S_mean
        
        RETURN
    END SUBROUTINE
    
    !It computes the histogram and normalizes it to unit area
    SUBROUTINE ComputePDF(ss, Nbin, probs, xmids, xmin, xmax)
        DOUBLE PRECISION, DIMENSION(:) :: ss
        INTEGER Nbin, ii, jj, kk
        DOUBLE PRECISION xmin, xmax, deltax, Integral
        DOUBLE PRECISION, DIMENSION(Nbin) :: xmids, probs, counts
        
        Nbin = MAX(1,Nbin)  !to avoid the stop of the program if Nbin=0
        deltax = ((xmax-xmin)/Nbin)
        
        DO kk = 1,Nbin
            xmids(kk) = deltax*(kk+0.5)
        ENDDO
        
        DO ii = 1,size(ss)
            DO jj = 1,Nbin           
                IF ( (ss(ii) .ge. (xmin+deltax*jj)) .and. (ss(ii) .lt. (xmin+deltax*(jj+1))) ) THEN
                    counts(jj) = counts(jj)+1
                    EXIT
                ENDIF
            ENDDO
        ENDDO
        
        Integral = SUM(counts*deltax)
        probs = counts / Integral        
    
        RETURN
    END SUBROUTINE
    
    !It prints the results in columns
    SUBROUTINE PrintColumnsOnFile(filename, xmids, probs)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:) :: xmids, probs
        INTEGER ii
        
        OPEN(unit=10 , file=filename, status='unknown')
        DO ii = 1,size(probs)
            WRITE(10,"(ES16.8,ES16.8)") xmids(ii), probs(ii)
        ENDDO
        CLOSE(10)
        
        RETURN
    END SUBROUTINE
    
    !It prints an array in column
    SUBROUTINE PrintColumnOnFile(filename, ss)
        CHARACTER(len=*) filename
        DOUBLE PRECISION, DIMENSION(:) :: ss
        INTEGER ii
        
        OPEN(unit=10 , file=filename, status='unknown')
        DO ii = 1,size(ss)
            WRITE(10,"(ES16.8)") ss(ii)
        ENDDO
        CLOSE(10)
        
        RETURN
    END SUBROUTINE
    
END MODULE
