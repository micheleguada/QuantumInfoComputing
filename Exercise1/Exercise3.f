      MODULE LoopMul
c     Only square matrices by now with this program
c     Msize is the parameter that controls the shape of the matrices
         INTEGER*2, PARAMETER :: Msize = 1000
         CONTAINS
c     Matrix multiplication with 1st loop order 
         FUNCTION MatLoopMul1(AA, BB) RESULT(CC)
            INTEGER*2 ii, jj, kk, somma
            INTEGER*2, DIMENSION(Msize,Msize) :: AA,BB,CC
            DO ii = 1,Msize
               DO jj = 1,Msize
                  somma = 0
                  DO kk = 1,Msize
                     somma = somma + AA(ii,kk)*BB(kk,jj)
                  ENDDO
                  CC(ii,jj) = somma
               ENDDO
            ENDDO            
            RETURN
         END FUNCTION
         
c     Matrix multiplication with 2nd loop order 
         FUNCTION MatLoopMul2(AA, BB) RESULT(CC)
            INTEGER*2 ii, jj, kk, somma
            INTEGER*2, DIMENSION(Msize,Msize) :: AA,BB,CC
            DO jj = 1,Msize
               DO ii = 1,Msize
                  somma = 0
                  DO kk = 1,Msize
                     somma = somma + AA(ii,kk)*BB(kk,jj)
                  ENDDO
                  CC(ii,jj) = somma
               ENDDO
            ENDDO            
            RETURN
         END FUNCTION
         
c     Subroutine for printing the matrices
         SUBROUTINE PrMat(DD)
            INTEGER*2 :: ll, hh
            INTEGER*2, DIMENSION(Msize,Msize) :: DD
            DO ll = 1,Msize
               DO hh = 1,Msize
                  WRITE(*,fmt='(I8)', advance="no") DD(ll,hh)
               ENDDO
               WRITE(*,*)
            ENDDO
            RETURN
         END SUBROUTINE
      
      END MODULE
      
      PROGRAM MatrixTest
         USE LoopMul
         IMPLICIT NONE
         REAL*4, DIMENSION(Msize,Msize) :: MM
         INTEGER*2, DIMENSION(Msize,Msize) :: AA
         INTEGER*2, DIMENSION(Msize,Msize) :: BB
         INTEGER*2, DIMENSION(Msize,Msize) :: CC
         REAL*4 start, fin
c     Initialization of the matrices
         CALL RANDOM_SEED()
         CALL RANDOM_NUMBER(MM)
         AA = FLOOR(10*MM)
         CALL RANDOM_NUMBER(MM)
         BB = FLOOR(10*MM)
c     Printing the matrices
         !PRINT *, 'Matrix A:'
         !CALL PrMat(AA)
         !PRINT *, 'Matrix B:'
         !CALL PrMat(BB)
c     First Loop order
         CALL CPU_TIME(start)
         CC = MatLoopMul1(AA,BB)
         CALL CPU_TIME(fin)
         PRINT *, 'Mat. C with 1st loop order:'
         !CALL PrMat(CC)
         PRINT *, 'Exec. time [s]: ',fin-start
c     Second Loop order
         CALL CPU_TIME(start)
         CC = MatLoopMul2(AA,BB)
         CALL CPU_TIME(fin)
         PRINT *, 'Mat. C with 2nd loop order:'
         !CALL PrMat(CC)
         PRINT *, 'Exec. time [s]: ',fin-start
c     Fortran intrinsic function
         CALL CPU_TIME(start)
         CC = MATMUL(AA,BB)
         CALL CPU_TIME(fin)
         PRINT *, 'Mat. C with MATMUL func:'
         !CALL PrMat(CC)
         PRINT *, 'Exec. time [s]: ',fin-start
         
         STOP
      END PROGRAM MatrixTest
         
