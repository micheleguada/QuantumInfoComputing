MODULE ComplexMatrix
   IMPLICIT NONE
   TYPE CMatrix
      INTEGER, DIMENSION(2) :: N
      DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: Elem, ElemAdj
      DOUBLE COMPLEX CTrace
      DOUBLE COMPLEX CDet
   END TYPE
   
   INTERFACE OPERATOR (.Tr.)
      MODULE PROCEDURE Trace
   END INTERFACE
   
   INTERFACE OPERATOR (.Adj.)
      MODULE PROCEDURE Adjoint
   END INTERFACE
   
   CONTAINS
   
   FUNCTION Trace(CMat) RESULT(DiagSum)
      TYPE(CMatrix), INTENT(IN) :: CMat
      DOUBLE COMPLEX DiagSum
      INTEGER kk
      
      PRINT *, "Computing the Trace of the matrix..."
      IF (CMat%N(1) .ne. CMat%N(2)) THEN
         PRINT *, "  Error: Trace is defined only for square matrices."
         RETURN
      ENDIF
      
      DiagSum = CMPLX(0.,0.)
      DO kk = 1,CMat%N(1)
         DiagSum = DiagSum + CMat%Elem(kk,kk)
      ENDDO
      
      RETURN
   END FUNCTION
   
!    FUNCTION Determinant(CMat) RESULT(Det)
!       TYPE(CMatrix), INTENT(IN) :: CMat
!       DOUBLE COMPLEX Det
!    
!       RETURN
!    END FUNCTION
   
   FUNCTION Adjoint(CMat) RESULT(AdjMat)
      TYPE(CMatrix), INTENT(IN) :: CMat
      DOUBLE COMPLEX, DIMENSION(CMat%N(2),CMat%N(1)) :: AdjMat
      
      PRINT *, "Computing the conjugate-transpose matrix..."
      AdjMat = CONJG(TRANSPOSE(CMat%Elem))
      
      RETURN
   END FUNCTION
   
   SUBROUTINE Init(CMat)
      IMPLICIT NONE
      TYPE(CMatrix) :: CMat
      DOUBLE PRECISION, DIMENSION(CMat%N(1),CMat%N(2)) :: Re, Im
      
      PRINT *, "Random initialization of the matrix..."
      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(Re)
      Re = Re*2 -1
      CALL RANDOM_NUMBER(Im)
      Im = Im*2 -1
      CMat%Elem = CMPLX(Re, Im)     
   
      RETURN
   END SUBROUTINE
   
   SUBROUTINE WriteOnFile(CMat, Outfile)
      IMPLICIT NONE
      TYPE(CMatrix) :: CMat
      CHARACTER(30) Outfile
      INTEGER ll
      
      PRINT *, "Writing Output on file..."
      
      Outfile = TRIM(Outfile)
      OPEN(UNIT=20, FILE=Outfile, STATUS='unknown')
      WRITE(20,'(A)') "### COMPLEX MATRIX DERIVED TYPE: ###"
      WRITE(20,"(A,I2,A,I2,A)") "MATRIX SIZE: (", CMat%N(1)," ,",CMat%N(2)," )"
      IF ( CMat%N(1) == CMat%N(2) ) THEN
         WRITE(20,"(A,('('sf6.3xspf6.3x'i)'))") "MATRIX TRACE: ", CMat%CTrace
      ENDIF 
      
      WRITE(20,*) NEW_LINE('A'), "MATRIX ELEMENTS: "
      DO ll = 1,CMat%N(2) 
         WRITE(20, "(*('('sf6.3xspf6.3x'i)':x))") CMat%Elem(:,ll)
      ENDDO
      
      WRITE(20,*) NEW_LINE('A'), "ADJOINT MATRIX: "
      DO ll = 1,CMat%N(1) 
         WRITE(20, "(*('('sf6.3xspf6.3x'i)':x))") CMat%ElemAdj(:,ll)
      ENDDO
      CLOSE(20)
      
      RETURN
   END SUBROUTINE
   
END MODULE ComplexMatrix
