PROGRAM TestComplexMatrix
   USE ComplexMatrix
   IMPLICIT NONE
   CHARACTER(30) filename
   CHARACTER(6) arg2, arg3
   INTEGER, DIMENSION(2) :: xy
   TYPE(CMatrix) :: AA 
   
   CALL get_command_argument(1, filename)
   IF (LEN_TRIM(filename) == 0) THEN
      filename = "TestOutput.txt"
   ENDIF
   !Setting matrix's size from arguments or with default values
   CALL get_command_argument(2, arg2)
   CALL get_command_argument(3, arg3)
   IF ((LEN_TRIM(arg2) == 0) .or. (LEN_TRIM(arg3) == 0)) THEN
      PRINT *, "Matrix size not provided; using default values (4,4)"
      xy = (4,4) !default values
   ELSE
      PRINT *, "Provided Matrix size: (",TRIM(arg2),",",TRIM(arg3),")"
      READ(arg2,*) xy(1)
      READ(arg3,*) xy(2)      
   ENDIF
   
   AA%N = xy
   ALLOCATE(AA%Elem(xy(1),xy(2)))
   !Randomly initializing the matrix
   CALL Init(AA)
   !Computing the trace 
   AA%CTrace = .Tr.AA
   !Computing the Adjoint
   ALLOCATE(AA%ElemAdj(xy(2),xy(1)))
   AA%ElemAdj = .Adj.AA
   
!    AA%CDet = Determinant(AA)

   !Write on file the CMatrix object
   CALL WriteOnFile(AA, filename)

   DEALLOCATE(AA%Elem)
   DEALLOCATE(AA%ElemAdj)
   STOP
END PROGRAM
