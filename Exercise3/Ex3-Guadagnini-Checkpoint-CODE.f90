MODULE Checkpoints
    IMPLICIT NONE    
    CONTAINS
    
    SUBROUTINE Checkpoint(Debug, Text, INumber, RNumber) 
        LOGICAL Debug
        CHARACTER(len=*), OPTIONAL :: Text
        INTEGER, OPTIONAL :: INumber    !optional line number variable
        REAL, OPTIONAL :: RNumber       !optional real variable
        
        IF (Debug .eqv. .TRUE.) THEN
            IF (PRESENT(Text) .and. PRESENT(INumber)) THEN
                WRITE(*,"(A,A,A,I4,A)") "--> ", Text, " [Line #", INumber, "]"
            ELSEIF (PRESENT(Text)) THEN
                WRITE(*,"(A,A)") "--> ", Text
            ENDIF
            
            IF (PRESENT(RNumber)) THEN
                WRITE(*,"(A,f8.5)") "--> ", RNumber
            ENDIF
        ENDIF
        
        RETURN
    END SUBROUTINE

END MODULE
