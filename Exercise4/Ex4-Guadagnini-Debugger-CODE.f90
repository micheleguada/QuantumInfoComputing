MODULE Debugger
    IMPLICIT NONE    
    CHARACTER(len=*), PARAMETER :: ErrorFile = "FATAL_ERRORS.txt"
    
    CONTAINS
    
    !This subroutine is called repeatedly during the execution as a checkpoint
    SUBROUTINE Checkpoint(Debug, Text, LNumber, RNumber) 
        LOGICAL Debug
        CHARACTER(len=*), OPTIONAL :: Text
        INTEGER, OPTIONAL :: LNumber    !optional line number variable
        REAL, OPTIONAL :: RNumber       !optional real variable
        
        IF (Debug .eqv. .TRUE.) THEN
            IF (PRESENT(Text) .and. PRESENT(LNumber)) THEN
                WRITE(*,"(A,A,A,I4,A)") "--> ", Text, " [Line # ", LNumber, "]"
            ELSEIF (PRESENT(Text)) THEN
                WRITE(*,"(A,A)") "--> ", Text
            ENDIF
            
            IF (PRESENT(RNumber)) THEN
                WRITE(*,"(A,f8.5)") "--> ", RNumber
            ENDIF
        ENDIF
        
        RETURN
    END SUBROUTINE
    
    !This subroutine is called when an error occurs and prints a message on a file.
    ! An external script can read the message and deal with the error.
    ! The logical variable "Fatal" allows to stop the program for fatal errors.
    SUBROUTINE CatchError(Text, Fatal, LNumber)
        CHARACTER(len=*) Text
        LOGICAL Fatal
        INTEGER, OPTIONAL :: LNumber       !optional line number variable
        INTEGER ios
        
        OPEN(unit=90, file=ErrorFile, status='unknown')
        
        IF (PRESENT(LNumber)) THEN
            WRITE(*,"(A,A,A,I4,A)") "ERROR: ", Text, " [Line # ", LNumber, "]"
            IF (Fatal) THEN
                WRITE(90,"(A,A,I4,A)") Text, " [Line # ", LNumber, "]"
                WRITE(*,"(A)") "FATAL error! Stopping..."
                STOP
            ENDIF
        ELSE
            WRITE(*,"(A,A)") "ERROR: ", Text
            IF (Fatal) THEN
                WRITE(90,"(A)") Text
                WRITE(*,"(A)") "FATAL error! Stopping..."
                STOP
            ENDIF
        ENDIF
        
        CLOSE(90)
        
        RETURN
    END SUBROUTINE
    
END MODULE
