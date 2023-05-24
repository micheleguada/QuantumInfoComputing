MODULE Debugger
    IMPLICIT NONE    
    CHARACTER(len=*), PARAMETER :: ErrorFile = "FATAL_ERRORS.txt"
    
    CONTAINS
    
    !This subroutine prints a message at the beginning of the program if Debug=.TRUE.
    SUBROUTINE InitDebugger(Debug)
        LOGICAL Debug
        
        IF (Debug .eqv. .TRUE.) THEN
            WRITE(*,"(A)")   "#### DEBUGGER: ####"
            WRITE(*,"(A)")   "    Debugger is active. It prints messages to help finding bugs and errors."
            WRITE(*,"(A)")   "    Lines starting with: '-->' are produced by debugger subroutines."
            WRITE(*,"(A,A)") "    Fatal error messages are printed in the file: ", ErrorFile
            WRITE(*,"(A)")   "###################"
        ENDIF
    
        RETURN
    END SUBROUTINE
        
    
    !This subroutine is called repeatedly during the execution as a checkpoint
    SUBROUTINE Checkpoint(Debug, Text, RNumber, LNumber) 
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
        
        IF (PRESENT(LNumber)) THEN
            WRITE(*,"(A,A,A,I4,A)") "ERROR: ", Text, " [Line # ", LNumber, "]"
            IF (Fatal) THEN
                OPEN(unit=90, file=ErrorFile, status='unknown')
                WRITE(90,"(A,A,I4,A)") Text, " [Line # ", LNumber, "]"
                WRITE(*,"(A)") "FATAL error! Stopping..."
                CLOSE(90)
                STOP
            ENDIF
        ELSE
            WRITE(*,"(A,A)") "ERROR: ", Text
            IF (Fatal) THEN
                OPEN(unit=90, file=ErrorFile, status='unknown')
                WRITE(90,"(A)") Text
                WRITE(*,"(A)") "FATAL error! Stopping..."
                CLOSE(90)
                STOP
            ENDIF
        ENDIF
        
        RETURN
    END SUBROUTINE
    
END MODULE
