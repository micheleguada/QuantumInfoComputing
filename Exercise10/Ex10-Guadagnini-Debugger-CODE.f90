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
        
        !resetting fatal error file
        OPEN(unit=90, file=ErrorFile, status='replace')
        CLOSE(90)
    
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
    SUBROUTINE CatchError(Condition, Text, Fatal, LNumber)
        CHARACTER(len=*) Text
        LOGICAL Fatal, Condition
        INTEGER, OPTIONAL :: LNumber       !optional int variable
    
        IF (Condition) THEN
            IF (PRESENT(LNumber)) THEN
                WRITE(*,"(A,A,A,I4)") "ERROR: ", Text, " ", LNumber
            ELSE
                WRITE(*,"(A,A)") "ERROR: ", Text
            ENDIF
            
            IF (Fatal) THEN
                OPEN(unit=90, file=ErrorFile, status='unknown')
                IF (PRESENT(LNumber)) THEN
                    WRITE(90,"(A,A,I4)") Text, " ", LNumber
                ELSE 
                    WRITE(90,"(A)") Text
                ENDIF
                CLOSE(90)
                WRITE(*,"(A)") "FATAL error! Stopping..."
                STOP
            ENDIF
        ENDIF
        
        RETURN
    END SUBROUTINE
    
END MODULE
