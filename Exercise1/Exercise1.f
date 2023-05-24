      MODULE Intro
         INTEGER*2 ii
         INTEGER*2 jj
      END MODULE Intro
  
      MODULE Advanced
         CONTAINS
         FUNCTION sumab(aa, bb) RESULT(res)
            INTEGER*2 aa
            INTEGER*2 bb
            INTEGER*2 res
            res = aa + bb   
            RETURN
         END FUNCTION
      END MODULE Advanced
  
      PROGRAM Exercise1
      USE Intro
      USE Advanced
      IMPLICIT NONE
         INTEGER*2 somma
         ii = 4
         jj = 3
         somma = sumab(ii, jj)
         PRINT *, 'Sum is ', somma
      STOP
      END PROGRAM Exercise1

c     Output of the program:
c      Sum is       7
