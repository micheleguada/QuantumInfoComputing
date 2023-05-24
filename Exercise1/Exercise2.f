      MODULE Vars
         INTEGER*2 ii
         INTEGER*2 jj
         INTEGER*4 aa
         INTEGER*4 bb
         REAL*4 pp
         REAL*4 qq
         REAL*8 tt
         REAL*8 ss 
         REAL*4, PARAMETER :: Pi4 = 3.14159265358979323846
         REAL*8, PARAMETER :: Pi8 = 3.14159265358979323846
      END MODULE

      PROGRAM SumNumbers
         USE Vars
         IMPLICIT NONE
         INTEGER*2 SumInt2
         INTEGER*4 SumInt4
         REAL*4 SumReal4
         REAL*8 SumReal8
c     Integers summing
         ii = 2000000               !overflow
         jj = 1
         SumInt2 = ii+jj
         PRINT *, 'Sum with Int*2 is ', SumInt2 
         aa = 2000000
         bb = 1
         SumInt4 = aa+bb
         PRINT *, 'Sum with Int*4 is ', SumInt4
c     Reals summing
         pp = Pi4*1e32
         qq = SQRT(2e0)*1e21
         SumReal4 = pp+qq
         PRINT *, 'Sum with Real*4 is ', SumReal4
         tt = Pi8*1d32
         ss = DSQRT(2d0)*1d21
         SumReal8 = tt+ss
         PRINT *, 'Sum with Real*8 is ', SumReal8
      STOP
      END PROGRAM

c     Output of the program:
c      Sum with Int*2 is  -31615
c      Sum with Int*4 is      2000001
c      Sum with Real*4 is    3.14159278E+32
c      Sum with Real*8 is    3.1415927410267153E+032

         


      
