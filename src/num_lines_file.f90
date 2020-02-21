MODULE NUM_LINES_FILE

IMPLICIT NONE
CONTAINS

    INTEGER FUNCTION count_lines(unit_of_file) RESULT(res)
         IMPLICIT NONE
         INTEGER, INTENT(IN) ::unit_of_file
         INTEGER :: io, n
         CHARACTER*200 :: line

         n = 0
         DO
           READ(unit_of_file,*,IOSTAT=io) line
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              EXIT
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              n = n + 1
           ENDIF
         ENDDO
         REWIND(unit_of_file)
         res=n
    END FUNCTION count_lines

END MODULE NUM_LINES_FILE
