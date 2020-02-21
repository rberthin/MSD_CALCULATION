MODULE DISTANCE_PBC

IMPLICIT NONE

CONTAINS
    DOUBLE PRECISION FUNCTION DISTANCE_PBC_OPT(coord_a,coord_b,cell_size) RESULT(res)

         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: coord_a
         DOUBLE PRECISION, INTENT(IN) :: coord_b
         DOUBLE PRECISION, INTENT(IN) :: cell_size 
         DOUBLE PRECISION :: dcf, halfboxrec
    
         halfboxrec = 2.0/cell_size
         dcf = coord_b-coord_a    
         ! MINIMAL DISTANCE
         dcf = dcf-cell_size*int(dcf*halfboxrec)
         res = dcf

    END FUNCTION DISTANCE_PBC_OPT 

END MODULE DISTANCE_PBC

