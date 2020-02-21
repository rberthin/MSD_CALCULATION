REAL(8) FUNCTION dist_w_PBC_opt(coor_a,coor_b,cell_size) RESULT(res)
         IMPLICIT NONE
         REAL(8), INTENT(IN) :: coor_a
         REAL(8), INTENT(IN) :: coor_b
         REAL(8), INTENT(IN) :: cell_size
 
         integer :: i
         REAL(8), DIMENSION(9,3) :: coor_b_rep
         REAL(8), DIMENSION(9) :: dists
         REAL(8) :: dxcf, dycf,dzcf, halfboxxrec, halfboxyrec
         REAL(8) :: halfboxzrec
 
         halfboxxrec=2.0/cell_size
 
         dxcf=coor_b-coor_a
 
         ! minimal distance convenction
         dxcf=dxcf-cell_size*int(dxcf*halfboxxrec)

         res=dxcf
END FUNCTION dist_w_PBC_opt 
 
INTEGER FUNCTION num_lines_file(num_file) RESULT(res)
        IMPLICIT NONE
        INTEGER, INTENT(IN) ::num_file
        INTEGER :: io
        CHARACTER*200 :: inputline
        INTEGER :: n_lin
        n_lin=0
        DO
           READ(num_file,*,IOSTAT=io)  inputline
           IF (io > 0) THEN
              WRITE(*,*) 'Check input.  Something was wrong'
              EXIT
           ELSE IF (io < 0) THEN
              EXIT
           ELSE
              n_lin = n_lin + 1
           ENDIF
        ENDDO
        REWIND(num_file)
        res=n_lin
END FUNCTION num_lines_file


PROGRAM DEPLACEMENTS 

IMPLICIT NONE

INTEGER, EXTERNAL :: num_lines_file
real(8), EXTERNAL :: dist_w_PBC_opt

CHARACTER*200 :: crap
CHARACTER(LEN=10) :: nom
INTEGER :: nb_lines, nb_species_tot, i, nb_atomes_par_conf
INTEGER :: w, j, element, id, n, nb_configs, k, m, q, p
INTEGER :: io
DOUBLE PRECISION :: bohr, x, y, z
DOUBLE PRECISION, DIMENSION(3) :: box_size, coord_1, coord_2
INTEGER, ALLOCATABLE, DIMENSION (:) :: nb_part
DOUBLE PRECISION, DIMENSION(3) :: distance
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: CY, CY2

bohr = 0.529177

! PARAMETRE DE LA BOITE DE SIMULATION
!------------------------------------
x = 38.834749574*bohr
y = 38.834749574*bohr
z = 38.834749574*bohr
nb_species_tot = 594
!------------------------------------

! NB LIGNES TRAJECTOIRE
!------------------------------------
nb_lines = 0
OPEN(unit = 10, file = "trajectories.xyz", status='old', iostat=io)
nb_lines = num_lines_file(10)
!------------------------------------
  
nb_configs = nb_lines/(nb_species_tot+2)
write(6,*) nb_configs
OPEN(unit = 20, file = "deplacements.out")
OPEN(unit = 21, file = "test.out")

! LECTURE DE LA TRAJECTOIRE
!------------------------------------

ALLOCATE(CY(95,3))
ALLOCATE(CY2(95,3))

READ(10,*) crap  
READ(10,*) crap
WRITE(21,*) crap

DO n = 1, 119
   READ(10,*) crap
ENDDO

DO m = 1, 95
   READ(10,*) nom, CY(m,1), CY(m,2), CY(m,3)
   !WRITE(21,*) CY(m,1), CY(m,2), CY(m,3)
ENDDO

DO p = 1,  (nb_species_tot-119-95)
   READ(10,*) crap
ENDDO

DO w = 1, nb_configs-1
   
   READ(10,*) crap
   READ(10,*) crap
   WRITE(21,*) crap

   write(6,*) 'STEP', w
   DO n = 1, 119
      READ(10,*)
   ENDDO
   DO m = 1, 95
      READ(10,*) nom, CY2(m,1), CY2(m,2), CY2(m,3)
      !WRITE(21,*) CY2(m,1), CY2(m,2), CY2(m,3)
      distance(1) = dist_w_PBC_opt(CY(m,1), CY2(m,1), x)
      distance(2) = dist_w_PBC_opt(CY(m,2), CY2(m,2), y)
      distance(3) = dist_w_PBC_opt(CY(m,3), CY2(m,3), z)
      WRITE(20,*) distance(1), distance(2), distance(3)
      CY(m,1) = CY2(m,1)
      CY(m,2) = CY2(m,2)
      CY(m,3) = CY2(m,3)           
   ENDDO
   DO p = 1, (nb_species_tot-119-95)  
      READ(10,*) crap
   ENDDO
ENDDO
WRITE(6,*) ' **** ALL DONE !! ****'
CLOSE(10)
CLOSE(20)

END PROGRAM
