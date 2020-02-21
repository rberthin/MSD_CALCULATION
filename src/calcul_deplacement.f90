PROGRAM CALCUL_DEPLACEMENT
   
   USE DISTANCE_PBC
   USE NUM_LINES_FILE
   USE CONSTANTS
   
   IMPLICIT NONE
   
   CHARACTER*200 :: crap
   CHARACTER(LEN=10) :: nom
   INTEGER :: nb_lines, nb_species_tot, i, nb_atomes_par_conf
   INTEGER :: w, j, element, id, n, nb_configs, k, m, q, p
   INTEGER :: io
   DOUBLE PRECISION :: x, y, z
   DOUBLE PRECISION, DIMENSION(3) :: box_size, coord_1, coord_2
   INTEGER, ALLOCATABLE, DIMENSION (:) :: nb_part
   DOUBLE PRECISION, DIMENSION(3) :: distance
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: CY, CY2
   
   
   ! PARAMETRE DE LA BOITE DE SIMULATION
   !------------------------------------
   x = 50.175232345*bohr
   y = 50.175232345*bohr
   z = 50.175232345*bohr
   nb_species_tot = 1296
   !------------------------------------
   
   ! NB LIGNES TRAJECTOIRE
   !------------------------------------
   nb_lines = 0
   OPEN(unit = 10, file = "trajectories.xyz", status='old', iostat=io)
   nb_lines = count_lines(10)
   !------------------------------------
     
   nb_configs = nb_lines/(nb_species_tot+2)
   write(6,*) nb_configs
   OPEN(unit = 20, file = "deplacements.out")
   
   ! LECTURE DE LA TRAJECTOIRE
   !------------------------------------
   
   ALLOCATE(CY(216,3))
   ALLOCATE(CY2(216,3))
   
   READ(10,*) crap  
   READ(10,*) crap
   
   DO n = 1, 216
      READ(10,*) crap
   ENDDO
   
   DO m = 1, 216
      READ(10,*) nom, CY(m,1), CY(m,2), CY(m,3)
   ENDDO
   
   DO p = 1,  (nb_species_tot-432)
      READ(10,*) crap
   ENDDO
   
   DO w = 1, nb_configs-1
      
      READ(10,*) crap
      READ(10,*) crap
   
      write(6,*) 'STEP', w
      DO n = 1, 216
         READ(10,*)
      ENDDO
      DO m = 1, 216
         READ(10,*) nom, CY2(m,1), CY2(m,2), CY2(m,3)
         distance(1) = DISTANCE_PBC_OPT(CY(m,1), CY2(m,1), x)
         distance(2) = DISTANCE_PBC_OPT(CY(m,2), CY2(m,2), y)
         distance(3) = DISTANCE_PBC_OPT(CY(m,3), CY2(m,3), z)
         WRITE(20,*) distance(1), distance(2), distance(3)
         CY(m,1) = CY2(m,1)
         CY(m,2) = CY2(m,2)
         CY(m,3) = CY2(m,3)           
      ENDDO
      DO p = 1, (nb_species_tot-432)  
         READ(10,*) crap
      ENDDO
   ENDDO
   WRITE(6,*) ' **** ALL DONE !! ****'
   CLOSE(10)
   CLOSE(20)

END PROGRAM CALCUL_DEPLACEMENT
