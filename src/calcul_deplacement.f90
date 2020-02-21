PROGRAM CALCUL_DEPLACEMENT
   
   USE DISTANCE_PBC
   USE NUM_LINES_FILE
   USE CONSTANTS
   
   IMPLICIT NONE
   
   CHARACTER(LEN=200) :: crap
   CHARACTER(LEN=10) :: nom
   INTEGER :: nb_lines, nb_species
   INTEGER :: w, j, n, k, m, q, p, i
   INTEGER :: io, element, id, nb_configs
   DOUBLE PRECISION :: boxx, boxy, boxz
   DOUBLE PRECISION, DIMENSION(3) :: distance
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: CY, CY2
   INTEGER, ALLOCATABLE, DIMENSION(:) :: num_atoms, ref_atom, num_molecules
 
   OPEN(UNIT = 11, FILE = "deplacement.inpt", STATUS='old', IOSTAT=io)

   ! READ DEPLACEMENT.INPT
   !-----------------------
   READ(11,*) boxx
   boxx = boxx*bohr
   READ(11,*) boxy
   boxy = boxy*bohr
   READ(11,*) boxz
   boxz = boxz*bohr
   READ(11,*) nb_configs
   READ(11,*) nb_species
   READ(11,*) !read comment
   ALLOCATE(num_atoms(nb_species))
   ALLOCATE(ref_atom(nb_species))
   ALLOCATE(num_molecules(nb_species))
   DO i = 1, nb_species
     READ(11,*) num_molecules(i)
     READ(11,*) num_atoms(i)
     READ(11,*) ref_atom(i)
   ENDDO
   !-----------------------
   ! NB LIGNES TRAJECTOIRE
   !------------------------------------
   nb_lines = 0
   OPEN(unit = 10, file = "trajectories.xyz", status='old', iostat=io)
   nb_lines = count_lines(10)
   !------------------------------------
   !nb_configs = nb_lines/(nb_species_tot+2)
   !write(6,*) nb_configs
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
         distance(1) = DISTANCE_PBC_OPT(CY(m,1), CY2(m,1), boxx)
         distance(2) = DISTANCE_PBC_OPT(CY(m,2), CY2(m,2), boxy)
         distance(3) = DISTANCE_PBC_OPT(CY(m,3), CY2(m,3), boxz)
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
