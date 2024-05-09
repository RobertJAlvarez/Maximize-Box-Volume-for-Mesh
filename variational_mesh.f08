MODULE variational_mesh
   USE element_data_type, ONLY: DBL, elements, makeCube
   USE cut_molecule, ONLY: find_best_permutation, cuts_with_best_permutation
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: variationalMesh

CONTAINS

   SUBROUTINE variationalMesh
      !Declare file variables
      CHARACTER(len=80) :: msg
      INTEGER :: IOStatus

      !Declare program variables
      TYPE(elements), ALLOCATABLE :: R(:) !"allocatable" is use when size is not known before the program runs.
      INTEGER :: nAtoms, i, arrayStatus

      !Open file and catch any error that may happen
      OPEN (UNIT=3, FILE='XMOL.DAT', STATUS='OLD', ACTION='READ', IOSTAT=IOStatus, IOMSG=msg)

      !Read file and save data
      IF (IOStatus == 0) THEN
         READ (3, *, IOSTAT=IOStatus) nAtoms
         READ (3, *, IOSTAT=IOStatus)

         ALLOCATE (R(nAtoms), STAT=arrayStatus, ERRMSG=msg) !Initialize arrays size

         IF (arrayStatus /= 0) THEN
            WRITE (*, 1020) arrayStatus
1020        FORMAT('VariationalMesh: Error allocating R, IOSTAT', I6)
            WRITE (*, *) TRIM(msg)
            STOP
         END IF

         DO i = 1, nAtoms
            BLOCK !Read and store atoms data
               CHARACTER(100) :: line_info
               CHARACTER(4) :: char_element = ''
               INTEGER :: int_element
               REAL(DBL) :: x, y, z

               READ (3, *, IOSTAT=IOStatus) int_element, x, y, z

               IF (IOStatus == 0) THEN                          !If the element was a number
                  CALL R(i)%set_molecule(int_element, x, y, z)
               ELSE IF (IOStatus > 0) THEN                      !If the element was a letter
                  BACKSPACE (UNIT=3)

                  READ (3, "(A100)", IOSTAT=IOStatus) line_info   !Read all data from this line
                  line_info = TRIM(line_info)                     !Delete firsts and lasts blank spaces

                  READ (line_info, *, IOSTAT=IOStatus) char_element, x, y, z
                  CALL R(i)%set_molecule(TRIM(char_element), x, y, z)
               ELSE IF (IOStatus < 0) THEN                          !End of file had been reach
                  WRITE (*, *) 'End of file reach before expecting' !Make something to handle this situation
               END IF
            END BLOCK

            CALL R(i)%set_radius()
            CALL R(i)%reset_walls()
         END DO
      ELSE
         WRITE (*, 1030) IOStatus
1030     FORMAT('Error opening file: IOSTAT = ', I6)
         WRITE (*, *) TRIM(msg)
         STOP
      END IF

      CLOSE (3)

      make_best_cuts: BLOCK
         INTEGER, PARAMETER :: max_n = 5  !"n" stands for the number of concatenations
         CHARACTER(max_n*3) :: best_n_permutation
         INTEGER :: n_concatenations
         REAL(DBL) :: best_n_volume

         best_n_volume = -HUGE(best_n_volume)

         !Find permutation with the best cuts
         DO n_concatenations = int(max_n/2.0), max_n
            CALL find_best_permutation(R, n_concatenations, best_n_volume, best_n_permutation)
         END DO

         best_n_permutation = TRIM(best_n_permutation)

         !Use best permutation to set the walls
         CALL cuts_with_best_permutation(R, best_n_permutation)
      END BLOCK make_best_cuts

      !Create cubes from boxes
      make_cubes: BLOCK
         TYPE(elements) :: cubes(SIZE(R))
         TYPE(elements) :: emptyBox(SIZE(R)*5)
         INTEGER :: nBoxes

         nBoxes = 1
         cubes = R

         DO i = 1, SIZE(cubes)
            CALL cubes(i)%makeCube(emptyBox, nBoxes)
         END DO
         nBoxes = nBoxes - 1

         DEALLOCATE (R)
         ALLOCATE (R(SIZE(cubes) + nBoxes), STAT=arrayStatus, ERRMSG=msg)

         IF (arrayStatus /= 0) THEN
            WRITE (*, 1040) arrayStatus
1040        FORMAT('VariationalMesh: Error allocating R in make_cubes block, IOSTAT', I6)
            WRITE (*, *) TRIM(msg)
            STOP
         END IF

         R(:SIZE(cubes)) = cubes
         R(SIZE(cubes) + 1:nBoxes) = emptyBox(:nBoxes)
      END BLOCK make_cubes

      OPEN (UNIT=4, FILE='BoxInfoOutput.DAT', STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOStatus, IOMSG=msg)

      IF (IOStatus == 0) THEN
         !Print molecule information
         WRITE (4, *) nAtoms, SIZE(R)  !Number of boxes
         DO i = 1, SIZE(R)
            CALL R(i)%print_info(4)
         END DO
      ELSE
         WRITE (*, 1050) IOStatus
1050     FORMAT('Error creating/writing on BoxInfoOutput.DAT: IOSTAT = ', I6)
         WRITE (*, *) TRIM(msg)
         STOP
      END IF

      CLOSE (4)

      DEALLOCATE (R)
   END SUBROUTINE variationalMesh
END MODULE variational_mesh
