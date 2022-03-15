! Programmed by Robert Alvarez
! Last modified: March 15th, 2022
!
! Library to perform all cuts needed on a molecule to try to maxmize the volume of the smallest box
MODULE cut_molecule
    USE element_data_type, ONLY: elements, calc_distance, calc_location ! Class that define the attributes and function of the elements in a box
    USE sort, ONLY: merge_sort                                          ! Modify merge sort to work with the data type pass above
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: find_best_permutation, cuts_with_best_permutation

    CONTAINS

    !Generate all possible combinations, test them and return the sequence of character that would generate the smaller box with the largest volume 
    SUBROUTINE find_best_permutation(all_elements, n, best_overall_volume, best_overall_permutation)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)           ! Array with all elements inside of the box
        INTEGER, INTENT(IN) :: n                                    ! Number of concatenations
        REAL*8, INTENT(INOUT) :: best_overall_volume                ! Store largest volume 
        CHARACTER(*), INTENT(INOUT) :: best_overall_permutation     ! Permutation that makes the largest volume

        CHARACTER(3*n) :: curr_perm, concatenations(6**n)           ! Temporarily save a concatenation // Save all different possible concatenations
        INTEGER, ALLOCATABLE :: depth_permutations(:)               ! Parameter for find_all_permutations
        REAL*8 :: temp_volume, volumes(SIZE(all_elements))          ! Temporarily save a volume // Save all volumes
        INTEGER :: i, j

        ALLOCATE(depth_permutations(n))

        i=1
        CALL find_all_permutations(concatenations, 1, depth_permutations, i)    ! Generate all possible concatenations

        DEALLOCATE (depth_permutations)

        ! Try all concatenations and get the best one
        DO i=1, SIZE(concatenations)
            curr_perm = concatenations(i)       !Make a copy of the concatenation

            ! Heap sort may be better for sorting with respect to element # because elements are partially sorted
            CALL merge_sort(all_elements, 'x')  !Merge in ascending order with respect to x
            CALL merge_sort(all_elements, 'n')  !Merge in ascending order with respect to the element number

            j=0
            CALL this_permutation(all_elements, curr_perm, j)    ! Perform all of the cuts

            !Calculate all boxes volumes
            DO j=1, SIZE(volumes)
                volumes(j) = all_elements(j)%calculate_volume()
            END DO

            !Reset atoms walls
            DO j=1, SIZE(all_elements)
                CALL all_elements(j)%reset_walls()
            END DO

            ! Get the smallest volume
            temp_volume = MINVAL(volumes)

            !If the smallest volume make by curr_perm is bigger than the previous smallest-biggest volume
            IF (temp_volume - best_overall_volume > 0.00001) THEN
                best_overall_volume = temp_volume               ! Save new smallest-biggest volume
                best_overall_permutation = concatenations(i)    ! Save permutation take makes best_overall_volume
            END IF
        END DO
    END SUBROUTINE find_best_permutation

    ! Generate all possible combinations by making size(depth) concatenations of the strings in base_permutations
    ! It return an array with all the combinations generated
    RECURSIVE SUBROUTINE find_all_permutations(all_conc, loop_n, depth, array_location)
        IMPLICIT NONE
        CHARACTER(*), INTENT(INOUT) :: all_conc(:)  ! Save all concatenations that would be generated
        INTEGER, INTENT(IN) :: loop_n               ! Current recursion number
        INTEGER, INTENT(INOUT) :: depth(:)          ! Save the indices that would be used from base_permutations to generate a permutation
        INTEGER, INTENT(INOUT) :: array_location    ! Index to store new concatenation in all_conc

        CHARACTER(3) :: base_permutations(6) = ['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx']
        INTEGER :: i

        IF (loop_n < SIZE(depth)) THEN  !For every depth column except the last one...
            DO i=1, 6
                depth(loop_n) = i   ! Save that in recursion loop_n we were in column i (column means the index of base_permutations)
                CALL find_all_permutations(all_conc, loop_n+1, depth, array_location)    !Make a recursion for each depth column
            END DO
        ELSE                            !For the last depth column...
            DO i=1, 6
                depth(loop_n) = i   ! Save in last column of depth the index i
                BLOCK
                    CHARACTER(3*SIZE(depth)) :: String  ! Store the string generated by the permutation in base_permutations base on indixes of depth
                    INTEGER :: j, x

                    DO j=1, SIZE(depth)
                        x = 3*(j-1)
                        String(1+x:3+x) = base_permutations(depth(j))    !Concatenate xyz 3-characters permutations to String base on depth column
                    END DO

                    all_conc(array_location) = String   ! Save String denerated with the concatenation in all_conc

                    array_location = array_location + 1 ! Update index in where to store the next concatenation
                END BLOCK
            END DO
        END IF
    END SUBROUTINE find_all_permutations

    ! Use the sequence of character given in cut_sequence to make cuts in the box where all_elements are at
    RECURSIVE SUBROUTINE this_permutation(all_elements, cut_sequence, no_cuts)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)   ! All elements inside of a box
        CHARACTER(*), INTENT(INOUT) :: cut_sequence         ! Sequence of 'xyz' characters to dictate where to make the cut
        INTEGER, INTENT(INOUT) :: no_cuts                   ! Keep track of how many times we haven't make a cut because we don't want to cut an element in half

        INTEGER, ALLOCATABLE :: best_cuts(:)                ! Save index of the element closest from the right to the cut for the new box
                                                            ! It is an array because mulitple good cut may be found

        IF (SIZE(all_elements) == 1) THEN
            !Everything is set
        ELSE
            ! Find the best location to make a cut (populate best_cuts)
            IF (no_cuts >= 5) THEN                       !5 no_cuts warranty that every coordinate has been tried at least once
                CALL get_best_cuts(all_elements,cut_sequence(1:1), best_cuts, .TRUE.)   ! Get a cut even if it is a bad one.
            ELSE
                CALL get_best_cuts(all_elements,cut_sequence(1:1), best_cuts)           ! Get the best_cuts
            END IF

            cut_sequence = next_sequence(cut_sequence)  ! Shift every character to the left and the first character to the end

            IF (ALLOCATED(best_cuts)) THEN              ! If a good cut was found
                no_cuts = 0
                CALL make_n_cuts(all_elements, best_cuts, cut_sequence, no_cuts)    ! Make those cuts
                DEALLOCATE(best_cuts)
            ELSE                                        ! If no good cut was found
                no_cuts = no_cuts + 1
                CALL this_permutation(all_elements, cut_sequence, no_cuts)          ! Try to find good cuts with the next coordinate
            END IF
        END IF
    END SUBROUTINE this_permutation

    ! Perform all cuts and send the new boxes to get more
    SUBROUTINE make_n_cuts(elements_array, cuts_locations, cut_seq, no_cut)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: elements_array(:) ! All elements inside of the box
        INTEGER, INTENT(IN) :: cuts_locations(:)            ! Index(es) of the right most element after the cut
        CHARACTER(*), INTENT(INOUT) :: cut_seq              ! Sequence of character use to cut the molecule
        INTEGER, INTENT(INOUT) :: no_cut                    ! Number of time cuts has not been perform

        CHARACTER :: coordinate
        INTEGER :: i, j
        REAL*8 :: cut_location

        coordinate = cut_seq(LEN(cut_seq):LEN(cut_seq))     ! Get the coordinate that was used to find the cuts_locations

        ! For every index found to make a molecule the right most molecule, make a cut
        DO i=1, SIZE(cuts_locations)
            cut_location = calc_location(elements_array(cuts_locations(i)+1), elements_array(cuts_locations(i)), coordinate)    ! Get coordinate value of the cut

            !Set walls and send new boxes to get more cuts
            IF (i+1 <= SIZE(cuts_locations)) THEN                                   ! More than one cut is left
                IF (i==1) THEN                                                          ! Cut one when cut 2 exist
                    ! Update walls of the elements in the first and second box from left to right
                    DO j=1, cuts_locations(i+1)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(1:cuts_locations(i)), cut_seq, no_cut) ! Send the left most box to get more cuts
                ELSE                                                                    ! Cut i when cut i+1 exist
                    ! Update walls of the elements in the box i-1, i and i+1
                    DO j=cuts_locations(i-1)+1, cuts_locations(i+1)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(cuts_locations(i-1)+1:cuts_locations(i)), cut_seq, no_cut) ! Send the box i-1 to get more cuts
                END IF
            ELSE                                                                    ! One cut is left
                IF (i==1) THEN                                                          ! Cut one when cut 2 not exist
                    ! Update walls for the only two boxes generated (only one cut)
                    DO j=1, SIZE(elements_array)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(1:cuts_locations(i)), cut_seq, no_cut)     ! Send box 1 to get more cuts
                    CALL this_permutation(elements_array(cuts_locations(i)+1:SIZE(elements_array)), cut_seq, no_cut)    ! Send box 2 to get more cuts
                ELSE                                                                    ! Last cut when cut i-1 exist
                    ! Update walls of elements in the right most box
                    DO j=cuts_locations(i-1)+1, SIZE(elements_array)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(cuts_locations(i-1)+1:cuts_locations(i)), cut_seq, no_cut)     ! Send the second to last box to get more cuts
                    CALL this_permutation(elements_array(cuts_locations(i)+1:SIZE(elements_array)), cut_seq, no_cut)    ! Send the last box to get more cuts
                END IF
            END IF
        END DO
    END SUBROUTINE make_n_cuts

    ! Find the largest gaps between elements in the coordinate given, save the index of the left element that make the largest gap found
    ! and save the index in best_location
    SUBROUTINE get_best_cuts(all_ele, coordinate, best_location, make_cut)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_ele(:)                ! All elements in the box
        CHARACTER, INTENT(IN) :: coordinate                         ! Coordinate that would be use as reference to find the largest gaps
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: best_location(:)     ! Save idexes for the molecuel that would becomes the right most molecule
                                                                    ! after a cut is made in the lagest gap
        LOGICAL, OPTIONAL, INTENT(IN) :: make_cut                   ! If optional parameter is sent as .True. a cut would be make regardless if it is a good cut or not

        INTEGER :: temp_best_location(SIZE(all_ele)-1)              ! Worst case scenario (in space complexity) all gapes are the largest gaps
        INTEGER :: i, n                                             ! n is the number of largest gaps with the same distance
        REAL*8 :: distances(SIZE(all_ele)-1)                        ! Save gap values for all neighbor atoms
        REAL*8 :: largest_distance                                  ! Save largest gap value

        !Calculate distances between neighbor atoms
        CAll merge_sort(all_ele, coordinate)

        ! Calculate all distances between neighbor atoms and find the largest distance
        largest_distance = -HUGE(largest_distance)
        DO i=1, SIZE(all_ele)-1
            distances(i) = calc_distance(all_ele(i+1), all_ele(i), coordinate)  ! Calculate distances between neighbor atoms
            IF (distances(i) > largest_distance) THEN   ! Find largest distance
                largest_distance = distances(i)
            END IF
        END DO

        ! Save atom positions in best_location if a good cut is found or if make_cut is passed
        IF (largest_distance > -0.8 .OR. PRESENT(make_cut)) THEN        ! largest_distance > -0.8 is to avoid to cut a molecule in half
            ! Look for symmetry, if any, to save multiple indixes
            n=0
            DO i=1, SIZE(distances)
                IF (largest_distance - distances(i) <= 0.0001) THEN     ! Consideration of symmetry
                    n = n+1
                    temp_best_location(n) = i
                END IF
            END DO

            ALLOCATE(best_location(n))      ! Make space in best_location to save indexes values

            ! Save indices in best_location
            DO i=1, n
                best_location(i) = temp_best_location(i)
            END DO
        END IF
    END SUBROUTINE get_best_cuts

    ! Receive a sequence of characters, shift all the characters one to the left, send the last one to the end of the sequence and return it
    PURE FUNCTION next_sequence(old_sequence)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: old_sequence        ! Pass string
        CHARACTER(LEN(old_sequence)) :: next_sequence   ! New string
        CHARACTER :: char_temp                          ! Save first character in old_squence
        INTEGER :: i

        char_temp = old_sequence(1:1)                   ! Save first character
        DO i=1, LEN(old_sequence)-1                     ! Copy all characters one shift to the left
            next_sequence(i:i) = old_sequence(i+1:i+1)
        END DO
        next_sequence(LEN(old_sequence):LEN(old_sequence)) = char_temp  ! Paste last character
    END FUNCTION next_sequence

    ! Perform the cut sequence with the best permutation
    SUBROUTINE cuts_with_best_permutation(all_elements, best_permutation)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)   ! All elements in the box
        CHARACTER(*), INTENT(INOUT) :: best_permutation     ! Best sequence of character to maximize the volume of the smaller box
        INTEGER :: i

        i=0
        CALL this_permutation(all_elements, best_permutation, i)   ! Perform cuts and with the best permutation sequence
    END SUBROUTINE cuts_with_best_permutation
END MODULE cut_molecule
