MODULE cut_molecule
    USE element_data_type, ONLY: elements, calc_distance, calc_location
    USE sort, ONLY: merge_sort
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: find_best_permutation, cuts_with_best_permutation

    CONTAINS

    SUBROUTINE find_best_permutation(all_elements, n, best_overall_volume, best_overall_permutation)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(INOUT) :: best_overall_volume
        CHARACTER(*), INTENT(INOUT) :: best_overall_permutation

        CHARACTER(3*n) :: temp_String, concatenations(6**n)
        INTEGER, ALLOCATABLE :: depth_permutations(:)
        REAL*8 :: temp_volume, volumes(SIZE(all_elements))
        INTEGER :: i, j

        ALLOCATE(depth_permutations(n))

        i=1
        CALL find_all_permutations(concatenations, 1, depth_permutations, i)

        DEALLOCATE (depth_permutations)

        DO i=1, SIZE(concatenations)
            temp_String = concatenations(i)     !Make a copy of the concatenation

            !Heap sort may be better for sorting with respect to element # because elements are partially sorted
            CALL merge_sort(all_elements, 'x')  !Merge in ascending order with respect to x
            CALL merge_sort(all_elements, 'n')  !Merge in ascending order with respect to the element number

            j=0
            CALL this_permutation(all_elements, temp_String, j)    !Do all of the cuts

            !Calculate all boxes volumes
            DO j=1, SIZE(volumes)
                volumes(j) = all_elements(j)%calculate_volume()
            END DO

            !Rest atoms walls
            DO j=1, SIZE(all_elements)
                CALL all_elements(j)%reset_walls()
            END DO

            ! Get the smallest volume
            temp_volume = MINVAL(volumes)

            !If the smallest volume of String is bigger than the previous smallest-biggest volume
            IF (temp_volume - best_overall_volume > 0.00001) THEN
                best_overall_volume = temp_volume
                best_overall_permutation = concatenations(i)
            END IF
        END DO
    END SUBROUTINE find_best_permutation

    RECURSIVE SUBROUTINE find_all_permutations(all_conc, loop_n, depth, array_location)
        IMPLICIT NONE
        CHARACTER(*), INTENT(INOUT) :: all_conc(:)
        INTEGER, INTENT(IN) :: loop_n
        INTEGER, INTENT(INOUT) :: depth(:)
        INTEGER, INTENT(INOUT) :: array_location

        CHARACTER(3) :: base_permutations(6) = ['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx']

        INTEGER :: i

        IF (loop_n < SIZE(depth)) THEN  !For every depth column except the last one...
            DO i=1, 6
                depth(loop_n) = i
                CALL find_all_permutations(all_conc, loop_n+1, depth, array_location)    !Make a recursion for each depth column
            END DO
        ELSE                            !For the last depth column...
            DO i=1, 6
                depth(loop_n) = i
                BLOCK
                    CHARACTER(3*SIZE(depth)) :: String
                    INTEGER :: j, x

                    DO j=1, SIZE(depth)
                        x = 3*(j-1)
                        String(1+x:3+x) = base_permutations(depth(j))    !Concatenate xyz 3-characters permutations to String base on depth column
                    END DO

                    all_conc(array_location) = String

                    array_location = array_location + 1
                END BLOCK
            END DO
        END IF
    END SUBROUTINE find_all_permutations

    RECURSIVE SUBROUTINE this_permutation(all_elements, cut_sequence, no_cuts)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)
        CHARACTER(*), INTENT(INOUT) :: cut_sequence
        INTEGER, INTENT(INOUT) :: no_cuts

        INTEGER, ALLOCATABLE :: best_cuts(:)

        IF (SIZE(all_elements) == 1) THEN
            !Everything is set
        ELSE
            !Get cuts locations
            IF (no_cuts >= 5) THEN                      !5 no_cuts warranty that every coordinate has been tried
                CALL get_best_cuts(all_elements,cut_sequence(1:1), best_cuts, .TRUE.)   !Get a cut even if it is a bad one.
            ELSE
                CALL get_best_cuts(all_elements,cut_sequence(1:1), best_cuts)           !Get the best_cuts
            END IF

            cut_sequence = next_sequence(cut_sequence)  !Shift every character to the left and the first character to the end

            IF (ALLOCATED(best_cuts)) THEN              !If a good cut was found
                no_cuts = 0
                CALL make_n_cuts(all_elements, best_cuts, cut_sequence, no_cuts)    !Make those cuts
                DEALLOCATE(best_cuts)
            ELSE                                        !If no good cut was found
                no_cuts = no_cuts + 1
                CALL this_permutation(all_elements, cut_sequence, no_cuts)          !Try to find good cuts with the next coordinate
            END IF
        END IF
    END SUBROUTINE this_permutation

    SUBROUTINE make_n_cuts(elements_array, cuts_locations, cut_seq, no_cut)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: elements_array(:)
        INTEGER, INTENT(IN) :: cuts_locations(:)
        CHARACTER(*), INTENT(INOUT) :: cut_seq
        INTEGER, INTENT(INOUT) :: no_cut

        CHARACTER :: coordinate
        INTEGER :: i, j
        REAL*8 :: cut_location

        coordinate = cut_seq(LEN(cut_seq):LEN(cut_seq))         !Get the coordinate that was used to find the cuts_locations

        DO i=1, SIZE(cuts_locations)
            cut_location = calc_location(elements_array(cuts_locations(i)+1), elements_array(cuts_locations(i)), coordinate)

            !Set walls and send new boxes to get more cuts
            IF (i+1 <= SIZE(cuts_locations)) THEN                                   !More than one cut is left
                IF (i==1) THEN                                                          !Cut one when cut 2 exist
                    DO j=1, cuts_locations(i+1)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(1:cuts_locations(i)), cut_seq, no_cut)
                ELSE                                                                    !Cut i when cut i+1 exist
                    DO j=cuts_locations(i-1)+1, cuts_locations(i+1)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(cuts_locations(i-1)+1:cuts_locations(i)), cut_seq, no_cut)
                END IF
            ELSE                                                                    !One cut is left
                IF (i==1) THEN                                                          !Cut one when cut 2 not exist
                    DO j=1, SIZE(elements_array)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(1:cuts_locations(i)), cut_seq, no_cut)
                    CALL this_permutation(elements_array(cuts_locations(i)+1:SIZE(elements_array)), cut_seq, no_cut)
                ELSE                                                                    !Last cut when cut i-1 exist
                    DO j=cuts_locations(i-1)+1, SIZE(elements_array)
                        CALL elements_array(j)%set_walls(coordinate, cut_location)
                    END DO
                    CALL this_permutation(elements_array(cuts_locations(i-1)+1:cuts_locations(i)), cut_seq, no_cut)
                    CALL this_permutation(elements_array(cuts_locations(i)+1:SIZE(elements_array)), cut_seq, no_cut)
                END IF
            END IF
        END DO
    END SUBROUTINE make_n_cuts

    SUBROUTINE get_best_cuts(all_ele, coordinate, best_location, make_cut)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_ele(:)
        CHARACTER, INTENT(IN) :: coordinate
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: best_location(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: make_cut

        INTEGER :: temp_best_location(SIZE(all_ele)-1)
        INTEGER :: i, n
        REAL*8 :: distances(SIZE(all_ele)-1)
        REAL*8 :: largest_distance

        !Calculate distances between neighbor atoms
        CAll merge_sort(all_ele, coordinate)

        DO i=1, SIZE(all_ele)-1
            distances(i) = calc_distance(all_ele(i+1), all_ele(i), coordinate)
        END DO

        largest_distance = MAXVAL(distances)

        !Get largest_distance
        IF (largest_distance > -0.8 .OR. PRESENT(make_cut)) THEN        !Avoid to cut a molecule in half
            !Find for symmetry, if any
            n=0
            DO i=1, SIZE(distances)
                IF (largest_distance - distances(i) <= 0.0001) THEN     !Consideration of symmetry
                    n = n+1
                    temp_best_location(n) = i
                END IF
            END DO

            ALLOCATE(best_location(n))

            DO i=1, n
                best_location(i) = temp_best_location(i)
            END DO
        END IF
    END SUBROUTINE get_best_cuts

    PURE FUNCTION next_sequence(old_sequence)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: old_sequence
        CHARACTER(LEN(old_sequence)) :: next_sequence
        CHARACTER :: char_temp
        INTEGER :: i

        char_temp = old_sequence(1:1)                   !Save first character
        DO i=1, LEN(old_sequence)-1                     !Copy all characters one shift to the left
            next_sequence(i:i) = old_sequence(i+1:i+1)
        END DO
        next_sequence(LEN(old_sequence):LEN(old_sequence)) = char_temp  !Paste last character
    END FUNCTION next_sequence

    SUBROUTINE cuts_with_best_permutation(all_elements, best_permutation)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: all_elements(:)
        CHARACTER(*), INTENT(INOUT) :: best_permutation

        INTEGER :: i

        i=0
        CALL this_permutation(all_elements, best_permutation, i)   !Set the walls with the best cut permutation
    END SUBROUTINE cuts_with_best_permutation
END MODULE cut_molecule
