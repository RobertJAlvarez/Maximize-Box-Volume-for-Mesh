! Programmed by Robert Alvarez
! Last modified: March 15th, 2022
!
! Modify merge sort to work with the data type defines as element_data_type
MODULE sort
   USE element_data_type, ONLY: elements, ASSIGNMENT(=)  ! Use class element_data_type and extend '=' meaning to work with array filled with element_data_type objects
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: merge_sort

CONTAINS

   RECURSIVE SUBROUTINE merge_sort(a, coordinate)
      IMPLICIT NONE
      CLASS(elements), INTENT(INOUT) :: a(:)
      CHARACTER, INTENT(IN) :: coordinate
      INTEGER :: high, mid

      high = UBOUND(a, 1)  ! Number of the last index of array a

      IF (1 < high) THEN
         mid = (high + 1)/2
         CALL merge_sort(a(:mid), coordinate)
         CALL merge_sort(a(mid + 1:), coordinate)
         a(:) = my_merge(a(:mid), a(mid + 1:), coordinate)
      END IF
   END SUBROUTINE merge_sort

   FUNCTION my_merge(a, b, coordinate)
      IMPLICIT NONE
      CLASS(elements), DIMENSION(:), INTENT(IN) :: a, b ! Pass pair of arrays
      CHARACTER, INTENT(IN) :: coordinate               ! Coordinate use as reference to sort elements

      TYPE(elements) :: my_merge(SIZE(a) + SIZE(b))     ! Sorted combination of the elements in arrays a and b
      INTEGER :: ai, a_high, bi, b_high, ci

      ai = 1
      a_high = SIZE(a)  ! Number of the last index of array a
      bi = 1
      b_high = SIZE(b)  ! Number of the last index of array b
      ci = 1

      DO WHILE (ai <= a_high .AND. bi <= b_high)
         IF (a(ai)%less_equal_than(b(bi), coordinate)) THEN  ! less_equal_than make a comparison base in the coordinate position
            my_merge(ci) = a(ai)
            ai = ai + 1
         ELSE
            my_merge(ci) = b(bi)
            bi = bi + 1
         END IF
         ci = ci + 1
      END DO

      IF (ai > a_high) THEN
         my_merge(ci:) = b(bi:)
      ELSE
         my_merge(ci:) = a(ai:)
      END IF
   END FUNCTION my_merge
END MODULE sort
