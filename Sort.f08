! Programmed by Robert Alvarez
! Last modified: March 15th, 2022
!
! Modify merge sort to work with the data type defines as element_data_type
MODULE sort
  USE element_data_type, ONLY: elements, ASSIGNMENT(=)	! Use class element_data_type and extend '=' meaning to work with array filled with element_data_type objects
  IMPLICIT NONE

  PRIVATE :: Merge

	CONTAINS

  RECURSIVE SUBROUTINE merge_sort(a, coordinate)
    IMPLICIT NONE
    CLASS(elements), INTENT(INOUT) :: a(:)
    CHARACTER, INTENT(IN) :: coordinate
    INTEGER :: low, high, mid

    low = LBOUND(a,1)		! Number of the first index of array a
    high = UBOUND(a,1)	! Number of the last index of array a

    IF (low < high) THEN
      mid = low + (high - low)/2
      CALL merge_sort(a(low:mid), coordinate)
      CALL merge_sort(a(mid+1:high), coordinate)
      a(low:high) = Merge(a(low:mid), a(mid+1:high), coordinate)
    END IF
  END SUBROUTINE merge_sort

  FUNCTION Merge(a, b, coordinate)
    IMPLICIT NONE
    CLASS(elements), DIMENSION(:), INTENT(INOUT) :: a, b	! Pass pair of arrays
    TYPE(elements), DIMENSION(SIZE(a)+SIZE(b)) :: Merge		! Sorted combination of the elements in arrays a and b
    CHARACTER, INTENT(IN) :: coordinate										! Coordinate use as reference to sort elements
    INTEGER a_ptr, a_high
    INTEGER b_ptr, b_high
    INTEGER c_ptr

    a_ptr = LBOUND(a,1)		! Number of the first index of array a
    a_high = UBOUND(a,1)	! Number of the last index of array a
    b_ptr = LBOUND(b,1)		! Number of the first index of array b
    b_high = UBOUND(b,1)	! Number of the last index of array b
    c_ptr = 1

    DO WHILE (a_ptr <= a_high .AND. b_ptr <= b_high)
      IF ( a(a_ptr)%less_equal_than(b(b_ptr), coordinate) ) THEN	! less_equal_than make a comparition base in the coordinate possition
        Merge(c_ptr) = a(a_ptr)
        a_ptr = a_ptr + 1
      ELSE
        Merge(c_ptr) = b(b_ptr)
        b_ptr = b_ptr + 1
      END IF
      c_ptr = c_ptr + 1
    END DO

    IF (a_ptr > a_high) THEN
      Merge(c_ptr:) = b(b_ptr:b_high)
    ELSE
      Merge(c_ptr:) = a(a_ptr:a_high)
    END IF
  END FUNCTION Merge
END MODULE sort
