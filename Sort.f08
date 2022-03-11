MODULE sort
    USE element_data_type, ONLY: elements, ASSIGNMENT(=)
    IMPLICIT NONE

    PRIVATE :: Merge

	CONTAINS

    RECURSIVE SUBROUTINE merge_sort(a, coordinate)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: a(:)
        CHARACTER, INTENT(IN) :: coordinate
        INTEGER :: low, high, mid

        low = LBOUND(a,1)
        high = UBOUND(a,1)

        IF (low < high) THEN
            mid = low + (high - low)/2
            CALL merge_sort(a(low:mid), coordinate)
            CALL merge_sort(a(mid+1:high), coordinate)
            a(low:high) = Merge(a(low:mid), a(mid+1:high), coordinate)
        END IF
    END SUBROUTINE merge_sort

    FUNCTION Merge(a, b, coordinate)
        IMPLICIT NONE
        CLASS(elements), DIMENSION(:), INTENT(INOUT) :: a, b
        TYPE(elements), DIMENSION(SIZE(a)+SIZE(b)) :: Merge
        CHARACTER, INTENT(IN) :: coordinate
        INTEGER a_ptr, a_high
        INTEGER b_ptr, b_high
        INTEGER c_ptr

        a_ptr = LBOUND(a,1)
        a_high = UBOUND(a,1)
        b_ptr = LBOUND(b,1)
        b_high = UBOUND(b,1)
        c_ptr = 1

        DO WHILE (a_ptr <= a_high .AND. b_ptr <= b_high)
            IF ( a(a_ptr)%less_equal_than(b(b_ptr), coordinate) ) THEN
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
