MODULE element_data_type
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: elements, ASSIGNMENT(=), calc_distance, calc_location, makeCube

    REAL*8, PARAMETER :: PI = 4.D0*DATAN(1.D0)    !*8 is hardcoded, use SELECTED_REAL_KIND() instead

    TYPE :: elements
        PRIVATE

        INTEGER :: element  !Proton #
        REAL*8 :: x, y, z   !Atom 3D location
        REAL :: radius      !Atomic radius
        REAL*8 :: x_low_wall, x_high_wall
        REAL*8 :: y_low_wall, y_high_wall
        REAL*8 :: z_low_wall, z_high_wall
    CONTAINS
        GENERIC :: set_molecule => set_molecule_integer, set_molecule_character
        PROCEDURE, PRIVATE, PASS :: set_molecule_integer
        PROCEDURE, PRIVATE, PASS :: set_molecule_character
        PROCEDURE, PASS :: set_radius
        PROCEDURE, PASS :: set_walls
        PROCEDURE, PASS :: reset_walls
        PROCEDURE, PASS :: less_equal_than
        PROCEDURE, PASS :: calculate_volume
        PROCEDURE, PASS :: print_info
        PROCEDURE, PASS :: makeCube
    END TYPE elements

    INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE copy_array
    END INTERFACE

    CONTAINS

    SUBROUTINE set_molecule_integer(this, ele_pass, x_pass, y_pass, z_pass)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: ele_pass
        REAL*8, INTENT(IN) :: x_pass, y_pass, z_pass

        this%element = ele_pass
        this%x = x_pass
        this%y = y_pass
        this%z = z_pass
    END SUBROUTINE set_molecule_integer

    SUBROUTINE set_molecule_character(this, ele_pass, x_pass, y_pass, z_pass)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: this
        CHARACTER(*), INTENT(IN) :: ele_pass
        REAL*8, INTENT(IN) :: x_pass, y_pass, z_pass

        this%element = get_element(TRIM(ele_pass))
        this%x = x_pass
        this%y = y_pass
        this%z = z_pass

        IF (this%element == -1) THEN
            WRITE(*,*) 'Element read from file does not exist'
        END IF
    END SUBROUTINE set_molecule_character

    SUBROUTINE set_radius(this)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: this
        REAL :: all_radius(85)
        INTEGER :: i
! 0.529
!            all_radius = [1.000, 1.400, 1.520, 1.113, 0.795, 1.700, 1.500, 1.400, 1.400, 1.500, 1.858, 1.599, 1.432 ,1.176, &
!            1.105, 1.800, 1.800, 1.800, 2.272, 1.974, 1.606, 1.448, 1.311, 1.249, 1.367, 1.241, 1.253, 1.246, 1.278, 1.335, &
!            1.221, 1.225, 1.245, 1.160, 2.000, 1.900, 2.475, 2.151, 1.776, 1.590, 1.429, 1.363, 1.352, 1.325, 1.345, 1.376, &
!            1.445, 1.489, 1.626, 1.405, 1.450, 1.432, 2.200, 2.100, 2.655, 2.174, 1.870, 1.825, 1.820, 1.814, 1.630 ,1.620, &
!            1.995, 1.787, 1.763, 1.752, 1.743, 1.734, 1.724, 1.940, 1.718, 1.564, 1.430, 1.370, 1.371, 1.338, 1.357, 1.371, &
!            1.442, 1.503, 1.700, 1.750, 1.545, 1.673, 1.450, 2.300]

        all_radius = [0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56, 0.48, 0.42, 0.38, 1.90, 1.45, 1.18, 1.11, 0.98, &
        0.88, 0.79, 0.71, 2.43, 1.94, 1.84, 1.76, 1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25, 1.14, &
        1.03, 0.94, 0.88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.90, 1.83, 1.78, 1.73, 1.69, 1.65, 1.55, 1.56, 1.45, 1.33, &
        1.23, 1.15, 1.08, 2.98, 2.53, 2.26, 2.10, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33, 2.25, 2.28, 2.26, 2.26, 2.22, &
        2.22, 2.17, 2.08, 2.00, 1.93, 1.88, 1.85, 1.80, 1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, 1.27]

        this%radius = all_radius(this%element)
    END SUBROUTINE set_radius

    SUBROUTINE set_walls(this, coordinate, location)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: this
        CHARACTER, INTENT(IN) :: coordinate
        REAL*8, INTENT(IN) :: location

        SELECT CASE(coordinate)
            CASE ('x')
                IF (this%x > location) THEN
                    IF (this%x_low_wall < location) THEN
                        this%x_low_wall = location
                    END IF
                ELSE
                    IF (this%x_high_wall > location) THEN
                        this%x_high_wall = location
                    END IF
                END IF
            CASE('y')
                IF (this%y > location) THEN
                    IF (this%y_low_wall < location) THEN
                        this%y_low_wall = location
                    END IF
                ELSE
                    IF (this%y_high_wall > location) THEN
                        this%y_high_wall = location
                    END IF
                END IF
            CASE('z')
                IF (this%z > location) THEN
                    IF (this%z_low_wall < location) THEN
                        this%z_low_wall = location
                    END IF
                ELSE
                    IF (this%z_high_wall > location) THEN
                        this%z_high_wall = location
                    END IF
                END IF
        END SELECT
    END SUBROUTINE set_walls

    SUBROUTINE reset_walls(ele_pass)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: ele_pass

        ele_pass%x_low_wall = -50.0
        ele_pass%x_high_wall = 50.0
        ele_pass%y_low_wall = -50.0
        ele_pass%y_high_wall = 50.0
        ele_pass%z_low_wall = -50.0
        ele_pass%z_high_wall = 50.0
    END SUBROUTINE reset_walls

    PURE LOGICAL FUNCTION less_equal_than(this, another_element, coordinate)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: this, another_element
        CHARACTER, INTENT(IN) :: coordinate

        less_equal_than = .FALSE.

        SELECT CASE(coordinate)
            CASE('x')                                    !When sorting with respect to x
                IF (this%x <= another_element%x) THEN
                    less_equal_than = .TRUE.
                END IF
            CASE('y')                                    !with respect to y
                IF (this%y <= another_element%y) THEN
                    less_equal_than = .TRUE.
                END IF
            CASE('z')                                    !with respect to z
                IF (this%z <= another_element%z) THEN
                    less_equal_than = .TRUE.
                END IF
            CASE DEFAULT
                IF (this%element <= another_element%element) THEN
                    less_equal_than = .TRUE.
                END IF
        END SELECT
    END FUNCTION less_equal_than

    PURE REAL*8 FUNCTION calculate_volume(this)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: this
        REAL*8 :: volume

        volume = 1.0D0 * (this%x_high_wall - this%x_low_wall)
        volume = volume * (this%y_high_wall - this%y_low_wall)
        volume = volume * (this%z_high_wall - this%z_low_wall)

        calculate_volume = ABS(volume)! - (4.0 * PI * this%radius**3)/3.0
    END FUNCTION calculate_volume

    PURE REAL*8 FUNCTION calc_distance(high_ele, low_ele, coordinate)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: high_ele, low_ele
        CHARACTER, INTENT(IN) :: coordinate

        SELECT CASE(coordinate)
            CASE('x')
                calc_distance = high_ele%x - low_ele%x
            CASE('y')
                calc_distance = high_ele%y - low_ele%y
            CASE('z')
                calc_distance = high_ele%z - low_ele%z
        END SELECT
        calc_distance = calc_distance - high_ele%radius - low_ele%radius
    END FUNCTION calc_distance

    PURE REAL*8 FUNCTION calc_location(high_ele, low_ele, coordinate)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: high_ele, low_ele
        CHARACTER, INTENT(IN) :: coordinate

        calc_location = calc_distance(high_ele, low_ele, coordinate) / 2.0D0

        SELECT CASE(coordinate)
            CASE('x')
                calc_location = calc_location + low_ele%x + low_ele%radius
            CASE('y')
                calc_location = calc_location + low_ele%y + low_ele%radius
            CASE('z')
                calc_location = calc_location + low_ele%z + low_ele%radius
        END SELECT
    END FUNCTION calc_location

    SUBROUTINE print_info(this, fileNum)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: fileNum

        IF (PRESENT(fileNum)) THEN
            WRITE(fileNum,*) this%x, this%x_low_wall, this%x_high_wall
            WRITE(fileNum,*) this%y, this%y_low_wall, this%y_high_wall
            WRITE(fileNum,*) this%z, this%z_low_wall, this%z_high_wall
            WRITE(fileNum,*) ' '
        ELSE
            WRITE(*,*) this%x, this%x_low_wall, this%x_high_wall
            WRITE(*,*) this%y, this%y_low_wall, this%y_high_wall
            WRITE(*,*) this%z, this%z_low_wall, this%z_high_wall
            WRITE(*,*) ' '
        END IF
    END SUBROUTINE print_info

    SUBROUTINE copy_array(array_result, array_pass)
        IMPLICIT NONE
        TYPE(elements), INTENT(OUT) :: array_result(:)
        TYPE(elements), INTENT(IN) :: array_pass(:)
        INTEGER :: k

        DO k=1, SIZE(array_pass)
            array_result(k) = array_pass(k)
        END DO
    END SUBROUTINE copy_array

    PURE INTEGER FUNCTION get_element(find_ele)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: find_ele

        CHARACTER(5) :: all_elements(118)
        INTEGER :: mid, low, high

        all_elements = ['Ac 89', 'Ag 47', 'Al 13', 'Am 95', 'Ar 18', 'As 33', 'At 85', 'Au 79', 'B   5', 'Ba 56', &
        'Be  4', 'Bh107', 'Bi 83', 'Bk 97', 'Br 35', 'C   6', 'Ca 20', 'Cd 48', 'Ce 58', 'Cf 98', 'Cl 17', 'Cm 96', &
        'Cn112', 'Co 27', 'Cr 24', 'Cs 55', 'Cu 29', 'Db105', 'Ds110', 'Dy 66', 'Ed 99', 'Er 68', 'Eu 63', 'F   9', &
        'Fe 26', 'Fi114', 'Fm100', 'Fr 87', 'Ga 31', 'Gd 64', 'Ge 32', 'H   1', 'He  2', 'Hf 72', 'Hg 80', 'Ho 67', &
        'Hs108', 'I  53', 'In 49', 'Ir 77', 'K  19', 'Kr 36', 'La 57', 'Li  3', 'Lr103', 'Lu 71', 'Lv116', 'Mc115', &
        'Md101', 'Mg 12', 'Mn 25', 'Mo 42', 'Mt109', 'N   7', 'Na 11', 'Nb 41', 'Nd 60', 'Ne 10', 'Nh113', 'Ni 28', &
        'No102', 'Np 93', 'O   8', 'Og118', 'Os 76', 'P  15', 'Pa 91', 'Pb 82', 'Pd 46', 'Pm 61', 'Po 84', 'Pr 59', &
        'Pt 78', 'Pu 94', 'Ra 88', 'Rb 37', 'Re 75', 'Rf104', 'Rg111', 'Rh 45', 'Rn 86', 'Ru 44', 'S  16', 'Sb 51', &
        'Sc 21', 'Se 34', 'Sg106', 'Si 14', 'Sm 62', 'Sn 50', 'Sr 38', 'Ta 73', 'Tb 65', 'Tc 43', 'Te 52', 'Th 90', &
        'Ti 22', 'Tl 81', 'Tm 69', 'Ts117', 'U  92', 'V  23', 'W  74', 'Xe 54', 'Y  39', 'Yb 70', 'Zn 30', 'Zr 40']

        low = 1
        high = SIZE(all_elements)
        get_element = -1

        !Binary search
        DO WHILE (high >= low)
            mid = (high + low) / 2
            IF (all_elements(mid)(1:2) < find_ele) THEN
                low = mid + 1
            ELSE IF (all_elements(mid)(1:2) > find_ele) THEN
                high = mid - 1
            ELSE
                READ(all_elements(mid)(3:5), *) get_element
                EXIT
            END IF
        END DO
    END fUNCTION get_element

    SUBROUTINE makeCube(this, arr, atP)
        IMPLICIT NONE
        CLASS(elements), INTENT(INOUT) :: this, arr(:)
        INTEGER, INTENT(INOUT) :: atP   !atP stands for atPosition

        REAL*8 :: Lx, Rx, Ly, Ry, Lz, Rz, shortestD, cutAt

        Rx = this%x_high_wall - this%x
        Lx = this%x - this%x_low_wall
        Ry = this%y_high_wall - this%y
        Ly = this%y - this%y_low_wall
        Rz = this%z_high_wall - this%z
        Lz = this%z - this%z_low_wall

        shortestD = MIN(Rx, Lx, Ry, Ly, Rz, Lz)

        IF (Lx /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%x - shortestD

            !Set new walls
            this%x_low_wall = cutAt
            arr(atP)%x_high_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF

        IF (Rx /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%x + shortestD

            !Set new walls
            this%x_high_wall = cutAt
            arr(atP)%x_low_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF

        IF (Ly /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%y - shortestD

            !Set new walls
            this%y_low_wall = cutAt
            arr(atP)%y_high_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF

        IF (Ry /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%y + shortestD

            !Set new walls
            this%y_high_wall = cutAt
            arr(atP)%y_low_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF

        IF (Lz /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%z - shortestD

            !Set new walls
            this%z_low_wall = cutAt
            arr(atP)%z_high_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF

        IF (Rz /= shortestD) THEN
            !Set variables
            CALL copyWalls(this, arr(atP))
            cutAt = this%z + shortestD

            !Set new walls
            this%z_high_wall = cutAt
            arr(atP)%z_low_wall = cutAt

            !Update position to store next empty box
            atP = atP + 1
        END IF
    END SUBROUTINE makeCube

    SUBROUTINE copyWalls (org, copy)
        CLASS(elements), INTENT(IN) :: org
        CLASS(elements), INTENT(OUT) :: copy

        copy%element = -1
        copy%x_low_wall = org%x_low_wall
        copy%x_high_wall = org%x_high_wall
        copy%y_low_wall = org%y_low_wall
        copy%y_high_wall = org%y_high_wall
        copy%z_low_wall = org%z_low_wall
        copy%z_high_wall = org%z_high_wall
    END SUBROUTINE copyWalls
END MODULE element_data_type
