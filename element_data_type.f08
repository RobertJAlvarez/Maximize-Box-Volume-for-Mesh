! Programmed by Robert Alvarez
! Last modified: March 15th 2022
!
! Creation of a class called elements to save an atom atribute and methods
MODULE element_data_type
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: elements, ASSIGNMENT(=), calc_distance, calc_location, makeCube

  REAL*8, PARAMETER :: PI = 4.D0*DATAN(1.D0)    !*8 is hardcoded, use SELECTED_REAL_KIND() instead

  TYPE :: elements
    PRIVATE ! Make atributes inaccessible outside of this module

    INTEGER :: element  ! Proton #
    REAL*8 :: x, y, z   ! Atom 3D location
    REAL :: radius      ! Atomic radius in radii units
    REAL*8 :: x_low_wall, x_high_wall   ! Walls, in the x-coordinate, of box containing the atom
    REAL*8 :: y_low_wall, y_high_wall   ! Walls, in the y-coordinate, of box containing the atom
    REAL*8 :: z_low_wall, z_high_wall   ! Walls, in the z-coordinate, of box containing the atom
  CONTAINS
    GENERIC :: set_molecule => set_molecule_integer, set_molecule_character ! Let the program decide for us which one to use base in the parameters
    PROCEDURE, PRIVATE, PASS :: set_molecule_integer        ! Set molecule information when an integer is pass for the number of atoms
    PROCEDURE, PRIVATE, PASS :: set_molecule_character      ! Set molecule information when an character is pass for the type of atoms
    PROCEDURE, PASS :: set_radius           ! Set the radius of the atom
    PROCEDURE, PASS :: set_walls            ! Update the atom walls
    PROCEDURE, PASS :: reset_walls          ! Set the walls to its maximum length
    PROCEDURE, PASS :: less_equal_than      ! Compare two atoms base in a coordinate
    PROCEDURE, PASS :: calculate_volume     ! Calculate volume made by the walls
    PROCEDURE, PASS :: print_info           ! Print atom information: element, location, and walls
    PROCEDURE, PASS :: makeCube             ! Convert the box into a cube
  END TYPE elements

  ! Extend the meaining of '=' to work with the elements derived data type
  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE copy_array
  END INTERFACE

  CONTAINS

  ! Set element type and x, y, and z position. Element type is pass as an integer representing the number of protons
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

  ! Set element type and x, y, and z position. Element type is pass as a character representing the abbreviation use on the periodic table
  SUBROUTINE set_molecule_character(this, ele_pass, x_pass, y_pass, z_pass)
    IMPLICIT NONE
    CLASS(elements), INTENT(INOUT) :: this
    CHARACTER(*), INTENT(IN) :: ele_pass
    REAL*8, INTENT(IN) :: x_pass, y_pass, z_pass

    this%element = get_element(TRIM(ele_pass))  ! Find the number of protons that the character represent base on the periodic table
    this%x = x_pass
    this%y = y_pass
    this%z = z_pass

    ! Check is character does represent an element in the periodic table
    IF (this%element == -1) WRITE(*,*) 'Element read from file does not exist'
  END SUBROUTINE set_molecule_character

  ! Use element type to set its radii
  SUBROUTINE set_radius(this)
    IMPLICIT NONE
    CLASS(elements), INTENT(INOUT) :: this
    REAL :: all_radius(85)
    INTEGER :: i

    all_radius = [0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56, 0.48, 0.42, 0.38, 1.90, 1.45, 1.18, 1.11, 0.98, &
    0.88, 0.79, 0.71, 2.43, 1.94, 1.84, 1.76, 1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25, 1.14, &
    1.03, 0.94, 0.88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.90, 1.83, 1.78, 1.73, 1.69, 1.65, 1.55, 1.56, 1.45, 1.33, &
    1.23, 1.15, 1.08, 2.98, 2.53, 2.26, 2.10, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33, 2.25, 2.28, 2.26, 2.26, 2.22, &
    2.22, 2.17, 2.08, 2.00, 1.93, 1.88, 1.85, 1.80, 1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, 1.27]

    this%radius = all_radius(this%element)
  END SUBROUTINE set_radius

  ! Update wall values base in the coordinate and location of the cut
  SUBROUTINE set_walls(this, coordinate, location)
    IMPLICIT NONE
    CLASS(elements), INTENT(INOUT) :: this  ! Element in consideration
    CHARACTER, INTENT(IN) :: coordinate     ! Coordinate to update walls
    REAL*8, INTENT(IN) :: location          ! Location of the cut

    SELECT CASE(coordinate)
      CASE ('x')  ! Cut was perform in the x-coordinate
        IF (this%x > location) THEN                     ! If cut was between the left wall and the atom location
          this%x_low_wall = location
        ELSE                                            ! If cut was between the right wall and the atom location
          this%x_high_wall = location
        END IF
      CASE('y')   ! Cut was perform in the y-coordinate
        IF (this%y > location) THEN                     ! If cut was between the left wall and the atom location
          this%y_low_wall = location
        ELSE                                            ! If cut was between the right wall and the atom location
          this%y_high_wall = location
        END IF
      CASE('z')   ! Cut was perform in the z-coordinate
        IF (this%z > location) THEN                     ! If cut was between the left wall and the atom location
          this%z_low_wall = location
        ELSE                                            ! If cut was between the right wall and the atom location
          this%z_high_wall = location
        END IF
    END SELECT
  END SUBROUTINE set_walls

  ! Set walls to the maximum range of possible interaction between the elements and the surrounding
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

  ! Check if the coordinate possition of 'this' is greater or not compare to the coordinate possition of 'another_element'
  PURE LOGICAL FUNCTION less_equal_than(this, another_element, coordinate)
    IMPLICIT NONE
    CLASS(elements), INTENT(IN) :: this, another_element    ! Element pass, and element to campare with
    CHARACTER, INTENT(IN) :: coordinate                     ! Coordinate axis of reference

    less_equal_than = .FALSE.

    SELECT CASE(coordinate)
      CASE('x')   ! When sorting with respect to x
        IF (this%x <= another_element%x) less_equal_than = .TRUE.
      CASE('y')   ! with respect to y
        IF (this%y <= another_element%y) less_equal_than = .TRUE.
      CASE('z')   ! with respect to z
        IF (this%z <= another_element%z) less_equal_than = .TRUE.
      CASE DEFAULT    ! with respect to the element number
        IF (this%element <= another_element%element) less_equal_than = .TRUE.
    END SELECT
  END FUNCTION less_equal_than

  ! Calculate the volume inside of the box generated by the walls of the element
  PURE REAL*8 FUNCTION calculate_volume(this)
    IMPLICIT NONE
    CLASS(elements), INTENT(IN) :: this ! Element in consideration

    calculate_volume = (this%x_high_wall - this%x_low_wall) * (this%y_high_wall - this%y_low_wall)
    calculate_volume = calculate_volume *  (this%z_high_wall - this%z_low_wall)
  END FUNCTION calculate_volume

    ! Calculate the distance between the surfaces of two atoms
    PURE REAL*8 FUNCTION calc_distance(high_ele, low_ele, coordinate)
        IMPLICIT NONE
        CLASS(elements), INTENT(IN) :: high_ele, low_ele    ! two elements
        CHARACTER, INTENT(IN) :: coordinate   ! Coordinate of consideration

        ! Calculate distance from the center of the atoms
        SELECT CASE(coordinate)
            CASE('x')     ! in the x-coordinate
                calc_distance = high_ele%x - low_ele%x
            CASE('y')     ! in the y-coordinate
                calc_distance = high_ele%y - low_ele%y
            CASE('z')     ! in the z-coordinate
                calc_distance = high_ele%z - low_ele%z
        END SELECT
        calc_distance = calc_distance - high_ele%radius - low_ele%radius  ! Subtract the radius to get the distance from the surfaces of the atoms
    END FUNCTION calc_distance

    ! 
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
