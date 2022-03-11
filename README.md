# Maximize-Box-Volume-for-Mesh
A molecule is given and we return the boxes made to maximize the smaller box in the molecule.

Pseuducode:

- Read molecules information and set box walls.

- Find all permutations possible with n concatenations between xyz, xzy, yxz, yzx, zxy, and zyx.

- Use every permutation, one by one, to cut the molecule.

  - IF a box have one atom, THEN go back (you are done)

  - ELSE

    - Sort the atoms with respect to the first character of the permutation.

    - Calculate all distances between atoms radius and find the largest distance (x2 – x1 – radius2 – radius1)

      - Also look for symmetry, more than one gap can be the biggest gap.

    - On the permutation just used, shift all the characters one to the left and move the first one to the end.

    - IF no good cut was found THEN go back to line 4

    - ELSE, handle cuts scenarios

      - Send surround atoms to modify its walls.
  
      - Send group of atoms to line 4 to get more cuts.
