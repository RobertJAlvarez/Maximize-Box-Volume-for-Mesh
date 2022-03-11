# Maximize-Box-Volume-for-Mesh
A molecule is given and we return the boxes made to maximize the smaller box in the molecule.
Pseuducode:
1. Read molecules information and set box walls.
2. Find all permutations possible with n concatenations between xyz, xzy, yxz, yzx, zxy, and zyx.
3. Use every permutation, one by one, to cut the molecule.
  a) IF a box have one atom, THEN go back (you are done)
  b) ELSE
    I) Sort the atoms with respect to the first character of the permutation.
    II) Calculate all distances between atoms radius and find the largest distance (x2 – x1 – radius2 – radius1)
      -) Also look for symmetry, more than one gap can be the biggest gap.
    III) On the permutation just used, shift all the characters one to the left and move the first one to the end.
    IV) IF no good cut was found THEN go back to line 4
    V) ELSE
      -) Handle cuts scenarios
        A) Send surround atoms to modify its walls.
        B) Send group of atoms to line 4 to get more cuts.
