This code can be used to search for signed diffrential characteristics for RIPEMD-160.

When you use the code, you need to edit three files in data/left or data/right: msgFile, diffFile, parFile.

# Edit msgFile

In msgFile, you need to input the signed difference of the message words.

For example, in the data/left/msgFile, we have

        00 =========n==================n===

        01 ================================

        02 ================================

        03 ================================

        04 ================================

        05 ================================

        06 ===u============u===============

        07 ================================

        08 ================================

        09 ================n============n==

        10 ================================

        11 ================================

        12 ================================

        13 ================================

        14 ================================

        15 ================================


which means we have differences on m0, m6 and m9 and the differences are represented by "u", "n" and "=".
Please make sure there is a blank after the number "00", ..., "15".


# Edit diffFile

In diffFile, we need to specify the targeted differential characteristic.

For example, in the data/left/diffFile, we have

        -5 ================================

        -4 ================================

        -3 ================================

        -2 ================================

        -1 ================================

        00 ????????????????????????????????

        01 ????????????????????????????????

        02 ????????????????????????????????

        03 ????????????????????????????????

        04 ????????????????????????????????

        05 ================================

        06 ================================

        07 ================================

        08 ================================

        09 ================================

This means that we want to find a solution of step 0 - step 4, given the differential characteristic for Step (-5) - (-1) and Step 5 - Step 9.


# Edit parFile
Finally, we need to specify the parameters for the search, which is recorded in parFile.

For example, in data/left/parFile, we have

        15

        0

        4

        0

        -1

        0

        0 0 0 0 0 0 0 0 0 0

        1 1 1 1 0 1 1 1 1 1

        0 0 0 0 0 0 0 0 0 0

1. The 1st number represents the total number of steps. As shown in data/left/diffFile, there are in total 15 states. So, we have 15 steps.



2. The 2nd number represents the starting position and it should be adjusted according to diffFile. 

      * In data/left/diffFile, if the first position is x, the starting position in parFile should be x+5.
        
      * For example, in data/left/diffFile, x = -5. So, the starting position is 0 in parFile.



3. The 3rd number represents the ending position of the unknown state diff.
 
      * In data/left/diffFile, the last unknown state is at step 04. So, we set 4 for the 3rd number.

      * Note that if you set 5, ..., 9 for the 3rd number, it does not matter but it may affect the performance. 

      * The third number is only used to control the strategy of the expansion model.

      * In addition, if the diff is fully specified (as in data/left/text_diffFile), the 3rd number can be any number.


4. The 4th number and 5th number are used to control which part of the differential characteristic should be optimized, i.e. of low hamming weight.

      * The 4th number means the starting position of the to-be-optimized part.

      * The 5th number means the length of the to-be-optimized part.

      * In this example, the 4th/5th numbers are 0 and -1, reps. This means that we do not optimize the HW of the differential characteristics.

      * In this example, if we set the 4th/5th numbers to 0 and 5, reps, it means that we need to optimize the HW of the differential characteristics at position 0, 1, 2, 3, 4.

5. The 6th number represents the branch. If it is 0/1, it means it is the left/right branch.

      * In this example, we are considering the left branch, so, it is 0.

6. Finally, we need to configure the searching parameters: isC, isF, isV, all the which should be of size y - 5, where y is the 1st number in parFile.

      * In this example, the first number is 15, so the size of isC, isF, isV are all 15 - 5 = 10.

      * Note that in the paper, isK is also a parameter, but as can be seen from the code, we control isK according to the 3rd number.

