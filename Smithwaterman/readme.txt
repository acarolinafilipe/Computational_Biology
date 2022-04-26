The main function is smithwaterman and is the one to be used. An auxilary function called find_path is used to trace back the arrow matrix and return the local alignemnts to the main function.
The Bioinformatics toolbox was used to facilitate the BLOSUM50 scoring.

3 inputs should be provided:
s1 and s2 - The sequences you wish to align. Write them in char format, ex: 'TECTEA'
d - The module of the linear gap penalty.

Two outputs are obtained:
score - The score of the optimal local alignment(s).
list_opt - A matrix in which each column is an optimal local alignment. The first row corresponds to the first sequence while the second row corresponds to the second sequence.


EXAMPLE:

Introduce the following command:

	[a,b] = smithwaterman('TECTEA','CCTEC',5)

You should get the following output:

a =

    24


b = 

  2×2 string array

    "TEC"    "CTE"
    "TEC"    "CTE"
 