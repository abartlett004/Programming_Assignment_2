# Programming_Assignment_2
Another programming assignment I completed.

## Task:
1. Write a function that takes a set of amino acids (say, {‘A’,’V’,’I’}) and returns the set of most efficient codons and the achieved efficiency 
For {A’,’V’,’I’} the answer would be
({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)

2. Write a function that takes a set of amino acids and:
- if they can be encoded by codons with 100% efficiency, returns the input set
- if they cannot be encoded by a codon with 100% efficiency, removes the minimal number of amino acids from the list such as the resulting list can be encoded with 100% efficiency. If there are multiple possibilities, return all of them.

I.e. from {‘A’,’V’} one would get {‘A’,‘V’} since they can be both encoded by the a number of codons with 100% efficiency (say, GYT)

From {‘A’,’I’,’V’} one would get {‘I’,’V’} and {‘A’,’V’} (both 100% efficient), but not {‘A’,’I’} since {‘A’,’I’} is only 50% efficient.

## Things to consider:
1. Write everything in Python2 or 3, make use of standard libraries such as collections and itertools. To save time, we’re providing a Python file with genetic code already typed in
2. In part (1), how do you want to handle situations when the resulting codon also encodes for stop? No single right answer here. 
3. Assume that this function will be called >100 times from a parent routine. Think about precomputing or caching things, i.e. some of the results you computed for position 1 might be reused in calculating things for position 2.

