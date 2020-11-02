# irreducible_polynom
This repository includes lists of irreducible polynomials over the Galois Field 2^n.

The first line (n N) in the irreducible polynomial lists represent:

n = the degree of the polynomials in that file

N = the total number of polynomials in that file, number of lines - 1

Polynomials are represented using their coefficients. The rightmost element represents the leading coefficient:

*1 0 1 1* is equivalent to *x^3 + x^2 + 1*.

----

The lists were generated using the irreducible.py module, which finds all irreducible polynomials 
of degree n in GF(2^n) via repeated squaring.

----

This module makes heavy use of pyGF2 (https://github.com/popcornell/pyGF2) and needs NumPy to run.

Tested with Python 3.6.9 and NumPy 1.16.6
