Each function accepts monomials and polynomials both in traditional
and in parsed form as input. 
Example structures:

Monomial in traditional form:    5*y
Monomial in parsed form:         m(5, 1, [v(1, y)])
Polynomial in traditional form:  3*x^2+3
Polynomial in parsed form:       poly([m(3, 0, []), m(3, 2, [v(2, x)])])

The functions as-monomial and as-polynomial required a monomial and a polynomial
in traditional form respectively. They do not support objects in parsed form. 
The input is parsed, normalised and ordered lexicographically. 

The functions polyplus, polyminus, polytimes accept any kind of input in a pairwise
combination, but only if both objects adhere to the same structure (both in traditional
form or both in parsed form). Unsupported combinations produce FALSE as result. 
Null values are interpreted as 0 (zero).

It is assumed that parsed monomial and polynomial given as input are correct, eventually
not normalised and unordered. These operations are executed automatically by each function.

The value 0 (zero) is represented stand-alone in the following form:
m(0, 0, [])
If the value 0 (zero) is included within a Polynomial, it is ignored. 

Any variable raised to the power 0 (zero) are automatically interpreted as 1.

Traditional form ---> Parsed form :
0   --->  poly([m(0, 0, [])])
x+0  --->  poly([m(1, 1, [v(1, x)])])
3*x-3*x+2*y  --->  poly([m(2, 1, [v(1, y)])])
3*x*0+2  --->  poly([m(2, 0, [])])
x^0  --->  poly([m(1, 0, [])])
3*x^0+2  --->  poly([m(5, 0, [])])
x^2*x^(-2)*a  --->  poly([m(1, 1, [v(1, a)])]
