
% coefficients(Poly, Coefficients)
% Returns TRUE if: Coefficients is the list of the coefficients in Poly 
% 				   ordered according to the position of the respective monomial in Poly
% 				   after the sorting and normalisation process.

coefficients(poly([]), []) :- !.
coefficients(poly(Mons), Rs) :-
	is_polynomial(poly(Mons)),
	merge_sort(Mons, MonsSorted),
	normalize_poly(MonsSorted, MonsN),
	MonsN = [m(C,_,_)|Ys],
        coefficients(poly(Ys), Ms),
	append([C], Ms, Rsm),
	Rsm == [0], !,
	Rs = [].
coefficients(poly(Mons), Rs) :-
	is_polynomial(poly(Mons)),
	merge_sort(Mons, MonsSorted),
	normalize_poly(MonsSorted, MonsN),
	MonsN = [m(C,_,_)|Ys],
        coefficients(poly(Ys), Ms),
	append([C], Ms, Rs), !.
coefficients(Mon, []) :-
	Mon = m(C, _Td, _Vps),
	is_monomial(Mon),
	C == 0, !.
coefficients(Mon, [C]) :-
	Mon = m(C, _, _),
	is_monomial(Mon).
coefficients(P, R) :-
	is_trad_poly(P),
	as_polynomial(P, Pt),
	coefficients(Pt, R).



% variables(Poly, Variables).
% Returns TRUE if: Variables is the ordered list of unique symbols in Poly

variables(poly([m(_,_,[])]), []) :- !.
variables(poly([]), []) :- !.
variables(poly(Mons), Zs) :-
	is_polynomial(poly(Mons)),
	merge_sort(Mons, MonsSorted),
	normalize_poly(MonsSorted, MonsN),
	MonsN = [m(C, Td, Vps)|Ys],
	Vps == [], !,
	variables(poly([m(C, Td, Vps)]), Ms),
	variables(poly(Ys), Ns),
	append(Ms, Ns, MNs),
	sort(MNs, Zs).
variables(poly(Mons), Zs) :-
	is_polynomial(poly(Mons)),
	merge_sort(Mons, MonsSorted),
	normalize_poly(MonsSorted, MonsN),
	MonsN = [m(C, Td, Vps)|Ys],
	Vps \= [], !,
	Vps = [v(_, X)|Ws],
	variables(poly([m(C, Td, Ws)]), Ms),
	variables(poly(Ys), Ns),
	append(Ms, Ns, MNs),
	append([X], MNs, MNXs),
	sort(MNXs, Zs).
variables(Mon, []) :-
	Mon = m(_C, _Td, Vps),
	is_monomial(Mon),
	merge_sort(Vps, VpsSorted),
	normalize_mons(VpsSorted, VpsN),
	VpsN == [].
variables(Mon, Zs) :-
	is_monomial(Mon),
	Mon = m(C, Td, Vps),
	merge_sort(Vps, VpsSorted),
	normalize_mons(VpsSorted, VpsN),
	VpsN = [v(_, X)|Ws],
	variables(m(C, Td, Ws), Ns),
	append([X], Ns, XNs),
	sort(XNs, Zs).
variables(M1, R) :-
	is_trad_poly(M1),
	as_polynomial(M1, Mp1),
	variables(Mp1, R).


% monomials(Poly, Monomials).
% Returns TRUE if: Monomials is the ordered list of the monomials in Poly

monomials(Mon, R) :-
	is_monomial(Mon),
	monomials(poly([Mon]), R).
monomials(poly([]), []).
monomials(poly(Mons), R) :-
	is_polynomial(poly(Mons)),
	merge_sort(Mons, Mt),
	normalize_poly(Mt, R).
monomials(P, R) :-
	is_trad_poly(P),
	as_polynomial(P, Pt),
	monomials(Pt, R).



% maxdegree(Poly, Degree).
% Returns TRUE if: Degree is the highest degree amongst the monomials in Poly

maxdegree([],0) :- !.
maxdegree(Mon, Degree) :-
	is_monomial(Mon), !,
	maxdegree(poly([Mon]), Degree).
maxdegree(poly(Mons), Degree) :-
	is_polynomial(poly(Mons)), !,
	merge_sort(Mons, MonsSorted),
	normalize_poly(MonsSorted, MonsN),
	find_degrees(MonsN, ListD),
	max_list(ListD, Degree).
maxdegree(P, Degree) :-
	is_trad_poly(P),
	as_polynomial(P, Pt),
	maxdegree(Pt, Degree).


find_degrees([],[]) :- !.
find_degrees([m(_,X,_)|Ys], [X|Ms]) :-
	find_degrees(Ys, Ms).

% min_degree(Poly, Degree).
% Returns TRUE if: Degree is the lowest degree amongst the monomials in Poly

mindegree([], 0) :- !.
mindegree(Mon, Degree) :-
	is_monomial(Mon), !,
	mindegree(poly([Mon]), Degree).
mindegree(poly(M), Degree) :-
	is_polynomial(poly(M)), !,
	find_degrees(M, ListD),
	min_list(ListD, Degree).
mindegree(P, Degree) :-
	is_trad_poly(P),
	as_polynomial(P, Pt),
	mindegree(Pt, Degree).


% polyplus(Poly1, Poly2, Result).
% Returns TRUE if: Result is equal to the sum between Poly1 and Poly2 

polyplus(poly([]), poly([]), poly([])) :- !.
polyplus(poly(Mons), poly([]), poly(Mons)) :-
	is_polynomial(poly(Mons)), !.
polyplus(poly([]), poly(Mons), poly(Mons)) :-
	is_polynomial(poly(Mons)), !.
polyplus(poly(Mons1), poly(Mons2), poly(R)) :-
	is_polynomial(poly(Mons1)),
	is_polynomial(poly(Mons2)), !,
	append(Mons1, Mons2, Mapp),
	merge_sort(Mapp, MSorted),
	normalize_poly(MSorted, R).
polyplus(Mon1, Mon2, R) :-
	is_monomial(Mon1),
	is_monomial(Mon2),
	polyplus(poly([Mon1]), poly([Mon2]), R).
polyplus(M1, M2, R) :-
	is_trad_poly(M1),
	is_trad_poly(M2),
	as_polynomial(M1, Mp1),
	as_polynomial(M2, Mp2),
	polyplus(Mp1, Mp2, R).


% polyminus(Poly1, Poly2, Result).
% Returns TRUE if: Result is equal to the difference between Poly1 and Poly2 

polyminus(poly([]), poly([]), poly([])) :- !.
polyminus(poly(Mons), poly([]), poly(Mons)) :-
	is_polynomial(poly(Mons)), !.
polyminus(poly([]), poly(Mons), poly(Mons_inv)) :-
	is_polynomial(poly(Mons)), !,
	invert_coefficient(Mons, Mons_inv).
polyminus(poly(Mons1), poly(Mons2), poly(R)) :-
	is_polynomial(poly(Mons1)),
	is_polynomial(poly(Mons2)), !,
	invert_coefficient(Mons2, Mons2inv),
	append(Mons1, Mons2inv, Mapp),
	merge_sort(Mapp, MSorted),
	normalize_poly(MSorted, R).
polyminus(Mon1, Mon2, R) :-
	is_monomial(Mon1),
	is_monomial(Mon2),
	polyminus(poly([Mon1]), poly([Mon2]), R).
polyminus(M1, M2, R) :-
	is_trad_poly(M1),
	is_trad_poly(M2),
	as_polynomial(M1, Mp1),
	as_polynomial(M2, Mp2),
	polyminus(Mp1, Mp2, R).


% polytimes(Poly1, Poly2, Result).
% Returns TRUE if: Result is equal to the product between Poly1 and Poly2 

polytimes(poly([]), poly([]), poly([])) :- !.
polytimes(poly(_Mons), poly([]), poly([])) :- !.
polytimes(poly([]), poly(_Mons), poly([])) :- !.
polytimes(poly([m(0, _Td, _Vps)]), poly(_Mons), poly([m(0, 0, [])])) :-  !.
polytimes(poly(_Mons), poly([m(0, _Td, _Vps)]), poly([m(0, 0, [])])) :-  !.
polytimes(poly(Mons1), poly(Mons2), poly(M)) :-
	is_polynomial(poly(Mons1)),
	is_polynomial(poly(Mons2)), !,
	multiply_poly(Mons1, Mons2, Mp),
	merge_sort(Mp, Msorted),
	normalize_poly(Msorted, M).
polytimes(Mon1, Mon2, R) :-
	is_monomial(Mon1),
	is_monomial(Mon2),
	polytimes(poly([Mon1]), poly([Mon2]), R).
polytimes(M1, M2, R) :-
	is_trad_poly(M1),
	is_trad_poly(M2), !,
	as_polynomial(M1, Mp1),
	as_polynomial(M2, Mp2),
	polytimes(Mp1, Mp2, R).

multiply_poly([], [], []) :- !.
multiply_poly(_M, [], []) :- !.
multiply_poly([], _M, []) :- !.
multiply_poly(Mons1, Mons2, R) :-
	Mons1 = [M|Ms],
	multiply_mons(M, Mons2, Mp1),
	multiply_poly(Ms, Mons2, Mp2),
	append(Mp1, Mp2, Mapp),
	merge_sort(Mapp, MSorted),
	normalize_poly(MSorted, R).

multiply_mons(_M, [], []) :- !.
multiply_mons(M1, [M2|Ms], [R|Rs]) :-
	M1 = m(C1, _Td1, _Vp1),
	M2 = m(C2, _Td2, _Vp2),
	Cf is C1*C2,
	Cf == 0, !,
	R = [],
	multiply_mons(M1, Ms, Rs).
multiply_mons(M1, [M2|Ms], [R|Rs]) :-
	M1 = m(C1, Td1, Vp1),
	M2 = m(C2, Td2, Vp2),
	Cf is C1*C2,
	Tdf is Td1 + Td2,
	append(Vp1, Vp2, VpfPartial),
	normalize_mons(VpfPartial, Vpf),
	R = m(Cf, Tdf, Vpf),
	multiply_mons(M1, Ms, Rs).




% as_monomial(Expression, Monomial)
% Returns: the ordered representation Monomial of the monomial Expression
%		   where Expression is in the form:
% 		   Monomial = m(Coefficient, TotalDegree, [VarsPower])


as_monomial(X, m(X, 0, [ ])) :- number(X).
as_monomial(-X, m(-1, 1, [v(1, X)])) :- atom(X).
as_monomial(X, m(1, 1, [v(1, X)])) :- atom(X).
as_monomial(^(X, Y), m(1, 0, [ ])) :-
	atom(X), number(Y), Y == 0, !.
as_monomial(-X^Y, m(-1, Y, [v(Y, X)])) :-
	atom(X), number(Y).
as_monomial(^(X, Y), m(1, Y, [v(Y, X)])) :-
	atom(X), number(Y).
as_monomial(*(X, Y), m(C, Td, Vs)) :-
	as_monomial(X, m(Ct1, _Tdt1, _Vpt1)),
	as_monomial(Y, m(Ct2, _Tdt2, _Vpt2)),
	C is Ct1 * Ct2,
	C == 0, !,
	Td is 0,
	Vs = [].
as_monomial(*(X, Y), m(C, Td, Vs)) :-
	as_monomial(X, m(Ct1, Tdt1, Vst1)),
	as_monomial(Y, m(Ct2, Tdt2, Vst2)),
	C is Ct1 * Ct2,
	Td is Tdt1 + Tdt2,
	append(Vst1, Vst2, Vst),
	merge_sort(Vst, VstReady),
	normalize_mons(VstReady, Vsnorm),
	sort(2, @<, Vsnorm, Vs).


normalize_mons([], []).
normalize_mons([v(G, T)], [v(G, T)]).
normalize_mons([v(G, T), v(G1, T1)|Vs], [v(G, T)|Ms]) :-
	T \== T1,!,
	normalize_mons([v(G1,T1)|Vs], Ms).
normalize_mons([v(G, T), v(G1, T1)|Vs], Ms) :-
	T == T1,
	Gf is G + G1,
	Gf == 0,
	normalize_mons(Vs, Ms).
normalize_mons([v(G, T), v(G1, T1)|Vs], Ms) :-
	T == T1,
	Gf is G + G1,
	normalize_mons([v(Gf,T)|Vs], Ms).




% as_polynomial(Expression, Polynomial)
% Returns: the ordered representation Polynomial of the polynomial Expression
%		   where Expression is in the form:
% 		   Polynomial = poly([Monomials])

as_polynomial(X, poly([M])) :-
	number(X), !, as_monomial(X, M).
as_polynomial(-X, poly([M])) :-
	atom(X), !, as_monomial(-1*X, M).
as_polynomial(X, poly([M])) :-
	atom(X), !,as_monomial(X, M).
as_polynomial(*(X, Y), poly([M])) :-
	as_monomial(*(X, Y), M), !.
as_polynomial(-X^Y, poly([M])) :-
	as_monomial(-X^Y, M), !.
as_polynomial(^(X, Y), poly([M])) :-
	as_monomial(^(X, Y), M), !.
as_polynomial(+(X, Y), poly(M)) :-
	as_polynomial(X, poly(Mt)),
	as_polynomial(Y, poly(Pt)), !,
	append(Mt, Pt , Ris),
	merge_sort(Ris, RisReady),
	normalize_poly(RisReady, M).
as_polynomial(-(X, Y), poly(M)) :-
	number(Y), !,
	Yinv is -Y,
	as_polynomial(X,poly(Mt)),
	as_polynomial(Yinv, poly(Pt)),
	append(Mt, Pt , Ris),
	merge_sort(Ris, RisReady),
	normalize_poly(RisReady, M).
as_polynomial(-(X, Y), poly(M)) :-
	invert_coefficient(Y, Yi),
	as_polynomial(X, poly(Mt)),
	append(Mt, [Yi] , Ris),
	merge_sort(Ris, RisReady),
	normalize_poly(RisReady, M).

normalize_poly([], []) :- !.
normalize_poly([m(0, Td, Vp)], [m(0, 0, [])]) :-
	Td \= 0,
	Vp \= [], !.
normalize_poly([m(C, Td, Vp)], [m(C, Td, Vpn)]) :-
	merge_sort(Vp, Vps),
	normalize_mons(Vps, Vpn), !.
normalize_poly([m(C1, _Td, _Vp)|Ms], Zs) :-
	C1 == 0, !,
	normalize_poly(Ms, Zs).
normalize_poly([m(C1, Td1, Vp1), m(C2, Td2, Vp2)|Ms], [M1|Zs]) :-
	merge_sort(Vp1, Vps1),
	merge_sort(Vp2, Vps2),
	normalize_mons(Vps1, Vpn1),
	normalize_mons(Vps2, Vpn2),
	not(equal_mons(Vpn1, Vpn2)), !,
	M1 = m(C1, Td1, Vp1),
	normalize_poly([m(C2, Td2, Vp2)|Ms], Zs).
normalize_poly([m(C1, Td1, Vp1), m(C2, Td2, Vp2)|Ms], Zs) :-
	merge_sort(Vp1, Vps1),
	merge_sort(Vp2, Vps2),
	normalize_mons(Vps1, Vpn1),
	normalize_mons(Vps2, Vpn2),
	equal_mons(Vpn1, Vpn2),
	Td1 == Td2,
	Cf is C1 + C2,
	Cf == 0, !,
	normalize_poly([m(0, 0, [])|Ms], Zs).
normalize_poly([m(C1, Td1, Vp1), m(C2, Td2, Vp2)|Ms], Zs) :-
	equal_mons(Vp1, Vp2),
	Td1 == Td2,
	Cf is C1 + C2,
	normalize_poly([m(Cf, Td1, Vp1)|Ms], Zs).


equal_mons([], []) :- !.
equal_mons([v(G, T)], [v(G, T)]) :- !.
equal_mons([v(G1, T1)|Vs1], [v(G2, T2)|Vs2]) :-
	G1 == G2,
	T1 == T2,
	equal_mons(Vs1, Vs2).



% polyval(Poly, VariablesValue, Value).
% Returns TRUE if: Value is equal to the value of Poly calculated
%				   in the n-dimentional point represented by VariablesValue
% It is assumed that the values in VariablesValue are ordered according
% to the lexicographical order of the variables in Poly


polyval([ ], [ ], 0) :- !.
polyval(m(C, 0, [ ]), [], C) :- !.
polyval(Mon, VV, Value) :-
	is_monomial(Mon),
	polyval(poly([Mon]), VV, Value).
polyval(poly([]), [], 0) :- !.
polyval(poly([m(C, 0, [])]), [], C) :- !.
polyval(poly(Mons), VV, Value) :-
	is_polynomial(poly(Mons)),
	normalize_poly(Mons, MonsN),
	variables(poly(MonsN), Vars),
	length(VV, LVV),
	length(Vars, LV),
	LV == LVV,
	evaluate_poly(MonsN, Vars, VV, Value).
polyval(P, VV, Value) :-
	is_trad_poly(P),
	as_polynomial(P, Pt),
	polyval(Pt, VV, Value).


evaluate_poly([], _, _, 0) :- !.
evaluate_poly([m(C, _, [ ])], _, _, C) :- !.
evaluate_poly(Mons, [V1|Vs], [VV1|VVs], Value) :-
	Mons = [m(C, Td, Vps)|Ms],
	Vps == [], !,
	evaluate_poly([m(C, Td, Vps)], Vs, VVs, Rp1),
	evaluate_poly(Ms, [V1|Vs], [VV1|VVs], Rp2),
	Value is Rp1 + Rp2.
evaluate_poly(Mons, [V1|Vs], [VV1|VVs], Value) :-
	Mons = [m(C, Td, [v(G, X)|Zs])|Ms],
	X == V1, !,
	evaluate_poly([m(C, Td, Zs)], Vs, VVs, Rp1),
	evaluate_poly(Ms, [V1|Vs], [VV1|VVs], Rp2),
	Value is Rp1*VV1^G + Rp2.
evaluate_poly(Mons, [V1|Vs], [VV1|VVs], Value) :-
	Mons = [m(C, Td, [v(G, X)|Zs])|Ms],
	X \= V1,
	evaluate_poly([m(C, Td, [v(G, X)|Zs])], Vs, VVs, Rp1),
	evaluate_poly(Ms, [V1|Vs], [VV1|VVs], Rp2),
	Value is Rp1 +Rp2.


% pprint_polynomial(Polynomial).
% Returns TRUE if: prints on screen the human-readable version of Polynomial

pprint_polynomial(poly([])) :- !.
pprint_polynomial(poly([m(0, 0, [])])) :-
	write(0), !.
pprint_polynomial(poly(Mons)) :-
	is_polynomial(poly(Mons)),
	normalize_poly(Mons, MonsN),
	MonsN = [m(C, _Td, Vps)|Ms],
	C > 0,
	C == 1, !,
	write(+),
	pprint_polynomial(Vps),
	pprint_polynomial(poly(Ms)).
pprint_polynomial(poly(Mons)) :-
	is_polynomial(poly(Mons)),
	normalize_poly(Mons, MonsN),
	MonsN = [m(C, _Td, Vps)|Ms],
	C > 0,
	C \= 1, !,
	write(+), write(C),
	pprint_polynomial(Vps),
	pprint_polynomial(poly(Ms)).
pprint_polynomial(poly(Mons)) :-
	is_polynomial(poly(Mons)),
	normalize_poly(Mons, MonsN),
	MonsN = [m(C, _Td, Vps)|Ms],
	C < 0, !,
	write(C),
	pprint_polynomial(Vps),
	pprint_polynomial(poly(Ms)).
pprint_polynomial([]) :-  !.
pprint_polynomial([v(1, X)|Vs]) :-
	write(X), !,
	pprint_polynomial(Vs).
pprint_polynomial([v(G, X)|Vs]) :-
	write(X^G), !,
	pprint_polynomial(Vs).
pprint_polynomial(Mon) :-
	is_monomial(Mon),
	pprint_polynomial(poly([Mon])).
pprint_polynomial(Poly) :-
	is_trad_poly(Poly),
	as_polynomial(Poly, P),
	pprint_polynomial(P).




%--------- MERGE SORT




merge_sort([], []) :- !.
merge_sort([X], [X]) :- !.
% algoritmo merge_sort dedicato ai monomi
merge_sort(List, Sorted) :-
	List = [v(_, _), v(_, _)|_], !,
	divide(List, L1, L2),
	merge_sort(L1, Sorted1),
	merge_sort(L2, Sorted2),
	merge(Sorted1, Sorted2, Sorted).
% algoritmo merge_sort dedicato ai polinomi
merge_sort(List, Sorted) :-
	List = [m(_, _, Vp1), m(_, _, Vp2)|_],
	Vp1 == [],
	Vp2 == [], !,
	divide(List, L1, L2),
	merge_sort(L1, Sorted1),
	merge_sort(L2, Sorted2),
	merge(Sorted1, Sorted2, Sorted).
merge_sort(List, Sorted) :-
	List = [m(_, _, Vp1), m(_, _, Vp2)|_],
	Vp1 \== [],
	Vp2 \== [], !,
	Vp1 = [v(_, _)|_],
	Vp2 = [v(_, _)|_],
	divide(List, L1, L2),
	merge_sort(L1, Sorted1),
	merge_sort(L2, Sorted2),
	merge(Sorted1, Sorted2, Sorted).
merge_sort(List, Sorted) :-
	List = [m(_, _, Vp1), m(_, _, Vp2)|_],
	Vp1 == [],
	Vp2 \== [], !,
	Vp2 = [v(_, _)|_],
	divide(List, L1, L2),
	merge_sort(L1, Sorted1),
	merge_sort(L2, Sorted2),
	merge(Sorted1, Sorted2, Sorted).
merge_sort(List, Sorted) :-
	List = [m(_, _, Vp1), m(_, _, Vp2)|_],
	Vp1 \== [],
	Vp2 == [],
	Vp1 = [v(_, _)|_],
	divide(List, L1, L2),
	merge_sort(L1, Sorted1),
	merge_sort(L2, Sorted2),
	merge(Sorted1, Sorted2, Sorted).


merge([], L, L) :- !.
merge(L, [], L) :-
	L \= [], !.
% algoritmo merge dedicato ai monomi
merge([v(G1, X)|T1], [v(G2, Y)|T2], [v(G1, X)|T]) :-
	X @=< Y, !,
	merge(T1, [v(G2, Y)|T2], T).
merge([v(G1, X)|T1], [v(G2, Y)|T2], [v(G2, Y)|T]) :-
	X @> Y, !,
	merge([v(G1,X)|T1], T2, T).
% algoritmo merge dedicato ai polinomi
merge([m(C1, Td1, Vp1)|Z1], [m(C2, Td2, Vp2)|Z2], [M1|Z]) :-
	Td1 < Td2, !,
	M1 = m(C1, Td1, Vp1),
	merge(Z1, [m(C2, Td2, Vp2)|Z2], Z).
merge([m(C1, Td1, Vp1)|Z1], [m(C2, Td2, Vp2)|Z2], [M2|Z]) :-
	Td1 > Td2, !,
	M2 = m(C2, Td2, Vp2),
	merge([m(C1, Td1, Vp1)|Z1], Z2,  Z).
merge([m(C1, Td1, Vp1)|Z1], [m(C2, Td2, Vp2)|Z2], [M1|Z]) :-
	Td1 == Td2,
	Vp1 == Vp2, !,
	M1 = m(C1, Td1, Vp1),
	merge(Z1, [m(C2, Td2, Vp2)|Z2], Z).
merge([m(C1, Td1, Vp1)|Z1], [m(C2, Td2, Vp2)|Z2], [M1|Z]) :-
	Td1 == Td2,
	comesfirst(Vp1, Vp2), !,
	M1 = m(C1, Td1, Vp1),
	merge(Z1, [m(C2, Td2, Vp2)|Z2], Z).
merge([m(C1, Td1, Vp1)|Z1], [m(C2, Td2, Vp2)|Z2], [M2|Z]) :-
	Td1 == Td2,
	not(comesfirst(Vp1, Vp2)),
	M2 = m(C2, Td2, Vp2),
	merge([m(C1, Td1, Vp1)|Z1], Z2, Z).


comesfirst([], []) :- !.
comesfirst([v(_C1, T1)|_Vs1], [v(_C2, T2)|_Vs2]) :-
	T1 @< T2, !.
comesfirst([v(C1, T1)|_Vs1], [v(C2, T2)|_Vs2]) :-
	T1 == T2,
	C1 < C2, !.
comesfirst([v(C1, T1)|Vs1], [v(C2, T2)|Vs2]) :-
	T1 == T2,
	C1 == C2,
	comesfirst(Vs1, Vs2).



divide(L, L1, L2) :- halve(L, L1, L2).

halve(L, A, B) :- hv(L, [], A, B).

hv(L, L, [], L) :- !.      % per liste di lunghezza pari
hv(L, [_|L], [], L) :- !.  % per liste di lunghezza dispari
hv([H|T], Acc, [H|L], B) :- hv(T, [_|Acc], L, B).



% SUPPORT FUNCTIONS

is_trad_poly(P) :-
	number(P), !.
is_trad_poly(-P) :-
	atom(P).
is_trad_poly(P) :-
	atom(P), !.
is_trad_poly(^(X, Y)) :-
	atom(X),
	number(Y), !.
is_trad_poly(-X^Y) :-
	atom(X),
	number(Y).
is_trad_poly(*(X, Y)) :-
	is_trad_poly(X),
	is_trad_poly(Y), !.
is_trad_poly(+(X, Y)) :-
	is_trad_poly(X),
	is_trad_poly(Y), !.
is_trad_poly(-(X, Y)) :-
	is_trad_poly(X),
	is_trad_poly(Y).


invert_coefficient([], []).
invert_coefficient(m(C, Td, Vps), m(-C, Td, Vps)).
invert_coefficient([Mon|Ms], [Moninv|Msi]) :-
	Mon = m(C, Td, Vps),
	Moninv = m(-C, Td, Vps),
	invert_coefficient(Ms, Msi).
invert_coefficient(M, Mi) :-
	as_monomial(M, Mp),
	Mp = m(C, Td, Vps),
	Mi = m(-C, Td, Vps).


is_monomial(m(_C,TD,VPs)) :-
	integer(TD),
	TD >= 0,
	is_list(VPs).

is_varpower(v(Power,VarSymbol)) :-
	integer(Power),
	Power >= 0,
	atom(VarSymbol).

is_polynomial(poly(Monomials)) :-
	is_list(Monomials),
	foreach(member(M, Monomials), is_monomial(M)).










