SetLogFile("main.log");
load "X0p_NiceModel.m";
load "auxiliary.m";
load "Chabauty_MWSieve_new.m";
load "Chabauty_MWSieve_PrimeSelector.m";
load "Saturation.m";
SetDebugOnError(true);

N := 43;

C := CuspForms(N);
printf "Dimension of CuspForms(%o) is: %o\n", N, Dimension(C);

//  Check rk J_0(N)(Q) = rk J_0(N)^+(Q)
if not IsRankOfALQuotEqual(N) then
	error "One needs rk J_0(N)(Q) = rk J_0(N)^+(Q) for our algorithm to work.";
else
	printf "rk J_0(N)(Q) = rk J_0(N)^+(Q).\n";
end if;

//we find models for X_0(N) and X_0(N)/w_N

ALN := AtkinLehnerOperator(C, N);
NN := Nullspace(Matrix(ALN - 1));

printf "Dimension of eigenspace lambda = +1 for w_%o is: %o\n", N, Dimension(NN);

NNc := Nullspace(Matrix(ALN + 1));

printf "Dimension of eigenspace lambda = -1 for w_%o is: %o\n", N, Dimension(NNc);

BN  := [&+[(Integers()!(2*Eltseq(Basis(NN )[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(NN)]];
BNc := [&+[(Integers()!(2*Eltseq(Basis(NNc)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(NNc)]];

XN, XN_Cusps := modformeqns(BNc, BN, N, 500, false);
printf "Nice model for X_0(%o) is: %o\n\n", N, XN;
RR<[u]> := CoordinateRing(AmbientSpace(XN));
n := Dimension(AmbientSpace(XN));

H := DiagonalMatrix(RationalField(), Dimension(C), [+1 : i in [1..Dimension(NNc)]] cat [-1 : i in [1..Dimension(NN)]]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]];
wN := iso<XN -> XN | rows, rows>;
printf "w_%o on X%o is given by: %o\n", N, N, wN;

printf "Genus of X_0(%o) is %o\n", N, Genus(XN);
printf "Genus of X_0(%o)/w_%o is %o\n", N, N, Dimension(NN);

printf "We have found these points on X_0(%o):\n%o\n", N, XN_Cusps;

// this can be used to compute a multiple of the torsion of J_0(N)(Q)
// TorsionBound(XN, PrimeDivisors(N));

l := 3;
r := rank_J0Nplus(N);
printf "rk J_0(N)^+(Q) = %o\n", r;
gens := [Divisor(c1) - Divisor(c2) : c1,c2 in XN_Cusps | c1 ne c2];
ptsQ := RationalPoints(XN : Bound := 100);
pts := [pt : pt in ptsQ | pt notin XN_Cusps]; // non-cuspidal Q-points
Append(~gens, Divisor(pts[1]) - Divisor(XN_Cusps[1]));
// or use NonCuspidalQRationalPoints(SmallModularCurve(43),43); for small N
printf "%o divides the index: %o\n", l, FindIntegerCoprimeToIndex(XN, l, r, gens : pprimes := PrimesUpTo(30));

//known degree 1 places
pls1 := [Place(XN_Cusps[i]) : i in [1..#XN_Cusps]];

//known degree 2 places, we conjecture all are pullbacks
pls2:=[];

deg2:=[];
deg2pb:=[];

for i in [1..#pls1] do
	for j in [i..#pls1] do
		Append(~deg2, 1*pls1[i] + 1*pls1[j]);
		if wN(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			Append(~deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;

printf "We have found %o points on X_0(%o)^2(Q).\n", #deg2, N;
printf "%o of them are pullbacks of rationals from X_0(%o)/w_%o.\n", #deg2pb, N, N;

//Finally, we do the sieve.
//A := AbelianGroup([order]);
//A := Ksub;
//divs := [Dtor];
//D := [Divisor(divsNew[i]) - Divisor(XN_Cusps[1]) : i in [1..#divsNew]];
//divs := [&+[coeffs[i] * D[i] : i in [1..#coeffs]] : coeffs in bas];

A, divs := GetTorsion(N, XN, XN_Cusps);
genusC := Dimension(NN);
bp := deg2pb[1];
wNMatrix := Matrix(wN);

primes := PrimesInInterval(3,50); // TODO: find suitable primes
B0, iA0 := sub<A | Generators(A)>;
W0 := {0*A.1};

B, iA, W := MWSieve(XN, wNMatrix, genusC, primes, A, divs, bp, B0, iA0, W0, deg2); // 
//MWSieve(XN, wNMatrix, genusC, primes, A, divs, I, bp, B0, iA0, W0, deg2); // I is a number

printf "\nFor unknown Q in X_0(%o)^2(\Q) we have (1 - w_%o)[Q - bp] is in a coset of %o represented by an element of %o.\n", N, N, B, W;
if #W eq 1 and IsIdentity(W[1]) then
	printf "It follows that if there is an unknown Q in X_0(%o)^2(\Q), then (1 - w_%o)[Q - bp] == 0.\n", N, N;
	printf "This implies that [Q - bp] is fixed by w_%o\n", N;
	printf "Then Q ~ w_%o(Q), which implies that Q = w_%o(Q) because X_0(%o) is not hyperelliptic.\n", N, N, N;
	printf "Then either Q is a pullback, or it is fixed by w_%o pointwise.\n", N;
	printf "If P = (X_i) is fixed by w_%o,\n", N;
	printf "either the first %o coordinates are 0 or the last %o coordinates are 0\n\n", Dimension(NN), Dimension(NNc);

	I := IdentityMatrix(Rationals(), Genus(XN));
	CR<[x]> := CoordinateRing(AmbientSpace(XN));
	// all coordinates where there is a -1 in w_N must be 0 for a point fixed by w_N
	J1 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..Genus(XN)]]> : v in Basis(Kernel(wNMatrix + I))];
	J2 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..Genus(XN)]]> : v in Basis(Kernel(wNMatrix - I))];

	Z1 := Scheme(XN, J1);
	Z2 := Scheme(XN, J2);
	Z := Union(Z1, Z2);
	printf "We find all such P by imposing these conditions and finding points on the scheme:\n%o\n\n", Z;

	pts := PointsOverSplittingField(Z);
	printf "All such P are:\n%o\n", pts;
	pts_of_deg_le_2 := [P : P in pts | forall{i : i in [1..Dimension(AmbientSpace(Z))] | Degree(MinimalPolynomial(P[i])) le 2}];
	ptsclstr := [Cluster(Z, P) : P in pts_of_deg_le_2];
	degrees := [Degree(P) : P in ptsclstr];
	printf "The degrees of these points are %o (or > 2).\n", degrees;
	if forall{deg: deg in degrees | deg ne 2} then
		printf "Hence there are no quadratic points on X_0(%o) not coming from pullbacks of rationals.\n", N;
	else
		error "TODO: Sieve worked, but we still need to analyze quadratic points (there are some).";
	end if;
else 
	error "TODO: Sieve did not prove what we wanted.";
end if;
