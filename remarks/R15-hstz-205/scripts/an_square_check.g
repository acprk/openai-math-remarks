# an_square_check.g  (GAP 4; uses only the generic character tables
#   CharacterTable("Symmetric",n) and CharacterTable("Alternating",n)).
# Written from scratch for remark R15.
#
# For every n in [NMIN..NMAX] and every self-conjugate lambda |- n:
#  * pick an A_n-irreducible chi_+ below Res chi^lambda (no use of GAP's labels),
#  * decompose chi_+^2 directly in the A_n table,
#  * compare with the closed formula of Proposition 1 (note.tex):
#       nu != nu'           : (g + eps*a)/2
#       nu = nu' != lambda  : (g + eps*a)/4      (each of the two constituents)
#       nu = lambda         : (g + 1 + 2*eps)/4 for chi_+,  (g + 1 - 2*eps)/4 for chi_-
#    where g = g(lambda,lambda,nu), a = chi^nu(hhat(lambda)), eps = (-1)^((n-d)/2),
#  * report whether chi^lambda is S_n-universal and whether chi_+ is A_n-universal.
# Run:  gap -q -o 8g an_square_check.g     (NMIN, NMAX set below)
NMIN := 5;; NMAX := 20;;
DiagHooks := function(lam)
  local d, i, lt;
  lt := AssociatedPartition(lam);
  d := Number([1..Length(lam)], i -> lam[i] >= i);
  return List([1..d], i -> (lam[i]-i) + (lt[i]-i) + 1);
end;;
totmis := 0;;
for n in [NMIN..NMAX] do
  tS := CharacterTable("Symmetric", n);; tA := CharacterTable("Alternating", n);;
  irrS := Irr(tS);; irrA := Irr(tA);; cent := SizesCentralizers(tS);;
  cpS := List(ClassParameters(tS), x -> x[2]);;
  chS := List(CharacterParameters(tS), x -> x[2]);;
  resS := List(irrS, x -> RestrictedClassFunction(x, tA));;
  # nu-label of each A_n irreducible: the S_n characters lying over it
  over := List(irrA, psi -> Filtered([1..Length(irrS)], i -> ScalarProduct(tA, resS[i], psi) <> 0));;
  for il in [1..Length(chS)] do
    lam := chS[il];
    if AssociatedPartition(lam) <> lam then continue; fi;
    hh := DiagHooks(lam);; d := Length(hh);; eps := (-1)^((n-d)/2);;
    c0 := Position(cpS, hh);;
    cands := Filtered([1..Length(irrA)], j -> il in over[j]);;
    if Length(cands) <> 2 then Error("split?"); fi;
    jp := cands[1];; jm := cands[2];;
    sq := List([1..Length(irrA[jp])], c -> irrA[jp][c]^2);;
    sq := ClassFunction(tA, sq);;
    mis := 0;; allpos := true;; Sall := true;; bad := [];;
    for j in [1..Length(irrA)] do
      obs := ScalarProduct(tA, sq, irrA[j]);
      inu := over[j][1]; nu := chS[inu];
      g := Sum([1..Length(cent)], c -> irrS[il][c]^2 * irrS[inu][c] / cent[c]);
      a := irrS[inu][c0];
      if AssociatedPartition(nu) <> nu then pred := (g + eps*a)/2;
      elif nu <> lam then pred := (g + eps*a)/4;
      elif j = jp then pred := (g + 1 + 2*eps)/4;
      else pred := (g + 1 - 2*eps)/4; fi;
      if obs <> pred then mis := mis + 1; fi;
      if obs = 0 then allpos := false; Add(bad, [nu, g, a]); fi;
    od;
    for inu in [1..Length(irrS)] do
      g := Sum([1..Length(cent)], c -> irrS[il][c]^2 * irrS[inu][c] / cent[c]);
      if g = 0 then Sall := false; fi;
    od;
    totmis := totmis + mis;
    Print("n=", n, " lambda=", lam, " hhat=", hh, " eps=", eps,
          " formula_mismatches=", mis, " S_n-universal=", Sall,
          " A_n-universal(chi+)=", allpos, " missing[nu,g,a]=", bad, "\n");
  od;
od;
Print("TOTAL formula mismatches: ", totmis, "\n");
QUIT;
