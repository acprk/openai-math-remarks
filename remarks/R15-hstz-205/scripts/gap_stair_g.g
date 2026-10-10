# gap_stair_g.g -- exact g(rho_m,rho_m,lambda) from the generic S_n character table,
# printed as "parts g a" (a = chi^lambda(hhat_m)), for comparison with kron_stair.
for m in [3..6] do
  n := m*(m+1)/2;; t := CharacterTable("Symmetric", n);;
  irr := Irr(t);; cent := SizesCentralizers(t);;
  cp := List(ClassParameters(t), x -> x[2]);; ch := List(CharacterParameters(t), x -> x[2]);;
  rho := Reversed([1..m]);; ir := Position(ch, rho);;
  hh := Filtered(List([1..m], i -> 2*m-1-4*(i-1)), x -> x > 0);; c0 := Position(cp, hh);;
  out := OutputTextFile(Concatenation("gap_g_m", String(m), ".txt"), false);;
  SetPrintFormattingStatus(out, false);
  for i in [1..Length(ch)] do
    g := Sum([1..Length(cent)], c -> irr[ir][c]^2 * irr[i][c] / cent[c]);
    PrintTo(out, JoinStringsWithSeparator(List(ch[i], String), ","), " ", g, " ", irr[i][c0], "\n");
  od;
  CloseStream(out);
od;
QUIT;
