# For each listed algebra: ModIsom canonical form and |Aut|.
Read("algs.g");
MakeAlg := function(n, entries)
    local T, e, F, k;
    F := GF(2);
    T := EmptySCTable(n, Zero(F));
    for e in entries do
        SetEntrySCTable(T, e[1], e[2], Concatenation(List(e[3], k -> [One(F), k])));
    od;
    return AlgebraByStructureConstants(F, T);
end;
for n in [1..Length(ALGS)] do
    canos := [];
    sizes := [];
    for entries in ALGS[n] do
        A := MakeAlg(n, entries);
        T := NilpotentTable(A);
        R := CanoFormWithAutGroupOfTable(T);
        C := ShallowCopy(R.cano);
        Unbind(C.iso);
        Add(canos, C);
        Add(sizes, R.auto.size);
    od;
    PrintTo(Concatenation("autsizes_", String(n), ".txt"), sizes, "\n");
    Print("dim ", n, ": listed ", Length(ALGS[n]), ", distinct canonical forms ",
          Length(Set(canos)), "\n");
od;
QUIT;
