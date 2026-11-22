#!/usr/bin/env python3
"""Write algs.g (GAP input for check_modisom.g) from nilpotent_F2_upto5.json."""
import json

data = json.load(open("nilpotent_F2_upto5.json"))
out = ['LoadPackage("modisom");;', "ALGS := [];;"]
for n in sorted(data, key=int):
    lst = []
    for A in data[n]:
        nn = len(A)
        entries = []
        for i in range(nn):
            for j in range(nn):
                ks = [k + 1 for k in range(nn) if (A[i][j] >> k) & 1]
                if ks:
                    entries.append("[%d,%d,[%s]]" % (i + 1, j + 1, ",".join(map(str, ks))))
        lst.append("[" + ",".join(entries) + "]")
    out.append("ALGS[%s] := [%s];;" % (n, ",\n".join(lst)))
open("algs.g", "w").write("\n".join(out) + "\n")
