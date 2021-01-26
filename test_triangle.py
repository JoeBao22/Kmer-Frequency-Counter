t = "CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT"
d = {}
for i in range(len(t) - 3 + 1):
    k = t[i: i+3]
    d[k] = d.get(k, 0) + 1
print(d)
