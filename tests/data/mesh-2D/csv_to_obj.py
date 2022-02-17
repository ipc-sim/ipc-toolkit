import pandas

V = pandas.read_csv("V_t0.csv", header=None).to_numpy()
E = pandas.read_csv("E.csv", header=None).to_numpy()

for v in V:
    print("v {:g} {:g} {:g}".format(*v))
for e in E:
    print("l {:d} {:d}".format(*(e)))
