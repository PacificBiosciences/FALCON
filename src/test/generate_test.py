import random
random.seed(42)

"""
For generating the `test.in` files:
```
python3 generate_test.py > test.in
```
"""

def sim_error(seq, pi=0.05, pd=0.05, ps=0.01):
    out_seq = []
    for c in seq:
        while 1:
            r = random.uniform(0,1)
            if r < pi:
                out_seq.append(random.choice(["A","C","G","T"]))
            else:
                break
        r -= pi
        if r < pd:
            continue
        r -= pd
        if r < ps:
            out_seq.append(random.choice(["A","C","G","T"]))
            continue
        out_seq.append(c)
    return "".join(out_seq)

for i in range(5):
    seq = [random.choice(["A","C","G","T"]) for _ in range(1000)]
    print("test{}".format(i), "".join(seq))
    for j in range(20):
        print("rnd{}".format(i), sim_error(seq))
    print("+ +")
print("- -")

