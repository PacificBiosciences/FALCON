from pbcore.io import FastaReader
import random
random.seed(42)


def get_simulated_seq( seq ):
    sim_seq = []
    ins_rate = 0.004
    del_rate = 0.004
    mis_rate = 0.002


    for c in seq:
        while 1:
            if random.uniform(0, 1) < ins_rate:
                sim_seq.append( random.choice(["A","C","G","T"]))
            else:
                break
        if random.uniform(0, 1) < del_rate:
            continue
        if random.uniform(0, 1) < mis_rate:
            sim_seq.append( random.choice(["A","C","G","T"]))
            sim_seq.append(c)
            continue
        sim_seq.append(c)
    sim_seq = "".join(sim_seq)
    return sim_seq



f = FastaReader("ecoli_k12_MG1655.fasta")
for r in f:
    seq = r.sequence
    break
seq = seq + seq[:40000] 
s = 0
rc_map = dict(zip("ACGT","TGCA"))
for i in range(12000):
    s += 1000
    if s + 40000 < len(seq):
        subs = seq[s:s+10000]
        if hash(subs) % 2 == 0:
            subs = "".join([rc_map[c] for c in subs[::-1]])
            print ">%05d_r" % i
            print subs
        else:
            print ">%05d" % i
            print subs

f = FastaReader("DH10B.fa")
for r in f:
    seq = r.sequence
    break
seq = seq + seq[:40000] 
s = 0
rc_map = dict(zip("ACGT","TGCA"))
i_map = dict(zip("ACGT","ACGT"))
for i in range(12000):
    s += 1000
    if s + 40000 < len(seq):
        subs = seq[s:s+10000]
        if hash(subs) % 2 == 0:
            subs = "".join([rc_map.get(c,"A") for c in subs[::-1]])
            print ">%05d_r" % (i + 12000)
            print subs
        else:
            subs = "".join([i_map.get(c,"A") for c in subs])
            print ">%05d" % (i + 12000)
            print subs
