import subprocess
from pbcore.io import FastaReader

def get_matches(seq0, seq1):
    with open("tmp_seq0.fa","w") as f:
        print >>f, ">seq0"
        print >>f, seq0
    with open("tmp_seq1.fa","w") as f:
        print >>f, ">seq1"
        print >>f, seq1
    mgaps_out=subprocess.check_output("mummer -maxmatch -c -b -l 24 tmp_seq0.fa tmp_seq1.fa | mgaps ", stderr = open("/dev/null"), shell=True)

    matches = []
    cluster = []
    for l in mgaps_out.split("\n"):
        l = l.strip().split()
        if len(l) == 0:
            continue
        if l[0] == ">":
            seq_id = l[1]
            
            if len(cluster) != 0:
                matches.append(cluster)
            
            cluster = []
            continue
        if l[0] == "#":
            if len(cluster) != 0:
                matches.append(cluster)            
            cluster = []
            continue
        len_ = int(l[2])
        r_s = int(l[0])
        q_s = int(l[1])
        r_e = r_s + len_
        q_e = q_s + len_
        cluster.append( ((r_s, r_e), (q_s, q_e)) )
    if len(cluster) != 0:
        matches.append(cluster)
    return matches


u_edges = {}
with open("./unit_edges.dat") as f:
    for l in f:
        v, w, path, seq = l.strip().split()
        u_edges.setdefault( (v, w), [] )
        u_edges[ (v, w) ].append( (path, seq) )
        

p_tig_path = {}
a_tig_path = {}
with open("primary_tigs_paths_c") as f:
    for l in f:
        l = l.strip().split()
        id_ = l[0][1:]
        path = l[1:]
        p_tig_path[id_] = path

with open("all_tigs_paths") as f:
    for l in f:
        l = l.strip().split()
        id_ = l[0][1:]
        path = l[1:]
        a_tig_path[id_] = path

p_tig_seqs = {}
for r in FastaReader("primary_tigs_c.fa"):
    p_tig_seqs[r.name] = r.sequence

a_tig_seqs = {}
for r in FastaReader("all_tigs.fa"):
    a_tig_seqs[r.name.split()[0]] = r.sequence

p_tig_to_node_pos = {}
node_pos = []
with open("primary_tigs_node_pos_c") as f:
    for l in f:
        l = l.strip().split()
        p_tig_to_node_pos.setdefault( l[0], [])
        p_tig_to_node_pos[l[0]].append( (l[1], int(l[2])))

duplicate_a_tigs = []
with open("a_nodup.fa","w") as out_f:
    for p_tig_id in p_tig_path:
        main_path = p_tig_path[p_tig_id]
        main_path_nodes = set(main_path[:])
        p_tig_seq = p_tig_seqs[p_tig_id]
        a_node = []
        a_node_range = []
        a_node_range_map = {}
        node_to_pos = dict( p_tig_to_node_pos[p_tig_id] )
        for id_ in a_tig_path:
            if id_[:4] != p_tig_id[:4]:
                continue
            if id_.split("-")[1] == "0000":
                continue
            
            a_path = a_tig_path[id_]
            if a_path[0] in main_path_nodes and a_path[-1] in main_path_nodes:
                #print p_tig_id, id_, a_path[0], a_path[-1]
                s, e = node_to_pos[a_path[0]], node_to_pos[a_path[-1]]
                p_seq = p_tig_seq[s:e]
                a_seq = a_tig_seqs[id_] 
                seq_match = get_matches(p_seq, a_seq)
                if len(seq_match) > 1:
                    print >>out_f, ">"+id_
                    print >>out_f,  a_seq
                    continue
                try:
                    r_s, r_e = seq_match[0][0][0][0], seq_match[0][-1][0][1]
                except:
                    print "XXX", seq_match
                if 1.0* (r_e - r_s) / (e - s) > 98:
                    print >>out_f, ">"+id_
                    print >>out_f, a_seq
                    continue
                duplicate_a_tigs.append(id_)

