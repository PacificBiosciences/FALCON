import networkx as nx
from pbcore.io import FastaReader

read_fasta = "preads4falcon.fasta"
edge_data_file = "sg_edges_list"
utg_data_file = "utg_data"
ctg_data_file = "ctg_paths"

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

if __name__ == "__main__":    

   

    reads_in_layout = set()
    with open(edge_data_file) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62"""
            v, w, rid, s, t, aln_score, idt, type_ = l
            if type_ != "G":
                continue
            r1 = v.split(":")[0]
            reads_in_layout.add(r1)
            r2 = w.split(":")[0]
            reads_in_layout.add(r2)

    seqs = {}
    # load all p-read name into memory
    f = FastaReader(read_fasta)
    for r in f:
        if r.name not in reads_in_layout:
            continue
        seqs[r.name] = r.sequence.upper()


    edge_data = {}
    with open(edge_data_file) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62"""
            v, w, rid, s, t, aln_score, idt, type_ = l

            if type_ != "G":
                continue
            r1 = v.split(":")[0]
            reads_in_layout.add(r1)
            r2 = w.split(":")[0]
            reads_in_layout.add(r2)

            s = int(s)
            t = int(t)
            aln_score = int(aln_score)
            idt = float(idt)

            if s < t:
                e_seq = seqs[ rid ][ s:t ]
            else:
                e_seq = "".join([ RCMAP[c] for c in seqs[ rid ][ s:t:-1 ] ])
            edge_data[ (v, w) ] = ( rid, s, t, aln_score, idt, e_seq )

    utg_data = {}
    with open(utg_data_file) as f:
        for l in f:
            l = l.strip().split()
            s, v, t, type_, is_spur, length, score, path_or_edges = l
            if is_spur == "+":
                continue
            length = int(length)
            score = int(score)
            if type_ == "path":
                path_or_edges = path_or_edges.split("~")
            else:
                path_or_edges = [ tuple(e.split("~")) for e in path_or_edges.split(".") ]
            utg_data[ (s,v,t) ] = type_, length, score, path_or_edges

    p_ctg_out = open("p_ctg.fa","w")
    a_ctg_out = open("a_ctg.fa","w")
    p_ctg_t_out = open("p_ctg_titling_path","w")
    a_ctg_t_out = open("a_ctg_titling_path","w")
    layout_ctg = set()

    with open(ctg_data_file) as f:
        for l in f:
            l = l.strip().split()
            ctg_id, c_type_, i_utig, t0, length, score, utgs = l
            ctg_id = ctg_id
            s0 = i_utig.split("~")[0]

            if (reverse_end(t0), reverse_end(s0)) in layout_ctg:
                continue
            else:
                layout_ctg.add( (s0, t0) ) 

            ctg_label = i_utig+"~"+t0
            length = int(length)
            utgs = utgs.split(".")
            one_path = []
            total_score = 0
            total_length =0
            a_ctg_data = []
            for utg in utgs:
                s,v,t  = utg.split("~")
                type_, length, score, path_or_edges = utg_data[ (s,v,t) ]
                total_score += score
                total_length += length
                if type_ == "path":
                    if len(one_path) != 0:
                        one_path.extend ( path_or_edges[1:] )
                    else:
                        one_path.extend ( path_or_edges )
                else:
                    c_graph = nx.DiGraph()

                    for vv, ww in path_or_edges:
                        c_graph.add_edge( vv, ww, score = 10000000 - edge_data[ (vv, ww) ][2]  )

                    shortest_path = nx.shortest_path( c_graph, s, t, weight = "score" )

                    if len(one_path) != 0:
                        one_path.extend ( shortest_path[1:] )
                    else:
                        one_path.extend ( shortest_path )

                    while 1:
                        if s == t:
                            break
                        n0 = shortest_path[0]
                        for n1 in shortest_path[1:]:
                            c_graph.remove_edge(n0, n1)
                            n0 = n1
                        try:
                            shortest_path = nx.shortest_path( c_graph, s, t, weight = "score" )
                            a_ctg_data.append( (s, t, shortest_path) )
                        except nx.exception.NetworkXNoPath:
                            break
                        if len(shortest_path) < 2:
                            break



            if len(one_path) == 0:
                continue

            one_path_edges = zip(one_path[:-1], one_path[1:])

            sub_seqs = []
            for vv, ww in one_path_edges:
                rid, s, t, aln_score, idt, e_seq = edge_data[ (vv, ww) ]
                sub_seqs.append( e_seq )
                print >> p_ctg_t_out, "%s %s %s %s %d %d %d %0.2f" % (ctg_id, vv, ww, rid, s, t, aln_score, idt) 
            print >> p_ctg_out, ">%s %s %s %d %d" % (ctg_id, ctg_label, c_type_, total_length, total_score)
            print >> p_ctg_out, "".join(sub_seqs)

            a_id = 1
            for v, w, atig_path in a_ctg_data:
                atig_path_edges = zip(atig_path[:-1], atig_path[1:])
                sub_seqs = []
                total_length = 0
                total_score = 0
                for vv, ww in atig_path_edges:
                    rid, s, t, aln_score, idt, e_seq = edge_data[ (vv, ww) ]
                    sub_seqs.append( e_seq )
                    total_length += abs(s-t)
                    total_score += aln_score
                    print >> a_ctg_t_out, "%s-%03d %s %s %s %d %d %d %0.2f" % (ctg_id, a_id, vv, ww, rid, s, t, aln_score, idt) 
                print >> a_ctg_out, ">%s-%03d %s %s %d %d" % (ctg_id, a_id, v, w, total_length, total_score )
                print >> a_ctg_out, "".join(sub_seqs)
                a_id += 1

    a_ctg_out.close()
    p_ctg_out.close()
    a_ctg_t_out.close()
    p_ctg_t_out.close()

