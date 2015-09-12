from falcon_kit import kup, falcon, DWA
from falcon_kit.fc_asm_graph import AsmGraph
import networkx as nx

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))
def rc(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])

def get_aln_data(t_seq, q_seq):
    aln_data = []
    K = 8
    seq0 = t_seq
    lk_ptr = kup.allocate_kmer_lookup( 1 << (K * 2) )
    sa_ptr = kup.allocate_seq( len(seq0) )
    sda_ptr = kup.allocate_seq_addr( len(seq0) )
    kup.add_sequence( 0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)
    q_id = "dummy"
    
    kmer_match_ptr = kup.find_kmer_pos_for_seq(q_seq, len(q_seq), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range_ptr = kup.find_best_aln_range(kmer_match_ptr, K, K*5, 12)
    aln_range = aln_range_ptr[0]
    x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count)] )
    kup.free_kmer_match(kmer_match_ptr)
    s1, e1, s2, e2 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2
    
    if e1 - s1 > 100:

        alignment = DWA.align(q_seq[s1:e1], e1-s1,
                              seq0[s2:e2], e2-s2,
                              1500,1)

        if alignment[0].aln_str_size > 100:
            aln_data.append( ( q_id, 0, s1, e1, len(q_seq), s2, e2, len(seq0), alignment[0].aln_str_size, alignment[0].dist ) )
            aln_str1 = alignment[0].q_aln_str
            aln_str0 = alignment[0].t_aln_str

        DWA.free_alignment(alignment)

    kup.free_kmer_lookup(lk_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_seq_addr_array(sda_ptr)
    return aln_data, x, y

def main(argv=None):
  G_asm = AsmGraph("sg_edges_list", "utg_data", "ctg_paths")
  G_asm.load_sg_seq("preads4falcon.fasta")

  utg_out = open("utgs.fa","w")


  for utg in G_asm.utg_data:
    s,t,v  = utg
    type_, length, score, path_or_edges = G_asm.utg_data[ (s,t,v) ]
    if type_ == "simple":
        path_or_edges = path_or_edges.split("~")
        seq = G_asm.get_seq_from_path( path_or_edges )
        print >> utg_out, ">%s~%s~%s-%d %d %d" % (s, v, t, 0, length, score ) 
        print >> utg_out, seq

    if type_ == "compound":

        c_graph = nx.DiGraph()

        all_alt_path = []
        path_or_edges = [ c.split("~") for c in path_or_edges.split("|")]
        for ss, vv, tt in path_or_edges:
            type_, length, score, sub_path = G_asm.utg_data[ (ss,tt,vv) ]
             
            sub_path = sub_path.split("~")
            v1 = sub_path[0]
            for v2 in sub_path[1:]:
                c_graph.add_edge( v1, v2, e_score = G_asm.sg_edges[ (v1, v2) ][1]  )
                v1 = v2
        
        shortest_path = nx.shortest_path( c_graph, s, t, "e_score" )
        score = nx.shortest_path_length( c_graph, s, t, "e_score" )
        all_alt_path.append( (score, shortest_path) )
        

        #a_ctg_data.append( (s, t, shortest_path) ) #first path is the same as the one used in the primary contig
        while 1:
            if s == t:
                break
            n0 = shortest_path[0]
            for n1 in shortest_path[1:]:
                c_graph.remove_edge(n0, n1)
                n0 = n1
            try:
                shortest_path = nx.shortest_path( c_graph, s, t, "e_score" )
                score = nx.shortest_path_length( c_graph, s, t, "e_score" )
                #a_ctg_data.append( (s, t, shortest_path) )
                all_alt_path.append( (score, shortest_path) )

            except nx.exception.NetworkXNoPath:
                break
            #if len(shortest_path) < 2:
            #    break

        all_alt_path.sort()
        all_alt_path.reverse()
        shortest_path = all_alt_path[0][1]

        
        score, atig_path = all_alt_path[0]

        atig_output = []

        atig_path_edges = zip(atig_path[:-1], atig_path[1:])
        sub_seqs = []
        total_length = 0
        total_score = 0
        for vv, ww in atig_path_edges:
            r, aln_score, idt, typs_  = G_asm.sg_edges[ (vv, ww) ]
            e_seq  = G_asm.sg_edge_seqs[ (vv, ww) ]
            rid, ss, tt = r
            sub_seqs.append( e_seq )
            total_length += abs(ss-tt)
            total_score += aln_score

        base_seq = "".join(sub_seqs)
        atig_output.append( (s, t, atig_path, total_length, total_score, base_seq, atig_path_edges, 1, 1) )


        duplicated = True
        for score, atig_path in all_alt_path[1:]:
            atig_path_edges = zip(atig_path[:-1], atig_path[1:])
            sub_seqs = []
            total_length = 0
            total_score = 0
            for vv, ww in atig_path_edges:
                r, aln_score, idt, type_ = G_asm.sg_edges[ (vv, ww) ]
                e_seq  = G_asm.sg_edge_seqs[ (vv, ww) ]
                rid, ss, tt = r
                sub_seqs.append( e_seq )
                total_length += abs(ss-tt)
                total_score += aln_score

            seq = "".join(sub_seqs)

            aln_data, x, y = get_aln_data(base_seq, seq)
            if len( aln_data ) != 0:
                idt =  1.0-1.0*aln_data[-1][-1] / aln_data[-1][-2]
                cov = 1.0*(aln_data[-1][3]-aln_data[-1][2])/aln_data[-1][4]
                if idt < 0.96 or cov < 0.98:
                    duplicated = False
                    atig_output.append( (s, t, atig_path, total_length, total_score, seq, atig_path_edges, idt, cov) )
            else:
                duplicated = False
                atig_output.append( (s, t, atig_path, total_length, total_score, seq, atig_path_edges, 0, 0) )

        #if len(atig_output) == 1:
        #    continue

        sub_id = 0
        for data in atig_output:
            v0, w0, tig_path, total_length, total_score, seq, atig_path_edges, a_idt, cov = data
            print >> utg_out, ">%s~%s~%s-%d %d %d" % (v0, "NA", w0, sub_id,  total_length, total_score ) 
            print >> utg_out, seq
            sub_id += 1
