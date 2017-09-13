import os
import sys

import networkx as nx
from fc_asm_graph import AsmGraph
from FastaReader import FastaReader

GFA_H_TAG = 'H'
GFA_S_TAG = 'S'
GFA_L_TAG = 'L'
GFA_P_TAG = 'P'
GFA_ORIENT_FWD = '+'
GFA_ORIENT_REV = '-'
GFA_SEQ_UNKNOWN = '*'
GFA_LINK_CIGAR_UNKNOWN = '*'

class GFAGraph:
    def __init__(self):
        self.paths = {}
        self.read_pairs = set()
        self.read_in_graph = set()
        self.edges = {}

    def add_tiling_path(self, ctg_path, ctg_name):
        """
        For a given tiling path, this method adds nodes
        and edges between every connected node.
        If an edge already exists, it will not be overwritten.
        Instead, it will be pdated with tiling path info
        (e.g. sequence tiling begin and end coordinate, overlap
        identity, etc.)
        ctg_path is a list containing a tiling path.
        ctg_name is the name of the contig that this tiling
        path describes.
        """
        self.paths[ctg_name] = ctg_path
        for v, w, b, e, l, idt, etype in ctg_path:
            self.add_read_from_node(v)
            self.add_read_from_node(w)
            self.add_or_update_edge(v, w, '*', l, idt, b, e, None, None, ctg_name, etype)

    def add_asm_graph(self, asm_graph):
        """
        Takes a Falcon AsmGraph object and adds all
        edges in the final assembly graph to `self`.
        It also adds nodes incident with those edges.
        """
        for v, w in asm_graph.sg_edges:
            edge_data = asm_graph.sg_edges[(v, w)]
            if edge_data[-1] == 'G':
                overlap_begin = int(edge_data[0][1])
                overlap_end = int(edge_data[0][2])
                overlap_length = int(edge_data[1])
                overlap_idt = float(edge_data[2])
                self.add_or_update_edge(v, w, '*', overlap_length, overlap_idt, overlap_begin, overlap_end, None, None, None, None)

    def add_nx_string_graph(self, nx_sg):
        """
        Takes an Networkx object specifying an
        assembly string graph and adds all edges and
        nodes to `self`.
        """
        for node in nx_sg.nodes():
            self.add_read_from_node(node)
        for v, w in nx_sg.edges():
            edata = nx_sg.get_edge_data(v, w);
            src_graph = edata.get('src', None)
            cross_phase = edata.get('cross_phase', None)
            self.add_or_update_edge(v, w, '*', None, None, None, None, cross_phase, src_graph, None, None)

    def write_gfa_v1(self, fp_out, preads_file, contig_files, write_reads, write_contigs):
        """
        Writes the nodes, edges and paths in a GFA-v1 format.
        preads_file is a string with a path to the preads4falcon.fasta file.
        contig_files is a list of strings with primary/associate contigs or haplotigs.
        If write_reads is True, read sequences will be output to GFA as well.
        If write_contigs is True, all contigs which have a corresponding tiling
        path in the GFAGraph object will have their sequences output.
        """

        # Output version.
        fp_out.write(GFA_H_TAG + '\tVN:Z:1.0\n')

        seq_len_map = {}

        # Output reads on the fly.
        # Perhaps solve this via yield, but we also want the sequence length map generated simultaneously as well.
        if self.read_in_graph:
            found_reads = set()
            f = FastaReader(preads_file)
            for r in f:
                rname = r.name.split()[0]
                seq_len_map[rname] = len(r.sequence)
                if rname in self.read_in_graph:
                    fp_out.write('\t'.join([GFA_S_TAG, rname, GFA_SEQ_UNKNOWN if (not write_reads) else r.sequence, 'LN:i:%s' % len(r.sequence)]) + '\n')
                    found_reads.add(rname)
            if len(found_reads) != len(self.read_in_graph):
                raise Exception('Not all reads were found in the specified preads file.')

        if write_contigs:
            for contig_fasta in contig_files:
                f = FastaReader(contig_fasta)
                for r in f:
                    rname = r.name.split()[0]
                    # Only output contigs which were added as paths.
                    if rname in self.paths:
                        fp_out.write('\t'.join([GFA_S_TAG, rname, r.sequence, 'LN:i:%d' % (len(r.sequence))]) + '\n');

        # Output links.
        for key, edge in self.edges.iteritems():
            link_line = self.format_gfa_v1_link_line(edge)
            fp_out.write(link_line + '\n')

        # Output contig paths.
        for ctg_name in sorted(self.paths.keys()):
            path_line = self.format_gfa_v1_path_line(ctg_name, self.paths[ctg_name], seq_len_map)
            fp_out.write(path_line + '\n')

    def add_read_from_node(self, node):
        """
        Takes a node formatted as "[0-9]+:[BE],
        extracts the read ID part and adds it to the GFAGraph.
        """
        r, rend = node.split(':')
        self.read_in_graph.add(r)

    def add_edge(self, v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_):
        """
        Adds an edge to the GFAGraph, but only if the edge
        between reads in v and w does not yet exist.
        GFA-1 format might not allow multiedges.
        """
        new_edge = [v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_]

        # GFA supposedly does not allow multiedges
        r1, r1end = v.split(':')
        r2, r2end = w.split(':')
        rp = [r1, r2]
        rp.sort()
        rp = tuple(rp)
        if rp in self.read_pairs:
            return
        self.read_pairs.add(rp)
        self.read_in_graph.add(r1)
        self.read_in_graph.add(r2)

        self.edges[(v, w)] = new_edge

    def update_edge(self, v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_):
        """
        If an edge between v and w already exists,
        this method will update only the values in the
        existing edge which are equal to None.
        This is useful if multiple sources of string graph
        edges are being loaded, and if they contain
        complementary info (e.g. tiling path and the sg.gexf
        graph).
        """
        if (v, w) not in self.edges:
            raise Exception('Edge "vwedge" does not exist and cannot be updated.'.format(vvedge=str((v, w))))

        # Update only non-None data.
        # Useful in case e.g. edges were loaded from the Networkx object which lacks some info,
        # while the tiling path might contain that info.
        # E.g. the tiling path specifies the contig where an edge came from, whereas
        # this info is not present in the string graph directly.
        # On the other hand, Unzip nx graph has info about phasing which is missing
        # from the tiling paths.
        curr_edge = self.edges[(v, w)]
        new_edge = [v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_]
        for i in xrange(len(curr_edge)):
            if curr_edge[i] == None:
                curr_edge[i] = new_edge[i]

    def add_or_update_edge(self, v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_):
        """
        A wrapper to add an edge if it does not exist
        or to update the edge if it already exists in
        the GFAGraph.
        """
        if (v, w) not in self.edges:
            self.add_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
        else:
            self.update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)

    def format_gfa_v1_link_line(self, edge):
        """
        Given an edge creates a GFA-1 formatted
        link ('L') string line.
        """
        v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_ = edge

        r1, r1end = v.split(':')
        r2, r2end = w.split(':')
        if r1end == 'E':
            o1 = GFA_ORIENT_FWD
        else:
            o1 = GFA_ORIENT_REV
        if r2end == 'E':
            o2 = GFA_ORIENT_FWD
        else:
            o2 = GFA_ORIENT_REV

        cig_str = '*' if not cigar else cigar

        vals = [GFA_L_TAG, r1, o1, r2, o2, cig_str]

        if overlap_len is not None:
            vals.append('ol:i:%d' % (int(overlap_len)))
        if overlap_idt is not None:
            vals.append('oi:f:%.1f' % (float(overlap_idt)))
        if overlap_begin is not None:
            vals.append('ob:i:%d' % (int(overlap_begin)))
        if overlap_end is not None:
            vals.append('oe:i:%d' % (int(overlap_end)))

        if src_graph is not None:
            vals.append('sg:Z:%s' % (str(src_graph)))
        if cross_phase is not None:
            vals.append('cp:Z:%s' % (str(cross_phase)))

        ctg_name_str = str(ctg_name) if ctg_name else 'NA'
        type_str = str(type_) if type_ else 'NA'
        vals.append('ci:Z:%s-%s' % (ctg_name_str, type_str))

        return '\t'.join(vals)

    def format_gfa_v1_path_line(self, ctg_name, tiling_path, seq_len_map):
        """
        Takes a tiling path and formulates a P line string for GFA-v1.
        seq_len_map is a dict with key equal to read header, and value
        the length of the read. If None, no CIGAR strings will be output
        in the path line.
        """
        if not tiling_path:
            return ''
        v, w, b, e, l, idt, etype = tiling_path[0]
        rname, rend = v.split(':')
        o = GFA_ORIENT_FWD if rend == 'E' else GFA_ORIENT_REV
        segs = [rname + o]
        segs_cigar = ['%dM' % (seq_len_map[rname])] if seq_len_map else ['*']
        for v, w, b, e, l, idt, etype in tiling_path:
            rname, rend = w.split(':')
            o = GFA_ORIENT_FWD if rend == 'E' else GFA_ORIENT_REV
            segs.append(rname + o)
            l = abs(b - e)
            if seq_len_map:
                segs_cigar.append('%dM' % (l))
            else:
                segs_cigar.append('*')
        out = [GFA_P_TAG, ctg_name, ','.join(segs), ','.join(segs_cigar)]
        return '\t'.join(out)
