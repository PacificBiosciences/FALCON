from __future__ import absolute_import
import os
import sys
import json

GFA_H_TAG = 'H'
GFA_S_TAG = 'S'
GFA_L_TAG = 'L'
GFA_P_TAG = 'P'
GFA_ORIENT_FWD = '+'
GFA_ORIENT_REV = '-'
GFA_SEQ_UNKNOWN = '*'
GFA_LINK_CIGAR_UNKNOWN = '*'
GFA2_E_TAG = 'E'

KW_NAME = 'name'
KW_TAGS = 'tags'
KW_LABELS = 'labels'

KW_NODE_SEQ = 'seq'
KW_NODE_LEN = 'len'

KW_EDGE_SOURCE = 'v'
KW_EDGE_SOURCE_ORIENT = 'v_orient'
KW_EDGE_SINK = 'w'
KW_EDGE_SINK_ORIENT = 'w_orient'
KW_EDGE_CIGAR = 'cigar'
KW_EDGE_SOURCE_START = 'v_start'
KW_EDGE_SOURCE_END = 'v_end'
KW_EDGE_SINK_START = 'w_start'
KW_EDGE_SINK_END = 'w_end'

KW_PATH_NODES = 'nodes'
KW_PATH_CIGARS = 'cigars'

"""
GFA-1:
- H line: line = '\t'.join([GFA_H_TAG, '\tVN:Z:1.0'])
- S line: line = '\t'.join([GFA_S_TAG, rname, GFA_SEQ_UNKNOWN if (not write_reads) else r.sequence, 'LN:i:%s' % len(r.sequence)])
- L line: line = '\t'.join([GFA_L_TAG, edge.sg_edge.v_name, edge.sg_edge.v_orient, edge.sg_edge.w_name, edge.sg_edge.w_orient, cig_str])
- P line: line = '\t'.join([GFA_P_TAG, ctg_name, ','.join(segs), ','.join(segs_cigar)]

GFA-2:
- H line: line = '\t'.join([GFA_H_TAG, '\tVN:Z:2.0'])
- S line: line = '\t'.join([GFA_S_TAG, rname, str(len(r.sequence)), GFA_SEQ_UNKNOWN if (not write_reads) else r.sequence])
- E line: line = '\t'.join([GFA2_E_TAG, edge_name, source_node, sink_node, source_start, source_end, sink_start, sink_end, cig_str])

"""

class GFAGraph:
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.paths = {}

        """
        Node: {KW_NAME: '01234', KW_NODE_SEQ: 'ACTG', 'len': 4}
        Node: {'name': '56789', KW_NODE_SEQ: 'CAGT', 'len': 4}
        Edge: {KW_NAME: 'edge1', 'source': '01234', 'sink': '56789', 'cigar': '*', 'source_start': 3, 'source_end': 4, 'sink_start': 0, 'sink_end': 1}
        Path: {KW_NAME: '000000F', 'nodes': ['01234', '56789'], '
        """

    def add_node(self, node_name, node_len, node_seq='*', tags={}, labels={}):
        if len(node_name) == 0:
            raise 'Node name should be a non-empty string.\n'
        if node_len < 0:
            raise 'Node length should be >= 0.\n'
        if len(node_seq) == 0:
            raise 'Node sequence should be a non-empty string. Use "*" instead.\n'
        if isinstance(tags, dict) == False:
            raise 'The tags object must be a dict.\n'
        if isinstance(labels, dict) == False:
            raise 'The labels object must be a dict.\n'

        self.nodes[node_name] = {
                                    KW_NAME: node_name,
                                    KW_NODE_LEN: node_len,
                                    KW_NODE_SEQ: node_seq,
                                    KW_TAGS: tags,
                                    KW_LABELS: labels
                                }

    def add_edge(self, edge_name, source, source_orient, sink, sink_orient, source_start, source_end, sink_start, sink_end, cigar, tags={}, labels={}):
        """
        source_orient   + if fwd, - otherwise.
        sink_orient   + if fwd, - otherwise.
        """
        if len(edge_name) == 0:
            raise 'Edge name should be a non-empty string.\n'
        if len(source) == 0:
            raise 'Source node not specified.\n'
        if len(sink) == 0:
            raise 'Sink node not specified.\n'
        if source_orient not in '+-':
            raise 'Source orientation should be either "+" or "-".\n'
        if sink_orient not in '+-':
            raise 'Sink orientation should be either "+" or "-".\n'
        if source_start < 0 or source_end < 0:
            raise 'Source coordinates should be >= 0.\n'
        if sink_start < 0 or sink_end < 0:
            raise 'Sink coordinates should be >= 0.\n'
        if len(cigar) == 0:
            raise 'Cigar string should not be empty. Use "*" instead.\n'
        if source_end < source_start:
            sys.stderr.write('ERROR with: source = %s, source_start = %s, source_end = %s, sink = %s, sink_start = %s, sink_end = %s\n' % (source, source_start, source_end, sink, sink_start, sink_end))
            raise 'Source end coordinate should be >= source start coordinate.\n'
        if sink_end < sink_start:
            raise 'Sink end coordinate should be >= sink start coordinate.\n'
        if isinstance(tags, dict) == False:
            raise 'The tags object must be a dict.\n'
        if isinstance(labels, dict) == False:
            raise 'The labels object must be a dict.\n'

        self.edges[str((source, sink))] = {
                                        KW_NAME: edge_name,
                                        KW_EDGE_SOURCE: source,
                                        KW_EDGE_SOURCE_ORIENT: source_orient,
                                        KW_EDGE_SINK: sink,
                                        KW_EDGE_SINK_ORIENT: sink_orient,
                                        KW_EDGE_SOURCE_START: source_start,
                                        KW_EDGE_SOURCE_END: source_end,
                                        KW_EDGE_SINK_START: sink_start,
                                        KW_EDGE_SINK_END: sink_end,
                                        KW_EDGE_CIGAR: cigar,
                                        KW_TAGS: tags,
                                        KW_LABELS: labels
                                    }

    def add_path(self, path_name, path_nodes, path_cigars, tags={}, labels={}):
        """
        path_nodes is a list of nodes which should be joined
        consecutively in a path.
        path_cigars is a list of CIGAR strings describing how the
        two neighboring nodes are joined.
        len(path_nodes) == len(path_cigars)
        """
        if len(path_name) == 0:
            raise 'Path name should be a non-empty string.\n'
        if len(path_nodes) == 0:
            raise 'Path nodes should be a non-empty list.\n'
        if len(path_cigars) == 0:
            raise 'Path cigars should be a non-empty list.\n'
        if isinstance(tags, dict) == False:
            raise 'The tags object must be a dict.\n'
        if isinstance(labels, dict) == False:
            raise 'The labels object must be a dict.\n'
        if len(path_nodes) != len(path_cigars):
            raise 'The path_nodes and path_cigars should have the same length.\n'

        self.paths[path_name] = {
                                    KW_NAME: path_name,
                                    KW_PATH_NODES: path_nodes,
                                    KW_PATH_CIGARS: path_cigars,
                                    KW_TAGS: tags,
                                    KW_LABELS: labels
                                }

    def write_gfa_v1(self, fp_out):
        # Header
        line = '\t'.join([GFA_H_TAG, 'VN:Z:1.0'])
        fp_out.write(line + '\n')

        # Sequences.
        for node_name, node_data in self.nodes.iteritems():
            line = '\t'.join([  GFA_S_TAG,
                                node_data[KW_NAME],
                                node_data[KW_NODE_SEQ],
                                'LN:i:%d' % node_data[KW_NODE_LEN]])
            fp_out.write(line + '\n')

        for edge, edge_data in self.edges.iteritems():
            cigar = edge_data[KW_EDGE_CIGAR] if edge_data[KW_EDGE_CIGAR] != '*' else '%dM' % (abs(edge_data[KW_EDGE_SINK_END] - edge_data[KW_EDGE_SINK_START]))

            line = '\t'.join([str(val) for val in
                                [  GFA_L_TAG,
                                    edge_data[KW_EDGE_SOURCE],
                                    edge_data[KW_EDGE_SOURCE_ORIENT],
                                    edge_data[KW_EDGE_SINK],
                                    edge_data[KW_EDGE_SINK_ORIENT],
                                    cigar
                                ]
                            ])
            fp_out.write(line + '\n')

        for path_name, path_data in self.paths.iteritems():
            line = '\t'.join([GFA_P_TAG, path_data[KW_NAME], ','.join(path_data[KW_PATH_NODES]), ','.join(path_data[KW_PATH_CIGARS])])
            fp_out.write(line + '\n')

    def write_gfa_v2(self, fp_out):
        # Header
        line = '\t'.join([GFA_H_TAG, 'VN:Z:2.0'])
        fp_out.write(line + '\n')

        # Sequences.
        for node_name, node_data in self.nodes.iteritems():
            line = '\t'.join([  GFA_S_TAG,
                                node_data[KW_NAME],
                                str(node_data[KW_NODE_LEN]),
                                node_data[KW_NODE_SEQ]])
            fp_out.write(line + '\n')

        for edge, edge_data in self.edges.iteritems():
            v = edge_data[KW_EDGE_SOURCE]
            w = edge_data[KW_EDGE_SINK]
            v_len = self.nodes[v][KW_NODE_LEN]
            w_len = self.nodes[w][KW_NODE_LEN]

            # GFA-2 specifies a special char '$' when a coordinate is the same as the sequence length.
            v_start = str(edge_data[KW_EDGE_SOURCE_START]) + ('$' if edge_data[KW_EDGE_SOURCE_START] == v_len else '')
            v_end = str(edge_data[KW_EDGE_SOURCE_END]) + ('$' if edge_data[KW_EDGE_SOURCE_END] == v_len else '')
            w_start = str(edge_data[KW_EDGE_SINK_START]) + ('$' if edge_data[KW_EDGE_SINK_START] == w_len else '')
            w_end = str(edge_data[KW_EDGE_SINK_END]) + ('$' if edge_data[KW_EDGE_SINK_END] == w_len else '')

            line = '\t'.join([str(val) for val in
                                [  GFA2_E_TAG, edge_data[KW_NAME],
                                    edge_data[KW_EDGE_SOURCE] + edge_data[KW_EDGE_SOURCE_ORIENT],
                                    edge_data[KW_EDGE_SINK] + edge_data[KW_EDGE_SINK_ORIENT],
                                    v_start, v_end,
                                    w_start, w_end,
                                    edge_data[KW_EDGE_CIGAR],
                                ]
                            ])
            fp_out.write(line + '\n')

def serialize_gfa(gfa_graph):
    gfa_dict = {}
    gfa_dict['nodes'] = gfa_graph.nodes
    gfa_dict['edges'] = gfa_graph.edges
    gfa_dict['paths'] = gfa_graph.paths
    return json.dumps(gfa_dict)

def deserialize_gfa(fp_in):
    gfa_dict = json.load(fp_in)
    gfa = GFAGraph()
    gfa.nodes = gfa_dict['nodes']
    gfa.edges = gfa_dict['edges']
    gfa.paths = gfa_dict['paths']
    return gfa
