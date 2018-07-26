class TilingPathEdge:

    def __init__(self, split_line = None):
        self.ctg_id, self.v, self.w, self.wrid, self.b, self.e, self.score, self.identity = \
                        None, None, None, None, None, None, None, None
        self.parsed = False
        if split_line:                          # pragma: no cover
            self.set_from(split_line)

    def set_from(self, split_line):
        assert(len(split_line) >= 8)
        self.parsed = False
        self.ctg_id = split_line[0]
        self.v = split_line[1]
        self.w = split_line[2]
        self.wrid = split_line[3]
        self.b = int(split_line[4])
        self.e = int(split_line[5])
        self.score = int(split_line[6])
        self.identity = float(split_line[7])
        self.parsed = True

    def get_split_line(self):
        return [str(val) for val in [self.ctg_id, self.v, self.w, self.wrid, self.b, self.e, self.score, self.identity]]

    # def __str__(self):
    #     return ' '.join(self.get_split_line())

class TilingPath:
    def __init__(self, tiling_edge_list, contig_sequence_len = None):
        self.edges = tiling_edge_list  # These are TilingPathEdge objects.
        self.v_to_edge = {}
        self.w_to_edge = {}
        self.coords = {}
        self.contig_len = 0
        self.first_node_offset = 0

        for i in xrange(1, len(tiling_edge_list)):
            assert(tiling_edge_list[i-1].w == tiling_edge_list[i].v)

        # If the total contig sequence len is known, use that to
        # calculate the length of the first read (in case proper
        # contigs are specified). This is needed to offset the coordinates
        # which can be calculated from the tiling path.
        if contig_sequence_len != None:
            _, tiling_len = calc_node_coords(tiling_edge_list)
            assert(contig_sequence_len >= tiling_len)
            self.first_node_offset = contig_sequence_len - tiling_len   # This is the length of the first node.

        # The self.coords is a dict: self.coords[v] = coordinate_on_contig
        self.coords, self.contig_len = calc_node_coords(tiling_edge_list, self.first_node_offset)

        # Sanity check.
        assert(contig_sequence_len == None or self.contig_len == contig_sequence_len)

        # Create a lookup of node to edge.
        self.v_to_edge = {}
        self.w_to_edge = {}
        for i in xrange(len(self.edges)):
            e = self.edges[i]
            self.v_to_edge[e.v] = i
            self.w_to_edge[e.w] = i

    def dump_as_split_lines(self):
        return [e.get_split_line() for e in self.edges]

    def get_subpath(self, start_coord, end_coord):
        """
        Given a TilingPath object, the method called `TilingPath.get_subpath() will
        attempt to extract a part of the tiling path which covers the specified
        coordinates.
        For example, user could specify alignment start and end positions, and
        provide the coordinates to this method, and the result would be a list
        of tiling path edges which correspond to the tiling in between the two
        coordinates (most likely slightly longer on both ends).

        Both start and end coordinates can be < 0 if the input contig was improper.

        Returns a list of split_line tiling path edges (not TilingEdge objects).
        """
        assert(self.edges)
        assert(start_coord <= end_coord)

        # end_coord -= 1  # Make the end inclusive.
        # start_node = None
        # end_node = None
        start_edge = None
        end_edge = None
        if start_coord < self.coords[self.edges[0].v]:
            start_edge = 0
        if end_coord <= self.coords[self.edges[0].v]:
            end_edge = 1
        for i in xrange(len(self.edges)):
            e = self.edges[i]
            if start_coord >= self.coords[e.v] and start_coord < self.coords[e.w]:
                start_edge = i
            if end_coord > self.coords[e.v] and end_coord <= self.coords[e.w]:
                end_edge = i + 1
        if end_coord >= self.coords[self.edges[-1].w]:
            end_edge = len(self.edges)
        assert(start_edge != None and end_edge != None)

        # Since the start_coord and end_coord can end within an edge,
        # we return the position in the final contigas.
        new_start_coord = start_coord - self.coords[self.edges[start_edge].v]
        new_end_coord = end_coord - self.coords[self.edges[start_edge].v]
        new_path = self.edges[start_edge:end_edge]
        new_path = [val.get_split_line() for val in new_path]
        return new_path, new_start_coord, new_end_coord

def calc_node_coords(tiling_edge_list, first_node_offset=0):
    """
    For a single tiling path (tiling_edge_list is a list
    of edges for a particular contig) calculates the
    genomic coordinate of every node in the path.
    In case there are cycles in the tiling path,
    the existing node's coordinate will be overwritten.
    `first_node_offset` refers to the length of the first node. If
    not specified, the contig length should not
    consider the length of the first node.
    """
    if not tiling_edge_list:
        return {}, 0
    coord_map = {}
    contig_len = 0
    edge0 = tiling_edge_list[0]
    coord_map[edge0.v] = first_node_offset
    for edge in tiling_edge_list:
        if edge.v not in coord_map:
            raise Exception(
                'Tiling path is not in sorted order. Node "{v!r}" does not yet have an assigned coordinate.'.format(v=edge.v))
        coord = coord_map[edge.v]
        coord += abs(int(edge.b) - int(edge.e))
        coord_map[edge.w] = coord
        contig_len = max(contig_len, coord)
    return coord_map, contig_len

def yield_split_line(fp_in):
    for line in fp_in:     # Example row: "0 000000007:B 000000005:B 000000005 9 0 1980 99.95"
        line = line.strip()
        if len(line) == 0: continue
        sl = line.split()
        yield sl

def load_tiling_paths(tp_file, contig_lens=None, whitelist_seqs=None):
    with open(tp_file) as fp_in:
        ret = load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)
    return ret

def load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None):
    split_lines = list(yield_split_line(fp_in))
    return load_tiling_paths_from_split_lines(split_lines, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

def load_tiling_paths_from_split_lines(split_lines, contig_lens=None, whitelist_seqs=None):
    """
    Parameters:
        contig_lens -   if a dict with contig sequence lengths is specified, the difference between the
                        contig len and the length obtained from the tiling path will be used to offset
                        the tiling path coordinates.
        whitelist_seqs - a dict or a set object containing contig IDs to load. If None, no filter will be applied.
    """
    tiling_path_edges = {}
    for sl in split_lines:     # Example row: "0 000000007:B 000000005:B 000000005 9 0 1980 99.95"
        new_edge = TilingPathEdge(sl)
        ctg_id = new_edge.ctg_id
        if whitelist_seqs != None and (ctg_id in whitelist_seqs) == False:
            continue
        tiling_path_edges.setdefault(ctg_id, [])
        tiling_path_edges[ctg_id].append(new_edge)

    # Convert the flat lists to TilingPath objects.
    # These keep track of
    tiling_paths = {}
    for ctg_id, edges in tiling_path_edges.iteritems():
        ctg_len = None
        if contig_lens != None and ctg_id in contig_lens:
            ctg_len = contig_lens[ctg_id]
        tiling_paths[ctg_id] = TilingPath(edges, ctg_len)

    return tiling_paths

def find_a_ctg_placement(p_paths, a_paths):
    """
    Determines placement coordinates for each a_ctg in a given a_paths dict of
    TilingPaths, and returns a dict of:
        placement[p_ctg_id][a_ctg_id] = (start, end, p_ctg_id, a_ctg_id, first_node, last_node)
    """
    placement = {}
    for a_ctg_id, a_tp in a_paths.iteritems():
        if len(a_tp.edges) == 0: continue       # pragma: no cover
        first_node = a_tp.edges[0].v
        last_node = a_tp.edges[-1].w
        p_ctg_id = a_ctg_id.split('-')[0].split('_')[0]
        p_tp = p_paths[p_ctg_id]
        start, end = p_tp.coords[first_node], p_tp.coords[last_node]
        placement.setdefault(p_ctg_id, {})
        placement[p_ctg_id][a_ctg_id] = (start, end, p_ctg_id, a_ctg_id, first_node, last_node)
    return placement
