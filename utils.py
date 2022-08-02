import os
import sys

import os
import sys
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)

import pprint
import copy
import math
import itertools
from collections import defaultdict
from pprint import pformat
import numpy as numpy
import numpy.linalg
from itertools import chain, product, combinations, combinations_with_replacement, permutations
import multiprocessing
import signal
import time
from sympy import Matrix, S, nsimplify, eye

try:
    import cvxpy
except:
    print("Error: could not import package cvxpy necessary for building the fixed deformation. Make sure it is installed.")


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


#############################################################################################################
# Topology Generator
#############################################################################################################

class TopologyGenerator(object):

    # Do not spam import warnings
    _HAS_ISSUED_IMPORT_WARNING = False

#   Below is tricky to make work
#    def __new__(cls, *args, **opts):
#        """ Factory creating an instance of the appropriate class for the inputs supplied for building the
#        loop topology. """
#        
#        if cls is TopologyGenerator:
#            if len(args)<1:
#                raise BaseException("The constructor of a TopologyGenerator instance requires at least one argumemt.")
#            if isinstance(args[0], dict) and all(isinstance(t,Propagator) for ts in args[0].values() for t in ts):
#                return super(TopologyGenerator, cls).__new__(TopologyGeneratorFromPropagators, *args, **opts)
#            else:
#                return super(TopologyGenerator, cls).__new__(cls, *args, **opts)
#        else:
#            return super(TopologyGenerator, cls).__new__(cls, *args, **opts)


    def __init__(self, edge_map_lin, powers=None):
        #if len(set(e[0] for e  in edge_map_lin)) != len(edge_map_lin):
        #    raise AssertionError("Every edge must have a unique name. Input: ", edge_map_lin)
        self.edge_map_lin = edge_map_lin
        self.edge_name_map = {name: i for (
            i, (name, _, _)) in enumerate(edge_map_lin)}
        self.edges = [(v1, v2) for (name, v1, v2) in edge_map_lin]
        vertices = [y for x in self.edges for y in x]

        self.num_vertices = len(set(vertices))

        self.ext = [i for i, x in enumerate(self.edges) if vertices.count(
            x[0]) == 1 or vertices.count(x[1]) == 1]
        self.vertices = vertices # with repetition
        # set the power for every propagator
        self.powers = {e[0]: powers[e[0]] if powers is not None and e[0] in powers
                else 1 for i, e in enumerate(edge_map_lin) if i not in self.ext}
        self.loop_momenta = None
        self.propagators = None

        self.n_loops = None
        self.spanning_trees = None

        # Not crucial, just a string representation of the denominators
        self.denominators_string = None

    def loop_momentum_bases(self):
        trees = []
        seen_states = set()
        self.generate_spanning_trees(trees, tree={self.edges[0][0]}, seen_state=seen_states)
        self.spanning_trees = trees
        self.n_loops = len(self.edge_map_lin) - len(trees[0])
        return [[i for i in range(len(self.edge_map_lin)) if i not in tree] for tree in trees]

    def merge_duplicate_edges(self, loop_momenta_names=None):
        """Create a new topology, merging edges with the same or opposite signature.
        If it is opposite, one of the lines will be flipped and the power will be increased.
        """
        if self.propagators is None:
            self.generate_momentum_flow(loop_momenta_names)

        duplicates = {}
        new_powers = {}
        if loop_momenta_names:
            duplicates = {tuple(self.propagators[self.edge_name_map[l]]): l for l in loop_momenta_names}
            new_powers = {e: p for e, p in self.powers.items() if e in loop_momenta_names}

        fuse_verts = set()
        new_edge_map_lin = []
        for p, e in zip(self.propagators, self.edge_map_lin):
            min_p = tuple([(x[0], not x[1]) for x in p])

            if p != () and min_p in duplicates:
                print('Flipping edge {} in deduplication process {} -> {}'.format(e[0], e[0], duplicates[min_p]))
                new_powers[duplicates[min_p]] += 1
                fuse_verts.add((e[1], e[2]))
            elif p != () and tuple(p) in duplicates and duplicates[tuple(p)] != e[0]:
                # an edge is a duplicate if it is internal with repeated signatures and it is not one of the loop momenta
                new_powers[duplicates[tuple(p)]] += 1
                fuse_verts.add((e[1], e[2]))
            else:
                new_edge_map_lin.append(e)
                if p != () and tuple(p) not in duplicates:
                    duplicates[tuple(p)] = e[0]
                    new_powers[e[0]] = self.powers[e[0]]

        for v1, v2 in fuse_verts:
            new_edge_map_lin = [(name, o1 if o1 != v2 else v1, o2 if o2 != v2 else v1) for name, o1, o2 in new_edge_map_lin]

        return TopologyGenerator(new_edge_map_lin, powers=new_powers)

    def generate_spanning_trees(self, result, tree, accum=[], excluded_edges=[], seen_state=None):
        """Compute all spanning trees of a graph component. Disconnected graphs
        are supported: only the component connected to the vertex in `tree` is considered.
        """
        if seen_state is not None:
            s = tuple(sorted(accum))
            if s in seen_state:
                return
            seen_state.add(s)

        # find all edges that connect the tree to a new node
        edges = [(i, e) for i, e in enumerate(
            self.edges) if e[0] != e[1] and len(set(e) & tree) == 1]
        edges_filtered = [(i, e) for i, e in edges if i not in excluded_edges]

        if len(edges_filtered) == 0:
            if len(edges) == 0:
                # no more new edges, so we are done
                s = list(sorted(accum))
                if s not in result:
                    result.append(s)
        else:
            # these edges can be added in any order, so only allow an order from lowest to highest
            for ei, (i, e) in enumerate(edges):
                excluded_edges.extend(j for j, _ in edges[:ei])
                new_v = e[0] if e[1] in tree else e[1]
                accum.append(i)
                tree.add(new_v)
                self.generate_spanning_trees(result, tree, accum, excluded_edges, seen_state)
                tree.remove(new_v)
                accum.pop()
                excluded_edges = excluded_edges[:-ei]

    def generate_spanning_tree(self):
        if len(self.edge_map_lin) == 0:
            return ()

        unprocessed_edges = list(enumerate(self.edges))
        vertices = {self.edges[0][0]}
        edges = []

        i = 0
        while len(unprocessed_edges) > 0:
            if i >= len(unprocessed_edges):
                i = 0

            e = unprocessed_edges[i][1]
            v1_in = e[0] in vertices
            v2_in = e[1] in vertices
            if v1_in and not v2_in or not v1_in and v2_in:
                vertices.add(e[0])
                vertices.add(e[1])
                edges.append(unprocessed_edges[i][0])
                unprocessed_edges.pop(i)
            elif v1_in and v2_in:
                unprocessed_edges.pop(i)

            i += 1

        return tuple(edges)


    def find_path(self, start, dest, excluding=set()):
        # find all paths from source to dest
        loop = start == dest
        start_check = 1 if loop else 0
        paths = [[(start, True)]]  # store direction
        if not loop:
            paths.append([(start, False)])
        res = []
        while True:
            newpaths = []
            for p in paths:
                if len(p) > 1 and p[-1][0] == dest:
                    res.append(p)
                    continue
                last_vertex = self.edges[p[-1][0]][1] if p[-1][1] else self.edges[p[-1][0]][0]
                for i, x in enumerate(self.edges):
                    if i not in excluding and all(i != pp[0] for pp in p[start_check:]):
                        if loop and i == start:
                            # if we need a loop, we need to enter from the right direction
                            if x[0] == last_vertex:
                                newpaths.append(p + [(i, True)])
                            continue

                        if x[0] == last_vertex and all(x[1] not in self.edges[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, True)])
                        if x[1] == last_vertex and all(x[0] not in self.edges[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, False)])
            paths = newpaths
            if len(paths) == 0:
                break
        return res

    def contruct_uv_forest(self, vertex_weights, edge_weights, UV_min_dod_to_subtract=0):
        """Construct the UV forest. The subgraphs in each spinney are order from smallest to largest, as this
        is the order in which they need to be Taylor expanded.
        """
        # construct all paths in a graph
        # we store edges instead of vertices, so that the procedure
        # supports duplicate edges
        cycles = set()
        unique_cycles = []
        for i, _ in enumerate(self.edge_map_lin):
            for path in self.find_path(i, i):
                path_rep = tuple(sorted(e for e, _ in path[1:]))
                if path_rep not in cycles:
                    unique_cycles.append([e for e, _ in path[1:]])
                    cycles.add(path_rep)

        trees = set() # all trees, also no subgraphs
        subgraphs = set()

        # now construct all subgraphs by taking the power set of cycles
        for subset in chain.from_iterable(combinations(unique_cycles, r) for r in range(len(unique_cycles) + 1)):
            # check if the set is vertex-wise connected
            vertex_subset = [set(v for e in cycle for v in self.edge_map_lin[e][1:]) for cycle in subset]

            subgraph_edges = tuple(sorted(set(e for c in subset for e in c)))
            if subgraph_edges in trees:
                continue
            
            trees.add(subgraph_edges)

            tree = set()
            while len(vertex_subset) > 0:
                for cycle in vertex_subset:
                    if len(tree) == 0 or len(tree & cycle) > 0:
                        tree |= cycle
                        vertex_subset.remove(cycle)
                        break
                else:
                    # not connected
                    break

            if len(vertex_subset) == 0:
                subgraphs.add(subgraph_edges)

        # filter for UV divergences
        div_subgraphs = []
        for s in subgraphs:
            dod = 0
            for e in s:
                dod += edge_weights[self.edge_map_lin[e][0]]

            vertices = set(v for e in s for v in self.edge_map_lin[e][1:])
            for v in vertices:
                dod += vertex_weights[v]

            loops = 0 if len(s) == 0 else len(s) - len(vertices) + 1
            dod += 4 * loops


            if dod >= UV_min_dod_to_subtract and loops > 0:
                div_subgraphs.append((s, dod))

        # create the forest from all UV subdivergences
        forest = []
        for subset in chain.from_iterable(combinations(div_subgraphs, r) for r in range(len(div_subgraphs) + 1)):
            vertex_subset = [set(v for e in cycle for v in self.edge_map_lin[e][1:]) for cycle, _ in subset]

            # check if the components are not overlapping (ie, they do not share a vertex or they are embedded in the other)
            for ((e1, _), v1), ((e2, _), v2) in combinations(zip(subset, vertex_subset), 2):
                fuse_len_v12 = len(v1 | v2)
                fuse_len_e12 = len(set(e1) | set(e2))
                if fuse_len_v12 != len(v1) and fuse_len_v12 != len(v2) and fuse_len_v12 != len(v1) + len(v2) or \
                    fuse_len_e12 != len(e1) and fuse_len_e12 != len(e2) and fuse_len_e12 != len(e1) + len(e2):
                    break
            else:
                forest.append(tuple(sorted(((tuple(sorted(self.edge_map_lin[e][0] for e in s)), dod) for s, dod in subset), key=lambda x: len(x[0]))))

        return forest

    def construct_uv_limits(self, vertex_weights, edge_weights, particle_ids=None, UV_min_dod_to_subtract=0, sg_name='N/A'):
        """Construct the uv limits by factorizing out all UV subgraphs from smallest to largest"""
        f = self.contruct_uv_forest(vertex_weights, edge_weights, UV_min_dod_to_subtract=UV_min_dod_to_subtract)
        #print('subgraph forest', f)

        new_graphs = []
        for spinney in f:
            gs = [{
                'uv_subgraphs': [],
                'uv_vertices' : [],
                'remaining_graph': copy.deepcopy(self),
                'spinney': spinney
                }]

            for graph_index, (subgraph_momenta_full, dod) in enumerate(spinney):
                # strip all subgraphs in the spinney from the current subgraph momenta
                # as these subgraphs have already been factored out by the Taylor expansion
                subgraph_momenta = [m for m in subgraph_momenta_full if not any(m in s for s, _ in spinney[:graph_index])]
                # strip indices of sub-sub-graphs
                subgraph_indices = [i for i, (m, _) in enumerate(spinney[:graph_index]) if len(set(m) & set(subgraph_momenta_full)) != 0 and
                    not any(len(m1) > len(m) and set(m).issubset(m1) for m1, _ in spinney[:graph_index])]

                new_gs = []
                for graph_info in gs:
                    g = graph_info['remaining_graph']
                    subgraph_vertices = set(v for m in subgraph_momenta for v in g.edge_map_lin[g.edge_name_map[m]][1:])
                    # find all external edges, supporting duplicate edges with one edge inside the subgraph and the other outside
                    subgraph_external_edges = [e[0] for i, e in enumerate(g.edge_map_lin) if len(set(e[1:]) & subgraph_vertices) > 0 and e[0] not in subgraph_momenta]
                    external_vertices = [v for v in subgraph_vertices if any(v in g.edge_map_lin[g.edge_name_map[e]][1:] for e in subgraph_external_edges)]

                    uv_subgraph_edges = [ee for ee in g.edge_map_lin if ee[0] in subgraph_momenta]
                    highest_vertex = max(v for e in uv_subgraph_edges for v in e[1:]) + 1
                    for c in subgraph_external_edges:
                        (_, cut_v1, cut_v2) = next(e for e in g.edge_map_lin if e[0] == c)
                        # add all external momenta to the UV subgraph
                        # for some configurations an external edge connects with the UV subgraph
                        # on both ends. In this case, we repeat the name
                        # TODO: check if this does not have unwanted side-effects
                        if cut_v1 in subgraph_vertices:
                            uv_subgraph_edges.append((c, cut_v1, highest_vertex))
                            highest_vertex += 1

                        if cut_v2 in subgraph_vertices:
                            uv_subgraph_edges.append((c, highest_vertex, cut_v2))
                            highest_vertex += 1

                    uv_subgraph = TopologyGenerator(uv_subgraph_edges, powers=g.powers)
                    uv_subgraph.inherit_loop_momentum_basis(self)
                    uv_subgraph.loop_momentum_bases()

                    # construct the remaining graph by shrinking the UV subgraph in g
                    remaining_graph_edges = [e for e in g.edge_map_lin if e[0] not in subgraph_momenta]
                    for ei, e in enumerate(remaining_graph_edges):
                        if e[1] in external_vertices[1:]:
                            remaining_graph_edges[ei] = (e[0], external_vertices[0], e[2])
                        if e[2] in external_vertices[1:]:
                            remaining_graph_edges[ei] = (e[0], remaining_graph_edges[ei][1], external_vertices[0])

                    # make sure the remaining diagram has the correct loop momentum basis
                    remaining_graph = TopologyGenerator(remaining_graph_edges, powers=g.powers)
                    remaining_graph.inherit_loop_momentum_basis(self)
                    
                    subgraph_loop_edges = [m for m in subgraph_momenta if len(uv_subgraph.propagators[uv_subgraph.edge_name_map[m]]) == 1 and 
                        uv_subgraph.propagators[uv_subgraph.edge_name_map[m]][0][0] not in uv_subgraph.ext]

                    # group subgraph momenta by their subgraph loop momentum signature
                    loop_lines = defaultdict(list)
                    for m in subgraph_momenta:
                        sig = tuple(uv_subgraph.edge_map_lin[k[0]][0] for k in uv_subgraph.propagators[uv_subgraph.edge_name_map[m]] 
                            if uv_subgraph.edge_map_lin[k[0]][0] in subgraph_loop_edges)
                        assert(len(sig) > 0)
                        loop_lines[sig].append(m)
                    loop_lines = [p for l, p in loop_lines.items()]

                    uv_vertices = []
                    n_loops = uv_subgraph.n_loops
                    for v, l in graph_info['uv_vertices']:
                        if v in subgraph_vertices:
                            n_loops += l
                        else:
                            uv_vertices.append((v, l))
                    uv_vertices.append((external_vertices[0], n_loops))

                    graph_configuration = {
                        'graph_index': graph_index,
                        'subgraph_indices': subgraph_indices,
                        'taylor_order': dod - UV_min_dod_to_subtract, # maximum taylor order
                        'external_edges': subgraph_external_edges,
                        'loop_edges': subgraph_loop_edges,
                        'graph': copy.deepcopy(uv_subgraph),
                        'subgraph_momenta': subgraph_momenta,
                        'full_subgraph_momenta': subgraph_momenta_full,
                    }

                    new_gs.append( {
                        'uv_subgraphs': [copy.deepcopy(x) for x in graph_info['uv_subgraphs']] + [graph_configuration],
                        'uv_vertices' : uv_vertices,
                        'remaining_graph': copy.deepcopy(remaining_graph),
                        'spinney': spinney,
                    })

                gs = new_gs
            new_graphs.extend(gs)

        #import pprint
        #pprint.pprint(new_graphs)
        return new_graphs

    def get_signature_map(self):
        # TODO: the order of self.ext should be consistent with Rust
        edge_map = {}
        for i, (prop, (edge_name, _, _)) in enumerate(zip(self.propagators, self.edge_map_lin)):
            signature = [[0]*len(self.loop_momenta), [0]*len(self.ext)]

            # For tadpole it may be that i is not in self.ext
            if prop == () and i in self.ext:
                signature[1][self.ext.index(i)] = 1 # FIXME: is it always +?

            for (mom, sign) in prop:
                s = 1 if sign else -1
                if mom not in self.ext:
                    signature[0][self.loop_momenta.index(mom)] = s
                else:
                    signature[1][self.ext.index(mom)] = s
            edge_map[edge_name] = signature
        return edge_map

    def inherit_loop_momentum_basis(self, supergraph, sink=None):
        lmbs = self.loop_momentum_bases()
        if self.n_loops > 0:
            # deterministicly select a basis based on the priority of the edges
            # this is based on whether they are part of the lmb and their name
            edge_names_sorted = list(sorted(supergraph.edge_map_lin[e][0] for e in supergraph.loop_momenta)) + \
                list(sorted(supergraph.edge_map_lin[e][0] for e in range(len(supergraph.edge_map_lin)) if e not in supergraph.loop_momenta))

            bases = [tuple(self.edge_map_lin[e][0] for e in b) for b in lmbs]
            bases = sorted(bases, key=lambda b: sum(2**edge_names_sorted.index(e) for e in b))
            self.generate_momentum_flow(loop_momenta=bases[0], sink=sink) 
        else:
            self.generate_momentum_flow(sink)

    def build_proto_topology(self, sub_graph, cuts, skip_shift=False):
        cut_momenta_names = [c['edge'] for c in cuts]

        sub_graph.inherit_loop_momentum_basis(self)

        # for each edge, determine the momentum map from full graph to subgraph and the shift
        # the shift will in general also have loop momentum dependence
        # TODO: the order of self.ext should be consistent with Rust
        edge_map = self.get_signature_map()
        sub_graph_edge_map = sub_graph.get_signature_map()
        loop_momentum_map = [[]]*sub_graph.n_loops
        param_shift = {}
        for (prop, (edge_name, _, _)) in zip(sub_graph.propagators, sub_graph.edge_map_lin):
            if prop == ():
                continue

            # determine the loop momentum map from the full loop momentum basis to the one of the subgraph
            if len(prop) == 1 and prop[0][0] not in sub_graph.ext:
                (mom, sign) = prop[0]
                s = 1 if sign else -1
                loop_momentum_map[sub_graph.loop_momenta.index(mom)] = [[s * a for a in i] for i in edge_map[edge_name]]

            # determine the shift in the subgraph in the cut basis
            signature = [[0]*len(cuts), [0]*len(self.ext)]

            if skip_shift:
                # keep the local subgraph map instead of translating it
                param_shift[edge_name] = ([0]*sub_graph.n_loops, sub_graph_edge_map[edge_name][1])
                continue

            # map the external momenta back to cuts and the external momenta of the full graph
            for ext_index, s in enumerate(sub_graph_edge_map[edge_name][1]):
                mom = sub_graph.edge_map_lin[sub_graph.ext[ext_index]][0]
                if s != 0:
                    if mom in cut_momenta_names:
                        signature[0][cut_momenta_names.index(mom)] = s
                    else:
                        # in the case of bubbles, we have propagators identical to a cut
                        # check if this is the case by comparing signatures
                        for mom1 in cut_momenta_names:
                            if edge_map[mom1] == edge_map[mom]:
                                signature[0][cut_momenta_names.index(mom1)] = s
                                break
                            if edge_map[mom1] == [[s * -1 for s in sig] for sig in edge_map[mom]]:
                                signature[0][cut_momenta_names.index(mom1)] = -s
                                break
                        else:
                            edge_index = next(i for i, e in enumerate(self.edge_map_lin) if e[0] == mom)
                            signature[1][self.ext.index(edge_index)] = s
            param_shift[edge_name] = signature

        return loop_momentum_map, param_shift

    def generate_momentum_flow(self, loop_momenta=None, sink=None):
        if loop_momenta is None:
            st = self.generate_spanning_tree()
            loop_momenta = [i for i in range(len(self.edge_map_lin)) if i not in st]
        else:
            self.n_loops = len(loop_momenta)
            if isinstance(loop_momenta[0], str):
                loop_momenta = [self.edge_name_map[edge_name]
                                for edge_name in loop_momenta]
            else:
                loop_momenta = loop_momenta

        self.loop_momenta = loop_momenta

        flows = []
        for l in loop_momenta:
            paths = self.find_path(l, l, excluding={lm for lm in loop_momenta if lm != l})
            assert(len(paths) == 1)
            flows.append(paths[0][:-1])

        # now route the external loop_momenta to the sink
        if sink is None:
            if len(self.ext) > 0:
                # vacuum bubbles don't have external momenta
                sink = self.ext[-1]
        else:
            sink = next(i for i, e in enumerate(self.edge_map_lin) if e[0] == sink)

        ext_flows = []
        for i, e in enumerate(self.ext):
            if e == sink:
                continue
            paths = self.find_path(e, sink, excluding=set(loop_momenta))
            assert(len(paths) == 1)
            ext_flows.append((i, paths[0]))

        # propagator momenta
        mom = []
        s = []
        for i, x in enumerate(self.edges):
            if i in self.ext:
                mom.append(())
                continue
            if i in loop_momenta:
                mom.append(((i, True),))
                #s.append("1/k{}^2".format(i))
                s.append("1/k{}^2".format(loop_momenta.index(i)))
            else:
                newmom = []
                s1 = "1/("
                for j, y in enumerate(flows):
                    for yy in y:
                        if yy[0] == i:
                            #s1 += ("+" if yy[1] else "-") + "k{}".format(loop_momenta[j])
                            s1 += ("+" if yy[1] else "-") + "k{}".format(j)
                            newmom.append((loop_momenta[j], yy[1]))
                            break
                for j, y in ext_flows:
                    overall_sign = 1 if y[0][1] else -1
                    for yy in y:
                        if yy[0] == i:
                            prop_sign = 1 if yy[1] else -1
                            s1 += ("+" if overall_sign * prop_sign == 1 else "-") + \
                                "{}".format(self.edge_map_lin[self.ext[j]][0])
                            newmom.append((self.ext[j], True if prop_sign * overall_sign == 1 else False))
                            break
                mom.append(tuple(newmom))
                s.append(s1 + ")^%d"%(2*self.powers[self.edge_map_lin[i][0]]))

        self.propagators = mom
        self.denominators_string = '(%s)'%('*'.join(s))
        #print("Constructed topology: {}".format('*'.join(s)))

    def create_loop_topology(self, name, ext_mom, mass_map={}, loop_momenta_names=None,
            contour_closure=None, analytic_result=None, fixed_deformation=None, constant_deformation=None,
            loop_momentum_map=None, cmb_indices=None, shift_map=None, numerator_tensor_coefficients=None,
            check_external_momenta_names=True):
        if loop_momentum_map is None:
            # FIXME: WHY?
            self.generate_momentum_flow(loop_momenta_names)

        # collect all loop lines and construct the signature
        loop_line_map = defaultdict(list)
        loop_line_vertex_map = defaultdict(list)

        # since edges could be flipped, we create an new shift map
        new_shift_map = copy.deepcopy(shift_map)
        powers = copy.deepcopy(self.powers)

        propagator_map = {}

        for prop, (edge_name, v1, v2) in zip(self.propagators, self.edge_map_lin):
            if prop == ():
                # external momentum
                continue

            mass = 0. if edge_name not in mass_map else mass_map[edge_name]

            # construct the signature
            signature = [0]*len(self.loop_momenta)
            q = [0., 0., 0., 0.]
            for (mom, sign) in prop:
                s = 1 if sign else -1
                if mom not in self.ext:
                    signature[self.loop_momenta.index(mom)] = s
                else:
                    q += numpy.array(ext_mom[self.edge_map_lin[mom][0]]) * s

            # we keep the direction of the loop momenta of the graph in the loop graph, so we need to flip
            # all edges that depend on one loop momentum and have negative sign
            should_flip = len([e for e in signature if e != 0]) == 1 and next(e == -1 for e in signature if e != 0)

            # flip the sign if the inverse exists
            # TODO: should this generate a minus sign when there is an odd-power numerator?
            alt_sig = tuple(s * -1 for s in signature)
            alt_prop = tuple((e, not d) for e, d in prop)
            if len(signature) > 0 and any(s != 0 for s in signature) and alt_sig in loop_line_map or should_flip:
                #print('warning: changing sign of propagator %s: %s -> %s' % (edge_name, tuple(signature), alt_sig) )
                if tuple(alt_sig) in loop_line_map:
                    for ll in loop_line_map[tuple(alt_sig)]:
                        # if the edge is duplicate, raise the power and don't add it to the map
                        if ll[3] == alt_prop and ll[2] == mass:
                            #print('Merging with flip', name, ll[0], edge_name)
                            propagator_map[edge_name] = ll[0]
                            powers[ll[0]] += powers[edge_name]
                            break
                    else:
                        loop_line_map[alt_sig].append((edge_name, [-a for a in q], mass, alt_prop))
                else:
                    loop_line_map[alt_sig].append((edge_name, [-a for a in q], mass, alt_prop))

                loop_line_vertex_map[alt_sig] += [(v2, v1)]

                if shift_map is not None:
                    new_shift_map[edge_name] = [[s * -1 for s in shift_map[edge_name][0]], [s * -1 for s in shift_map[edge_name][1]]]
            else:
                if tuple(signature) in loop_line_map:
                    for ll in loop_line_map[tuple(signature)]:
                        # if the edge is duplicate, raise the power and don't add it to the map
                        if ll[3] == prop and ll[2] == mass:
                            #print('Merging', name, ll[0], edge_name)
                            propagator_map[edge_name] = ll[0]
                            powers[ll[0]] += powers[edge_name]
                            break
                    else:
                        loop_line_map[tuple(signature)].append((edge_name, q, mass, prop))
                else:
                    loop_line_map[tuple(signature)].append((edge_name, q, mass, prop))

                loop_line_vertex_map[tuple(signature)] += [(v1, v2)]

        # vertices that are fused may again be fused with another vertex
        def multifuse(fuse_map, v):
            if v in fuse_map:
                return multifuse(fuse_map, fuse_map[v])
            else:
                return v

        # shrink all edges expect for the first per loop line
        for sig, edges in loop_line_vertex_map.items():
            # shrink all edges for signature 0 loop lines
            # this will fuse disjoint graphs for non-1PI graphs
            shrink_start = 0 if all(s == 0 for s in sig) else 1
            fuse_map = {e[0]: e[1] for e in edges[shrink_start:] if e[0] != e[1]}
            loop_line_vertex_map[sig] = [(edges[0][0], edges[0][1])]
            # duplicate keys are not allowed, unless the signature is 0
            # a duplicate key signals an orientation issue
            assert(all(s == 0 for s in sig) or len(fuse_map) == len(edges) - 1)
            
            # apply the fuse map to all loop lines
            for edges2 in loop_line_vertex_map.values():
                for i, edge in enumerate(edges2):
                    if edge[0] in fuse_map:
                        edges2[i] = (multifuse(fuse_map, edge[0]), edge[1])
                    if edge[1] in fuse_map:
                        edges2[i] = (edges2[i][0], multifuse(fuse_map, edge[1]))

        for sig, edges in loop_line_vertex_map.items():
            loop_line_vertex_map[sig] = edges[0]
        
        loop_line_list = list(loop_line_map.items())
        ll = [LoopLine(
            start_node=loop_line_vertex_map[signature][0],
            end_node=loop_line_vertex_map[signature][1],
            signature=signature,
            propagators=tuple(
                Propagator(q=q, m_squared=mass**2,
                            power=powers[edge_name],
                            parametric_shift=new_shift_map[edge_name] if shift_map is not None else None,
                            name=edge_name)
                for (edge_name, q, mass, _) in propagators)) for signature, propagators in loop_line_list]

        # the external kinematics are stored in the right order
        try:
            external_kinematics = [ext_mom["q%d"%n] for n in sorted([int(qi.replace("q","")) for qi in ext_mom.keys()])]
        except ValueError:
            if check_external_momenta_names:
                print("External kinematics are not labels as q1,...,qn. Their order will be random.")
            external_kinematics = list(ext_mom.values())
        
        loop_topology = LoopTopology(name=name, n_loops=len(self.loop_momenta), external_kinematics=external_kinematics,
            ltd_cut_structure=[], loop_lines=ll, analytic_result = analytic_result, fixed_deformation = fixed_deformation,
            constant_deformation = constant_deformation, loop_momentum_map=loop_momentum_map, cmb_indices=cmb_indices,
            propagator_map=propagator_map)

        if analytic_result is None:
            loop_topology.analytic_result = self.guess_analytical_result(self.loop_momenta, ext_mom, mass_map)

        loop_topology.numerator_tensor_coefficients = numerator_tensor_coefficients


        return loop_topology

class LoopTopology(object):
    """ A simple container for describing a loop topology."""

    _cvxpy_threshold = 1.0e-12
    # The existence condition itself is inaccurate because of the kinematics being not exactly on-shell
    # and squaring + squareroots decreases the accuracy down to less than 16/2 digits.
    _existence_threshold = 1.0e-7
    def __init__(self, ltd_cut_structure, loop_lines, external_kinematics, n_loops=1, name=None, analytic_result=None,
        fixed_deformation=None, constant_deformation=None, maximum_ratio_expansion_threshold=None, loop_momentum_map=None,
        cmb_indices=None, numerator_tensor_coefficients=None, propagator_map={}, **opts):
        """
            loop_lines          : A tuple of loop lines instances corresponding to each edge of the directed
                                  graph of this topology.
            ltd_cut_structure   : A tuple of tuples instructing what cuts result from the iterative
                                  application of LTD. Example: ((0,1,-1),(1,-1,0),...)
                                  Each entry of the tuple is either 0 (NO_CUT), -1 (NEGATIVE_CUT sol.)
                                  or +1 (POSTIVIVE_CUT sol.).
            n_loops             : Number of loops in this topology
            name                : Name of this topology
        """
        self.external_kinematics = external_kinematics
        self.loop_lines          = loop_lines
        self.ltd_cut_structure   = ltd_cut_structure
        self.n_loops             = n_loops 
        self.name                = name
        if callable(analytic_result):
            self.analytic_result   = analytic_result(self.external_kinematics)
        else:
            self.analytic_result   = analytic_result

        self.fixed_deformation = fixed_deformation
        self.constant_deformation = constant_deformation
        self.maximum_ratio_expansion_threshold = maximum_ratio_expansion_threshold
        self.loop_momentum_map = loop_momentum_map
        self.cmb_indices = cmb_indices
        self.numerator_tensor_coefficients = numerator_tensor_coefficients
        self.propagator_map = propagator_map

    def __str__(self):
        return pformat(self.to_flat_format())
       



class LoopLine(object):
    """ A simple container for describing a loop line."""

    NO_CUT          = 0
    POSITIVE_CUT    = 1
    NEGATIVE_CUT    = -1

    def __init__(self, signature, propagators, start_node, end_node, **opts):
        """
            signature           : The signature of the loop line specifies its dependence on the loop
                                  momenta. It is a tuple of length n_loops and with entries which are
                                  either 0 (no dependence), +1 (positive dependence) or -1 (negative dependence).
                                  For instance (1,0,-1) means a loop momenta dependency reading 1*k + 0*l -1*m.
            propagators         : List of Propagator instances specifying all propagators making up this loop line.
            start_node          : Integer labeling the starting node of this directed edge.
            end_node            : Integer labeling the end node of this directed edge (where the arrow points).
        """
        self.signature      = signature
        self.propagators    = propagators
        # Automatically forward the signature attribute to the loop propagators if not specified by the user.
        for propagator in self.propagators:
            propagator.set_signature(self.signature, force=False)
        self.start_node     = start_node
        self.end_node       = end_node

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse of this loop line with the provided list loop momenta, given as a list of LorentzVector."""

        result = 1.0
        for prop in self.propagators:
            result *= prop.evaluate_inverse(loop_momenta)
        return result

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return LoopLine(
            signature   =   tuple(flat_dict['signature']),
            propagators =   tuple([Propagator.from_flat_format(p) for p in flat_dict['propagators']]),
            start_node  =   flat_dict['start_node'],
            end_node    =   flat_dict['end_node']             
        ) 


    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['signature'] = list(self.signature)
        res['propagators'] = [p.to_flat_format() for p in self.propagators]
        res['start_node'] = self.start_node
        res['end_node'] = self.end_node

        return res

class Propagator(object):
    """ A simple container for describing a loop propagator."""

    def __init__(self, q, m_squared, power=1, signature = None, name = None, parametric_shift = None, uv = False, **opts):
        self.name       = name
        self.q          = q
        self.m_squared  = m_squared
        self.power      = power
        self.uv         = uv
        # Note that this signature member is not strictly speaking necessary as it should always be the
        # same as the signature attribute of the LoopLine containing it. We forward it here for convenience however.
        self.signature  = signature
        self.parametric_shift = parametric_shift
        self.name = name

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        res={}

        res['q'] = [float(v) for v in self.q]
        res['m_squared'] = self.m_squared
        res['power'] = self.power
        res['uv'] = self.uv

        if self.parametric_shift is not None:
            res['parametric_shift'] = self.parametric_shift

        if self.name is not None:
            res['name'] = self.name

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return Propagator(
            q           =   vectors.LorentzVector(flat_dict['q']),
            m_squared   =   flat_dict['m_squared'],
            power       =   flat_dict['power'],
            name        =   flat_dict['name'],
            uv          =   flat_dict['uv'],
            parametric_shift = flat_dict['parametric_shift'] if 'parametric_shift' in flat_dict else None
        ) 

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse propagator with the provided list loop momenta, given as a list of LorentzVector.""" 
        return (sum(wgt*loop_momentum for wgt, loop_momentum in zip(self.signature, loop_momenta)) + self.q).square() - self.m_squared

    def set_signature(self, signature, force=True):
        """ Overwrite (if force=True) the signature attribute of this propagator."""
        if (self.signature is None) or force:
            self.signature = signature