#!/usr/bin/env python3

import copy
import logging
import os
from pathlib import Path
from pprint import pprint, pformat
import math
import time
import numpy as np
from itertools import combinations_with_replacement
from collections import OrderedDict
import glob as glob_module

import importlib
import progressbar
from itertools import chain, product
import sys
import functools
import subprocess
import argparse
import shutil
import py_compile
from warnings import catch_warnings
import os
import stat

from sympy import Matrix, diag

pjoin = os.path.join

import re
import multiprocessing

import utils

logger = logging.getLogger('rqft')

try:
    import yaml
    from yaml import Loader, Dumper
except ImportError:
    raise BaseException("Install yaml python module in order to import/export topologies from/to yaml.")

if __name__ == "__main__":
    logging.basicConfig()
    logger.setLevel(logging.INFO)

class ForestGenerator:
    def __init__(self, edges, name, incoming_momentum_names,
        loop_momenta_names=None, masses={}, powers=None, particle_ids={}, vertex_weights={}, edge_weights={}):

        self.name = name
        self.topo = utils.TopologyGenerator(edges, powers)
        self.topo.generate_momentum_flow(loop_momenta_names)

        vw = {v: vertex_weights[v] if v in vertex_weights else 0 for e in self.topo.edge_map_lin for v in e[1:]}
        ew = {e: edge_weights[e] if e in edge_weights else -2 for e, _, _ in self.topo.edge_map_lin}
        
        uv_limits = self.topo.construct_uv_limits(vw, ew, particle_ids=particle_ids, sg_name=self.name)
        # only keep the limit of the complete graph
        uv_limits = [uv_structure for uv_structure in uv_limits if uv_structure['remaining_graph'].n_loops == 0 and len(uv_structure['uv_subgraphs']) > 0]

        subgraph_counter = {}
        graph_counter = 0
        for uv_limit in uv_limits:
            # give every subdiagram a globally unique id
            for uv_sg in uv_limit['uv_subgraphs']:
                r = tuple(uv_sg['subgraph_momenta'])
                if r in subgraph_counter:
                    uv_sg['first_occurrence_id'] = subgraph_counter[r]
                else:
                    uv_sg['first_occurrence_id'] = graph_counter
                    subgraph_counter[r] = graph_counter

                uv_sg['id'] = graph_counter
                graph_counter += 1

        uv = [{
                'uv_subgraphs': uv_limit['uv_subgraphs'],
                'uv_spinney': [[list(g), dod] for g, dod in uv_limit['spinney']],
                'uv_vertices': [x for x in uv_limit['uv_vertices']],
                'uv_propagators': [m for g, _ in uv_limit['spinney'] for m in g],
            } for uv_limit in uv_limits]

        edge_map = self.topo.get_signature_map()

        uv_representative_graphs = {}

        for uv_structure in uv:
            # create the LTD representation of the derived UV graph
            forest_to_lmb = []
            for uv_subgraph in uv_structure['uv_subgraphs']:
                (loop_mom_map, shift_map) = self.topo.build_proto_topology(uv_subgraph['graph'], [], skip_shift=True)

                # remove all shifts from the loop momentum map as they will be expanded in
                basis_shift_map = [[[0]*len(lmp), shp] for lmp, shp in loop_mom_map]
                loop_mom_map = [[lmp, [0]*len(shp)] for lmp, shp in loop_mom_map]

                forest_to_lmb.extend([x[0] for x in loop_mom_map])

                # note: the shift map signs may get swapped when edges switch orientation
                # the basis_shift_map should be unaffected since defining edges are not swapped
                loop_topo = uv_subgraph['graph'].create_loop_topology(name,
                    # provide dummy external momenta
                    ext_mom={edge_name: [0, 0, 0, 0] for (edge_name, _, _) in self.topo.edge_map_lin},
                    fixed_deformation=False,
                    mass_map=masses,
                    loop_momentum_map=loop_mom_map,
                    shift_map=shift_map,
                    check_external_momenta_names=False,
                    analytic_result=0)
                loop_topo.external_kinematics = []

                # take the UV limit of the diagram, add the mass and set the parametric shift to 0
                # collect the parametric shifts of the loop lines such that it can be used to Taylor expand
                # the UV subgraph
                uv_loop_lines = []
                for ll in loop_topo.loop_lines:
                    derived_ll_power = sum(pp.power for pp in ll.propagators)
                    # determine the power of the loop line of the non-derived graph
                    orig_ll_power = sum(uv_subgraph['graph'].powers[pp.name]
                        + len([1 for _,mn in loop_topo.propagator_map.items() if mn == pp.name])
                        for pp in ll.propagators)
                    prop_pow = 1

                    # note: we assume that every propagator has power 1 and that only merging of identical edges raises it
                    uv_loop_lines.append((ll.signature, [(p.name, p.parametric_shift,
                        (prop_pow if ii == 0 else 1) + len([1 for _,mn in loop_topo.propagator_map.items() if mn == p.name])) for ii, p in enumerate(ll.propagators)], derived_ll_power - orig_ll_power))

                    prop = ll.propagators[0]
                    prop.uv = True
                    prop.m_squared = 100 # some UV mass, unused
                    prop.power = derived_ll_power + prop_pow - 1
                    prop.parametric_shift = [[], [0 for _ in range(len(incoming_momentum_names) * 2)]]
                    ll.propagators = [prop]

                loop_topo.uv_loop_lines = (uv_loop_lines, basis_shift_map)
                uv_subgraph['loop_topo'] = loop_topo

                # check if raised loop lines have a propagator with a shift, otherwise this configuration is impossible
                uv_subgraph['skip_pf'] = False
                for i, (ll_sig, propagators, raised_power) in enumerate(loop_topo.uv_loop_lines[0]):
                    if raised_power > 0 and all(all(s == 0 for s in param_shift[1]) for _, param_shift, _ in propagators):
                        if all(all(x == 0 for y in lm_shift for x in y) for s, lm_shift in zip(ll_sig, loop_topo.uv_loop_lines[1]) if s != 0):
                            uv_subgraph['skip_pf'] = True
                            break


            if forest_to_lmb != []:
                fmb_in_cb = Matrix(forest_to_lmb).tolist()
                loops = [r for r in fmb_in_cb]
                loopsi = Matrix(loops)**-1
                extshift = Matrix([[0]*len(self.topo.ext)]*(len(loops)))
                extshiti = loopsi * extshift

                uv_structure['forest_to_cb_matrix'] = (loopsi.tolist(), [[]]*len(forest_to_lmb), extshiti.tolist(), loops, [[]]*len(forest_to_lmb), extshift.tolist())
            else:
                uv_structure['forest_to_cb_matrix'] = ([[]], [[]], [[]], [[]], [[]], [[]])

        #pprint(uv)
        self.uv = uv



class FORMGraph(object):
    _include_momentum_routing_in_rendering=False
    _include_edge_name_in_rendering=True
    _rendering_size = (1.0*(11.0*60),1.0*(8.5*60)) # 1.0 prefactor should be about 1 landscape A4 format per graph
    # Choose graph layout strategy. Interesting options are in comment.
    _graph_layout_strategy = '{"PackingLayout"->"ClosestPacking"}' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringElectricalEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredEmbedding", "Orientation" -> Left, "RootVertex" -> "I1"}' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}'
    
    # '{"SpringEmbedding"}' gives interesting results too.

    def __init__(self, *args,
        call_identifier=None,
        name=None,
        edges=None,
        nodes=None,
        overall_factor="1",
        multiplicity=1,
        benchmark_result=0.0,
        default_kinematics=None,
        effective_vertex_id=None
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.is_zero = False
        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        self.multiplicity = multiplicity
        # Give the possibility of specifying a benchmark result to the output yaml file for future ref.
        self.benchmark_result = benchmark_result
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

        self.default_kinematics = default_kinematics
        self.effective_vertex_id = effective_vertex_id
        self.topology = None
        self.configurations = None

 
    def get_mathematica_rendering_code(self, FORM_id=None, lmb_id=None):
        """ Generate mathematica expression for drawing this graph."""

        repl_dict = {
            'edge_font_size'     : 10,
            'vertex_font_size'    : 10,
            'width'              : self._rendering_size[0],
            'height'             : self._rendering_size[1],
            'graph_layout_strategy' : self._graph_layout_strategy
        }
        if self.call_identifier and all(k in self.call_identifier for k in ['proc_id','left_diagram_id','right_diagram_id']):
            graph_name='MG: %s'%('P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%self.call_identifier)
        else:
            graph_name='MG: %s'%self.name
        if FORM_id is not None:
            graph_name +=' | FORM: #%d'%FORM_id
        if lmb_id is not None:
            graph_name +=' | LMB: #%d'%lmb_id
        repl_dict['graph_name'] = graph_name

        # Special name rendering rules
        def get_part_name(pdg):
            quark_names = {1:"d",2:"u",3:"s",4:"c",5:"b",6:"t"}
            if pdg==11:
                return r"\!\(\*SuperscriptBox[\(e\), \(+\)]\)"
            elif pdg==-11:
                return r"\!\(\*SuperscriptBox[\(e\), \(-\)]\)"
            elif pdg==22:
                return r"\[Gamma]"
            elif pdg in [-1,-2,-3,-4,-5,-6]:
                return r"\!\(\*OverscriptBox[\(%s\), \(_\)]\)"%quark_names[abs(pdg)]
            else:
                return {1: 'd', 21: 'g', 82: 'h'}[pdg]

        node_key_to_node_name = {}
        def get_node_label(node_key):
            if node_key in node_key_to_node_name:
                return node_key_to_node_name[node_key]
            node = self.nodes[node_key]
            if 'renormalisation_vertex_n_loops' in node:
                label = 'UV%dL'%node['renormalisation_vertex_n_loops']
                label += '%dP'%node['renormalisation_vertex_n_shrunk_edges']
                # If that label is not unique, then prefix it the actual key.
                if label in node_key_to_node_name.values():
                    label = '%s_%s'%(str(node_key),label)
            else:
                label = str(node_key)

            node_key_to_node_name[node_key] = label
            return label

        # Generate edge list
        edge_template = """Labeled[Style[CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]%(edge_style)s,%(edge_color)s,Thickness[%(thickness)f]],"%(edge_label)s"]"""
        all_edge_definitions = []
        all_edge_shape_functions = []
        for edge_key, edge_data in self.edges.items():
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key)
            edge_repl_dict = {}
            # is_LMB = ('name' in edge_data and 'LMB' in str(edge_data['name']).upper())
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            is_LMB = (sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0)
            if is_LMB:
                edge_repl_dict['thickness'] = 0.005
            else:
                edge_repl_dict['thickness'] = 0.002
            edge_repl_dict['in_node'] = get_node_label(edge_key[0])
            edge_repl_dict['out_node'] = get_node_label(edge_key[1])
            edge_repl_dict['edge_key'] = edge_key[2]
            edge_repl_dict['arrow_style'] = 'Arrow' if not is_LMB else 'HalfFilledDoubleArrow'
            edge_repl_dict['arrow_size'] = 0.015 if not is_LMB else 0.025
            all_edge_shape_functions.append(
                'CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]->GraphElementData["%(arrow_style)s", "ArrowSize" -> %(arrow_size)f]'%edge_repl_dict
            )
            if 'name' in edge_data and 'CUT' in str(edge_data['name']).upper():
                edge_repl_dict['edge_style'] = ",Dashed"
            else:
                edge_repl_dict['edge_style'] = ""
            if edge_data['PDG'] in [-1,-2,-3,-4,-5,1,2,3,4,5]:
                color = "Cyan"
            elif edge_data['PDG'] in [-6,6]:
                color = "Blue"
            elif edge_data['PDG'] in [21,]:
                color = "Red"
            elif edge_data['PDG'] in [82,-82]:
                color = "Pink"
            elif edge_data['PDG'] in [25,]:
                color = "Green"
            else:
                color = "Gray"
            edge_repl_dict['edge_color'] = color
            edge_label_pieces = [get_part_name(edge_data['PDG']),]
            if 'name' in edge_data and self._include_edge_name_in_rendering:
                edge_label_pieces.append(edge_data['name'])
            if self._include_momentum_routing_in_rendering:
                edge_label_pieces.append(edge_data['momentum'])
            if is_LMB:
                edge_label_pieces.append('#%d'%(abs_sig[0].index(1)))
            edge_label = "|".join(edge_label_pieces)
            edge_repl_dict['edge_label'] = edge_label
            all_edge_definitions.append(edge_template%edge_repl_dict)

        repl_dict['edge_lists'] = ',\n'.join(all_edge_definitions)
        repl_dict['edge_shape_definitions'] = ',\n'.join(all_edge_shape_functions)
        return \
"""Labeled[GraphClass[{
%(edge_lists)s
},
EdgeShapeFunction -> {
%(edge_shape_definitions)s
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(edge_font_size)d, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(vertex_font_size)d, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> %(graph_layout_strategy)s,
ImageSize -> {%(width)f, %(height)f}
],"%(graph_name)s"]"""%repl_dict
    
    def draw(self, output_dir,FORM_id=None, lmb_id=None):
        """ Outputs the mathematica code for rendering this FORMGraph."""
        
        if FORM_id is not None:
            file_name = 'Graph_%04d'%FORM_id
        else:
            file_name = 'Graph_%s'%self.name
        
        if lmb_id is not None:
            file_name += '_LMB_%04d'%lmb_id

        MM_code = \
"""GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=%s;
"""%self.get_mathematica_rendering_code(FORM_id=FORM_id, lmb_id=lmb_id)
        # Export to PDF in landscape format. One graph per page for now.
        # The 1.2 multiplier accounts for margins
        MM_code += 'Export["%s.pdf", aGraph];'%(
                    file_name)
        open(pjoin(output_dir,'%s.m'%file_name),'w').write(MM_code)

    def generate_numerator_form_input(self, additional_overall_factor='', only_algebra=False):
        # create the input file for FORM
        form_diag = self.overall_factor+additional_overall_factor
        for node in self.nodes.values():
            if node['vertex_id'] < 0:
                continue

            form_diag += '*\n vx({},{},{})'.format(
                ','.join(str(p) for p in node['PDGs']),
                ','.join(node['momenta']),
                ','.join(str(i) for i in node['indices'])
            )

        for edge in self.edges.values():
            form_diag += '*\n prop({},{},{},{})'.format(
                edge['PDG'],
                edge['type'],
                edge['momentum'],
                ','.join(str(i) for i in edge['indices'])
            )

        return form_diag

    def get_topo_generator(self, specified_LMB=None):
        """ Returns a topology generator for that FORMGraph."""

        topo_edges = copy.deepcopy(self.edges)

        original_LMB = {}
        external_edges = []
        other_edges = []
        edge_name_to_key = {}
        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in topo_edges.items():
            # Fix for QGRAF pipeline
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key) 
          
            if edge_data['type'] == 'virtual':
                if not edge_data['name'].startswith('p'):
                    edge_data['name'] = 'p%s'%edge_data['name']
                other_edges.append((edge_data['name'],edge_key[0],edge_key[1]))
            else:
                if not edge_data['name'].startswith('q'):
                    edge_data['name'] = 'q%s'%edge_data['name'][1:]
                external_edges.append((edge_data['name'],edge_data['vertices'][0],edge_data['vertices'][1]))
            edge_name_to_key[edge_data['name']]=edge_key

            # Test if it is a defining edge of the lmb
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            if sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0:
                original_LMB[abs_sig[0].index(1)]=edge_data['name']

        topo_edges = external_edges+other_edges

        # Set the LMB to a sorted one
        original_LMB = sorted(list(original_LMB.items()),key=lambda e: e[0])
        assert(all(oLMBe[0]==i for i,oLMBe in enumerate(original_LMB)))
        original_LMB = [oLMBe[1] for oLMBe in original_LMB]

        topo_generator = utils.TopologyGenerator(topo_edges)

        topo_generator.generate_momentum_flow( loop_momenta = (original_LMB if specified_LMB is None else specified_LMB) )
        original_LMB = [edge_name_to_key[oLMBe] for oLMBe in original_LMB]

        return topo_generator, edge_name_to_key, original_LMB

    def derive_signatures(self):
        n_incoming = sum([1 for edge in self.edges.values() if edge['type'] == 'in'])
        n_loops = len(self.edges) - len(self.nodes) + 1

        # parse momentum
        p = re.compile(r'(^|\+|-)(k|p)(\d*)')
        for edge in self.edges.values():
            parsed = [e.groups() for e in p.finditer(edge['momentum'])]
            signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
            for pa in parsed:
                if pa[1] == 'p':
                    signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                else:
                    signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1

            edge['signature'] = signature

    def impose_signatures(self):

        for eid, e in self.edges.items():
            e['momentum'] = FORMGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
            neg_mom =  FORMGraph.momenta_decomposition_to_string([[-s for s in sp] for sp in e['signature']], set_outgoing_equal_to_incoming=False)
            for i, vi in enumerate(e['vertices']):
                e_index = self.nodes[vi]['edge_ids'].index(eid)
                mom = list(self.nodes[vi]['momenta'])
                mom[e_index] = neg_mom if i == 0 else e['momentum']
                self.nodes[vi]['momenta'] = mom

    @classmethod
    def momenta_decomposition_to_string(cls, momenta_decomposition, set_outgoing_equal_to_incoming=True, overall_sign=1, leading_sign=False):
        """ Turns ((1,0,0,-1),(1,1)) into 'k1-k4+p1+p2'"""

        res = ""
        first=True
        # The outgoing momenta are set element-wise equal to the incoming ones.
        if set_outgoing_equal_to_incoming:
            momenta_decomposition = [
                momenta_decomposition[0],
                [inp+outp for inp, outp in zip(
                    momenta_decomposition[1][:len(momenta_decomposition[1])//2],
                    momenta_decomposition[1][len(momenta_decomposition[1])//2:]
                )]
            ]
        # Fuse outgoing and incoming
        for symbol, mom_decomposition in zip(('k','p'),momenta_decomposition):
            for i_k, wgt in enumerate(mom_decomposition):
                wgt = wgt * overall_sign
                if wgt!=0:
                    if first and not leading_sign:
                        if wgt<0:
                            res+="-"
                        first=False
                    else:
                        if wgt<0:
                            res+="-"
                        else:
                            res+="+"
                    if abs(wgt)!=1:
                        res+="%d*"%abs(wgt)
                    res+="%s%d"%(symbol,(i_k+1))

        return res

    @classmethod
    def from_dict(cls, file_path):
        """ Creates a FORMGraph from a Python dict file path."""

        # TODO: Creates an instance from a Python dict dump.
        pass

    def to_dict(self, file_path=None):
        """ Outputs the FORMGraph self to a Python dict file path."""

        dict_to_dump = {
            'edges' : self.edges,
            'nodes' : self.nodes,
            'overall_factor' : self.overall_factor,
            'multiplicity' : self.multiplicity
        }
        if file_path:
            open(file_path,'w').write(pformat(dict_to_dump))
        else:
            return dict_to_dump
     

    def get_edge_scaling(self, pdg):
        # all scalings that deviate from -2
        scalings = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 11: -1, 12: -1, 13: -1}
        return scalings[abs(pdg)] if abs(pdg) in scalings else -2

    def get_node_scaling(self, pdgs):
        # only the triple gluon vertex and the ghost gluon vertex have a non-zero scaling
        if pdgs == (25, 21, 21):
            return 2
        elif pdgs == (25, 21, 21, 21) or pdgs == (21, 21, 21) or pdgs == (-82, 21, 82) or pdgs == (210, 21, 21) or pdgs == (-82, 210, 82):
            return 1
        else:
            return 0

    def generate_topology_files(self, process_name, index):
        if self.is_zero:
            return False

        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in self.edges.items():
            edge_data['name'] = 'p' + edge_data['name'] if edge_data['type'] == 'virtual' else 'q' + edge_data['name'][1:]

        # TODO: sort such that the first 4 entries are external (it seems to happen by chance now every time)
        edge_map_lin = [(e['name'], e['vertices'][0], e['vertices'][1]) for e in self.edges.values()]
        assert(e[0] != 'q' or int(e[1:]) < 5 for e in edge_map_lin)

        particle_ids = { e['name']: e['PDG'] for e in self.edges.values() }

        num_incoming = sum(1 for e in edge_map_lin if e[0][0] == 'q') // 2

        if num_incoming == 1:
            external_momenta = {'q1': [500., 0., 0., 0.], 'q2': [500., 0., 0., 0.]}
            #external_momenta = {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]}
            p = np.array(external_momenta['q1'])
        else:
            external_momenta = { 'q1': [500., 0., 0., 500.], 'q2': [500., 0., 0., -500.], 'q3': [500., 0., 0., 500.], 'q4': [500., 0., 0., -500.] }
            #external_momenta = {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]}
            p = np.array(external_momenta['q1']) + np.array(external_momenta['q2'])

        loop_momenta = []
        n_loops = len(self.edges) - len(self.nodes) + 1
        for loop_var in range(n_loops):
            lm = next((ee['name'], ee['signature'][0][loop_var]) for ee in self.edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] == 1)
            loop_momenta.append(lm)

        topo = ForestGenerator(edge_map_lin,
            process_name + '_' + str(index),
            ['q1', 'q2'][:num_incoming],
            loop_momenta_names=tuple([l for l,s in loop_momenta]),
            masses={},
            edge_weights={e['name']: self.get_edge_scaling(e['PDG']) for e in self.edges.values()},
            vertex_weights={nv: self.get_node_scaling( n['PDGs']) for nv, n in self.nodes.items()},
        )

        self.topology = topo

        self.generate_form_input(process_name, topo)

        return True


    def generate_form_input(self, process_name, topo):
        configurations = []
        uv_diagrams = []
        uv_forest = []
        topo_map = '#procedure uvmap()\n'

        def sign_prefix(s):
            return '+' if s == 1 else ('-' if s == -1 else str(s) + '*')

        def strip_plus(s):
            return str(s[1:]) if len(s) > 1 and s[0] == '+' else s

        pure_forest_counter = 0
        diag_set_uv_conf = []
        diag_momenta = []

        n_loops = len(self.edges) - len(self.nodes) + 1
        cmb_offset = 0

        conf = []
        bubble_uv_derivative = ''
        for uv_index, uv_structure in enumerate(topo.uv):
            forest_element = []

            # construct the map from the cmb/lmb to the forest basis
            forest_to_cb = []
            for lmb_index, (r, aff, extshift, _, _, _) in enumerate(zip(*uv_structure['forest_to_cb_matrix'])):
                if all(x == 0 for x in r):
                    assert(all(x == 0 for x in aff))
                    continue
                mom = ''.join('{}fmb{}'.format(sign_prefix(a), forest_index + 1) for forest_index, a in enumerate(r) if a != 0)
                # the shift should be subtracted
                shift = ''
                for cmb_index, a in enumerate(aff):
                    if a == 0:
                        continue

                    # also subtract the external momenta
                    d = self.momenta_decomposition_to_string(([0] * n_loops, cut['cuts'][cmb_index]['signature'][1]), False, a, True)
                    shift += '{}c{}{}'.format(sign_prefix(-a), cmb_index + 1, d)

                shift += self.momenta_decomposition_to_string(([0] * n_loops, extshift), True, -1, True)
                m = 'k{},{}{}'.format(lmb_index + cmb_offset + 1, strip_plus(mom), shift)
                forest_to_cb.append(m)
            if len(forest_to_cb) > 0:
                forest_element.append('cbtofmb({})'.format(','.join(forest_to_cb)))

            for uv_subgraph in uv_structure['uv_subgraphs']:
                if uv_subgraph['first_occurrence_id'] == uv_subgraph['id']:
                    if uv_subgraph['skip_pf']:
                        continue
                    rp = '*'.join('t{}'.format(i) if raised_power == 1 else 't{}^{}'.format(i, raised_power)
                            for i, (_,_,raised_power) in enumerate(uv_subgraph['loop_topo'].uv_loop_lines[0]) if raised_power != 0)
                    topo_map += '\tid uvtopo({},{},k1?,...,k{}?) = diag({},{},k1,...,k{});\n'.format(uv_subgraph['id'],
                        '1' if rp == '' else rp, uv_subgraph['graph'].n_loops, 0, uv_subgraph['id'], uv_subgraph['graph'].n_loops)

                # construct the vertex structure of the UV subgraph
                uv_loop_graph = uv_subgraph['loop_topo'] # all derived graphs have the same topo
                uv_diag_moms = ','.join(self.momenta_decomposition_to_string(lmm, False) for lmm in uv_loop_graph.loop_momentum_map)
                vertex_structure = []
                subgraph_vertices = set(v for ll in uv_loop_graph.loop_lines for v in (ll.start_node, ll.end_node))
                for v in subgraph_vertices:
                    vertex = []
                    for ll in uv_loop_graph.loop_lines:
                        for (dv, outgoing) in ((ll.start_node, 1), (ll.end_node, -1)):
                            if v != dv:
                                continue
                            loop_mom_sig = ''
                            for s, lmm in zip(ll.signature, uv_loop_graph.loop_momentum_map):
                                if s != 0:
                                    loop_mom_sig += self.momenta_decomposition_to_string(lmm, False, s * outgoing,True)
                            vertex.append(strip_plus(loop_mom_sig))
                    vertex_structure.append('vxs({})'.format(','.join(vertex)))

                uv_props = []
                for i, (ll_sig, propagators, _raised_power) in enumerate(uv_loop_graph.uv_loop_lines[0]):
                    loop_mom_sig = ''
                    loop_mom_shift = ''
                    for s, lmm, lm_shift in zip(ll_sig, uv_loop_graph.loop_momentum_map, uv_loop_graph.uv_loop_lines[1]):
                        if s != 0:
                            loop_mom_sig += self.momenta_decomposition_to_string(lmm, False, s, True)
                            loop_mom_shift += self.momenta_decomposition_to_string(lm_shift, False, s, True)

                    # the parametric shift is given in terms of external momenta of the subgraph
                    # translate the signature and param_shift to momenta of the supergraph
                    for (edge_name, param_shift, power) in propagators:
                        ext_mom_sig = ''
                        edge_mass = 'masses({})'.format(next(ee for ee in self.edges.values() if ee['name'] == edge_name)['PDG'])

                        if all(s == 0 for s in param_shift[1]):
                            for _ in range(power):
                                if loop_mom_shift == '':
                                    uv_props.append('uvprop({},t{},0,{})'.format(strip_plus(loop_mom_sig), i, edge_mass))
                                else:
                                    uv_props.append('uvprop({},t{},{},{})'.format(strip_plus(loop_mom_sig), i, strip_plus(loop_mom_shift), edge_mass))
                            continue

                        ext_mom_sig = ''
                        for (ext_index, s) in enumerate(param_shift[1]):
                            if s != 0:
                                ext_mom = uv_subgraph['graph'].edge_map_lin[uv_subgraph['graph'].ext[ext_index]][0]
                                ext_edge = next(ee for ee in self.edges.values() if ee['name'] == ext_mom)
                                sig_l, sig_ext = np.array(ext_edge['signature'][0], dtype=int), np.array(ext_edge['signature'][1], dtype=int)
                                if ext_mom_sig == '':
                                    ext_mom_sig = (s * sig_l, s * sig_ext)
                                else:
                                    ext_mom_sig = (ext_mom_sig[0] + s * sig_l, ext_mom_sig[1] + s * sig_ext)

                        ext_mom_sig = (([x for x in ext_mom_sig[0]]), ([x for x in ext_mom_sig[1]]))
                        ext_mom_sig = self.momenta_decomposition_to_string(ext_mom_sig, False, 1, True)

                        for _ in range(power):
                            uv_props.append('uvprop({},t{},{},{})'.format(loop_mom_sig, i, strip_plus(loop_mom_shift + ext_mom_sig), edge_mass))
                # it could be that there are no propagators with external momentum dependence when pinching duplicate edges
                if uv_props == []:
                    uv_props = ['1']

                if uv_diag_moms == '':
                #    # should never happen!
                #    logger.warn("No diag moms in UV graph")
                    uv_diag = 'uvtopo({})'.format(uv_subgraph['first_occurrence_id'])
                else:
                    uv_diag = 'uvtopo({},{})'.format(uv_subgraph['first_occurrence_id'], uv_diag_moms)


                uv_diag += '*{}'.format('*'.join(vertex_structure))

                uv_conf_diag = '-tmax^{}*{}*{}'.format(uv_subgraph['taylor_order'],'*'.join(uv_props),uv_diag)
                if uv_conf_diag not in uv_diagrams:
                    uv_diagrams.append(uv_conf_diag)

                uv_conf = 'uvdiag({})'.format(uv_diagrams.index(uv_conf_diag))

                sg_call = 'subgraph({}{},{},{})'.format(uv_subgraph['graph_index'],
                    (',' if len(uv_subgraph['subgraph_indices']) > 0 else '') + ','.join(str(si) for si in uv_subgraph['subgraph_indices']),
                    uv_conf, uv_diag_moms)

                forest_element.append(sg_call)

            conf.append('*'.join(forest_element))

        #cmb_offset += len(diag_info['uv'][0]['forest_to_cb_matrix'][0][0]) # add the amplitude loop count to the cmb start
        uv_forest.append('+\n\t'.join(conf))

        if len(configurations) > 0:
            configurations[-1] += '\n'

        topo_map += '#endprocedure\n'

        os.makedirs(process_name, exist_ok = True)
        with open(pjoin(process_name, '%s_in.h'%topo.name), 'w') as f:
            f.write('CTable uvdiag(0:{});\n'.format(len(uv_diagrams)))
            f.write('{}\n\n'.format('\n'.join('Fill uvdiag({}) = {};'.format(i, uv) for i,uv in enumerate(uv_diagrams))))
            f.write('CTable forest(0:{});\n'.format(len(uv_forest)))
            f.write('{}\n\n'.format('\n'.join('Fill forest({}) = {};'.format(i, uv) for i,uv in enumerate(uv_forest))))
            f.write('L F = {}\n;'.format(self.generate_numerator_form_input('*forestid(0,{})'.format(','.join('k' + str(i) for i in range(1,topo.topo.n_loops + 1))))))

    @classmethod
    def from_dict(cls, dict_file_path, first=None, merge_isomorphic_graphs=False, verbose=False, model = None, workspace=None, cuts=None):
        """ Creates a FORMGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))

        # compile the file first before importing it
        # with the same optimization flag as MG5
        # this avoids that memory isn't freed after compiling
        # when using the __import__ directly
        logger.info("Compiling imported graphs.")
        subprocess.run([sys.executable, '-O', '-m', p.stem], cwd=p.parent)
        m = __import__(p.stem)
        # Reload to avoid border effects if this is the second process generated in this python session.
        importlib.reload(m)
        logger.info("Imported {} graphs.".format(len(m.graphs)))
        sys.path.pop(0)


        # Now convert the vertex names to be integers according to QGRAF format:
        for i, g in enumerate(m.graphs):
            if isinstance(list(g['nodes'].keys())[0], str):
                new_nodes ={}
                node_names={}
                for i_node, (n_key, n) in enumerate(g['nodes'].items()):
                    node_names[n_key]=i_node+1
                    new_nodes[i_node+1] = n
                for n in g['nodes'].values():
                    n['edge_ids'] = tuple([ (node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2]) for edge_key in n['edge_ids'] ])
                new_edges = {}
                for edge_key, edge in g['edges'].items():
                    edge['vertices']=(node_names[edge_key[0]],node_names[edge_key[1]])
                    new_edges[(node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2])]=edge
                g['nodes'] = new_nodes
                g['edges'] = new_edges
            else:
                # if the diagram comes from QGRAF then we want to rectify the flow of the loop momenta such 
                # that there is consistency within the same loop line
                topo_generator = utils.TopologyGenerator([(e['name'], e['vertices'][0], e['vertices'][1]) for e in g['edges'].values()])
                topo_generator.generate_momentum_flow()
                smap = topo_generator.get_signature_map()
                for e_key, e in g['edges'].items():
                    #print("\t",e['momentum'], " -> ", FORMGraph.momenta_decomposition_to_string(smap[e['name']], set_outgoing_equal_to_incoming=True))
                    e['momentum'] = FORMGraph.momenta_decomposition_to_string(smap[e['name']], set_outgoing_equal_to_incoming=True)
                    #print(m.graphs[i]['edges'][e_key])
                for n_key, n in g['nodes'].items():
                    if n['vertex_id'] < 0:
                        continue
                    new_moms = []
                    for idx, eid in zip(n['indices'], n['edge_ids']):
                        if g['edges'][eid]['type'] == 'in':
                            sgn = 1
                        elif g['edges'][eid]['type'] == 'out':
                            sgn = -1
                        else:
                            sgn = 2*g['edges'][eid]['indices'].index(idx)-1
                    
                        if not sgn in [1, -1]:
                            raise ValueError
                        
                        signature = [sgn * x for x in smap[g['edges'][eid]['name']][0]]
                        shifts = [sgn * x for x in smap[g['edges'][eid]['name']][1]]
                        new_moms += [FORMGraph.momenta_decomposition_to_string([signature, shifts], set_outgoing_equal_to_incoming=True)]
                    #print("\t\t",n['momenta'])
                    #print("\t\t",tuple(new_moms))
                    n['momenta'] = tuple(new_moms)
                    #print("\t\t",m.graphs[i]['nodes'][n_key]['momenta'])

        full_graph_list = []
        for i, g in enumerate(m.graphs):
            if hasattr(m,'graph_names'):
                graph_name=m.graph_names[i]
            else:
                graph_name=p.stem + '_' + str(i)
            # convert to FORM supergraph
            form_graph = FORMGraph(name=graph_name, edges = g['edges'], nodes=g['nodes'], 
                        overall_factor=g['overall_factor'], multiplicity = g.get('multiplicity',1) )
            form_graph.derive_signatures()
            full_graph_list.append(form_graph)


        # TODO: do isomorphisms
        return full_graph_list


if __name__ == "__main__":
    process_name = sys.argv[1]
    graphs = FORMGraph.from_dict(process_name + '.py')

    for i, g in enumerate(graphs):
        g.generate_topology_files(process_name, i)
        g.draw(process_name)

    # copy all FORM files to the directory
    for filename in os.listdir('rqft'):
        shutil.copy(pjoin('rqft', filename), pjoin(process_name, filename))

    # generate run script        
    with open(pjoin(process_name, 'run_process.py'), 'w') as f:
        f.write("""
import sys
import os.path
import subprocess
from multiprocessing import Pool

def run_form(i):
    subprocess.run(["form", "-D", "DIAG={{}}_{{}}".format("{}", i), "-D", "MAXPOLE=3", "rqft.frm"], stdout=subprocess.DEVNULL)
    print('Done with graph {{}}'.format(i))

if __name__ == '__main__':
    with Pool() as pool:
        pool.map(run_form, range({}))

    """.format(process_name, len(graphs)))

    # generate draw script
    with open(pjoin(process_name, 'Makefile'), 'w') as f:
        f.write("""
MATH=$(shell readlink `which math`)

GRAPHS_SRC=$(wildcard Graph_*.m)
GRAPHS_PDF=$(patsubst %.m,%.pdf,$(wildcard Graph_*.m))

ALLGRAPHS=all_graphs.pdf

all: $(ALLGRAPHS)

Graph_%.pdf: Graph_%.m
	@echo "Generating graph $<..."
	$(MATH) -run < $< > /dev/null 

$(ALLGRAPHS): $(sort $(GRAPHS_PDF))
	pdfunite $^ $@

clean:
	rm -f *.pdf
    """)