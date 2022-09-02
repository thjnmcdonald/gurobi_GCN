# coding: utf-8
# !/usr/bin/env python

from os import chdir, getcwd
import torch
import pandas as pd
import numpy as np
import math
import gurobipy as gp
import argparse
import csv

from gurobipy import GRB
from smiles_to_molecular_graphs.single_molecule_conversion import process
from torch_geometric.utils import to_dense_adj
from GCN.optimization.nn_milp_builder import *
from optimization.nn_milp_builder_tom_version_for_GCN import *
from optimization.nn_milp_builder_tom_version import *
from optimization.bound_tightening.bt_lp import bt_lrr as bt_lrr_ann
from GCN.optimization.bound_tightening.bt_lp import bt_lrr as bt_lrr_gcn

parser = argparse.ArgumentParser()

parser.add_argument('--location', default='2022-08-31_12:00:14-GCN_16_1')
parser.add_argument('--time_lim', default=1000)
parser.add_argument('--mol_len', default=4)

args = parser.parse_args()

path = str(args.location)
time_lim = float(args.time_lim)
find_mol_of_length = int(args.mol_len)

if int(path[-1]) == 1:
    model_type_multi_layer_bool = False
else:
    model_type_multi_layer_bool = True

model_name = path[-8:]
gnn = model_name[0:3]
neurons = model_name[4:6]
layers = model_name[-1]

print(f'{model_name=}')
print(f'{gnn=}')
print(f'{neurons=}')
print(f'{layers=}')

n = find_mol_of_length
F = 14
d_max = 4

path_fetch = f'trained_models/{path}/{model_name}'


state_dict = torch.load(path_fetch + ".tar")
print(f'{state_dict.keys()=}')

def hard_coded(mol):
    hardcode, mol_dict = process(mol, mol_dict = True)
    edge_index = hardcode.edge_index
    adjacency_matrix = to_dense_adj(edge_index)
    hardcode_A = adjacency_matrix
    hardcode_features = hardcode.x
    rel_data = [edge_index, hardcode_A, hardcode_features]
    return rel_data, mol_dict

rel_data, mol_dict = hard_coded('CCCCl')   

def normalization_term(d_max):
    g = np.zeros((d_max + 1)**2)
    for i in range(d_max+1):
        for j in range(d_max+1):
            if (i == 0) or (j == 0):
                g[i*(d_max + 1) + j] = 0
            else:
                g[i*(d_max + 1) + j] = 1/((math.sqrt(i))*math.sqrt(j))
    return g

g = normalization_term(d_max)


def make_input_constraints(m: gp.Model, rel_data, n, F, feature_map, hard_coded = False, hard_coded_self = False):
    A = m.addVars(n, n, vtype=GRB.BINARY, name="A")
    x = m.addVars(n, F, vtype = GRB.BINARY, name = "x")
    m.update()

    if hard_coded:
        if hard_coded_self:

            edge_matrix = np.array( [[0., 1., 1., 1.],
                                    [1., 0., 1., 1.],
                                    [1., 1., 0., 1.],
                                    [1., 1., 1., 0.]])

            feature_vectors = np.array( [[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
                                        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
                                        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
                                        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]])

        else:
            features = rel_data[2]
            edge_matrix = rel_data[1][0]
            feature_vectors = rel_data[-1]
            
        for i in range(n):
            for j in range(n):
                if i == j:
                    m.addConstr(A[i, j] == 1)
                else:
                    m.addConstr((A[i, j] == int(edge_matrix[i][j])))

        m.update()

        for i in range(n):
            for j in range(F):
                m.addConstr(x[i, j] == feature_vectors[i, j])

        m.update()

    else:
        num_atoms = feature_map.count('atom_type')
        num_properties = feature_map.count('properties')
        num_hybridization = feature_map.count('hybridization')
        num_neighbours = feature_map.count('neighbours')
        num_hydrogen = feature_map.count('hydrogen')
        
        features = [num_atoms, num_properties, num_hybridization, num_neighbours, num_hydrogen]        
        features_simple = [num_atoms, num_neighbours, num_hydrogen]    
        feature_cumsum = np.cumsum(features_simple)

        m.addConstr(1 == A[0,0], name = 'min constr')

        m.addConstr(1 == A[1,1], name = 'min constr')

        m.addConstr(1 == A[1,0], name = 'min constr')
        
        m.addConstrs(((A[i,i] >= A[i + 1, i + 1])
            for i in range(n-1)), name = '8a')

        m.addConstrs((A[i,i] == gp.quicksum(x[i, s] for s in range(feature_cumsum[0])) for i in range(n)), 
        name = 'one_atom')

        # m.addConstrs((A[i,i] == gp.quicksum(x[i, s] for s in range(feature_cumsum[1], feature_cumsum[2])) for i in range(n)), 
        # name = 'one_hybrid')

        m.addConstrs((A[i,i] == gp.quicksum(x[i, s] for s in range(feature_cumsum[0], feature_cumsum[1])) for i in range(n)), 
        name = 'one_nbor')

        m.addConstrs((A[i,i] == gp.quicksum(x[i, s] for s in range(feature_cumsum[1], feature_cumsum[2])) for i in range(n)), 
        name = 'one_hydro')

        # change this later to automatically update atom types and their covalence. 
        # m.addConstrs(( (4*x[i,0] + 2*x[i,1] + 1*x[i,2] + 1*x[i,3] ==
        #     gp.quicksum( (s-11)*x[i, s] for s in range(feature_cumsum[2], feature_cumsum[3]))
        #     + gp.quicksum((t - 16) * x[i, t] for t in range(feature_cumsum[3], feature_cumsum[4]) ))
        #       for i in range(n)), name = '8f')

        m.addConstrs(( (4*x[i,0] + 2*x[i,1] + 1*x[i,2] + 1*x[i,3] ==
            gp.quicksum( (s-4)*x[i, s] for s in range(feature_cumsum[0], feature_cumsum[1]))
            + gp.quicksum((t - 9) * x[i, t] for t in range(feature_cumsum[1], feature_cumsum[2]) ))
              for i in range(n)), name = '8f')

        M_4 = n + 1
        
        # big-M value, equal to amount of different feature types in feature vectors
        M_5 = len(set(feature_map))

        list_indices = list(range(n))
        sum_numbers = [[x for x in list_indices if x != i] for i in range(n)]
        
        m.addConstrs(((M_4*A[i,i] >= gp.quicksum(A[i,j] for j in sum_numbers[i]))
             for i in range(n)), name = '8g')

        # forces neighbour feature to be equal to the actual outdegree of A
        m.addConstrs(
            ( gp.quicksum(A[i,s] for s in sum_numbers[i]) == gp.quicksum((t-4)*x[i,t] for t in range(feature_cumsum[0], feature_cumsum[1])) 
                for i in range(n)), name = 'out_degree_feature'
        )

        m.addConstrs(((A[i,i] <= gp.quicksum(A[i,j] for j in sum_numbers[i]))
             for i in range(n)), name = '8h')

        m.addConstrs((A[i,j] == A[j,i] for i in range(n) for j in range(n)), name='symmetry_breaker')
        
        m.addConstrs((A[i,i]*M_5 >= gp.quicksum(x[i,s] for s in range(len(feature_map))) for i in range(n)), name = 'close_feature')

        # we need to have connected graph
        # since A[0,0] and A[1,1] are one, we only need this for i >= 1
        sum_numbers_2 = [[x for x in list_indices if x < i] for i in range(n)]
        m.addConstrs((A[i,i] <= gp.quicksum(A[i,j] for j in sum_numbers_2[i]) for i in range(2, n)), name = 'connected_graph')

    return m, A, x


def make_norm_layers(m: gp.Model, d_max, F, n, A, x, hard_coded = False):
    M_1 = 1
    M_2 = 1
    g = normalization_term(d_max)

    c = m.addVars(n, n, (d_max + 1)**2, vtype=GRB.BINARY, name="c")
    b = m.addVars(n, n, F, vtype=GRB.BINARY, name="b")
    sp = m.addVars(n, n, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="sp")
    d_plus = m.addVars(n, lb=0, ub=(n+1), vtype=GRB.INTEGER, name="d_plus")
    p = m.addVars(n, n, ub=(d_max + 1)**2, vtype=GRB.INTEGER, name="p")
    t = m.addVars(n, n, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="t")
    m.update()

    m.addConstrs(((d_plus[i] == gp.quicksum(A[i, j] for j in range(n)))
                 for i in range(n)), name='6a')

    m.addConstrs(((p[i, k] == d_plus[i]*(d_max + 1) + d_plus[k])
                 for i in range(n) for k in range(n)), name='6b')

    m.addConstrs(((0 == p[i, k] - (gp.quicksum(t*c[i, k, t] for t in range((d_max + 1)**2))))
                 for i in range(n) for k in range(n)), name='6c')

    m.addConstrs(((1 == (gp.quicksum(c[i, k, t] for t in range((d_max + 1)**2))))
                 for i in range(n) for k in range(n)), name='6d')

    m.addConstrs(((sp[i, k] == (gp.quicksum(g[t]*c[i, k, t] for t in range((d_max + 1)**2))))
                 for i in range(n) for k in range(n)), name='6e')

    # for bi-linear case (17.g and 17.h)
    m.addConstrs(((sp[i, k] - M_1*(1 - A[i, k]) <= t[i, k])
                 for i in range(n) for k in range(n)), name='17g-a')
    m.addConstrs(((t[i, k] <= M_1*(1 - A[i, k]) + sp[i, k])
                 for i in range(n) for k in range(n)), name='17g-b')

    m.addConstrs(((-M_1*(A[i, k]) <= t[i, k])
                 for i in range(n) for k in range(n)), name='17h-a')
    m.addConstrs(((t[i, k] <= M_1*(A[i, k]))
                 for i in range(n) for k in range(n)), name='17h-b')

    m.update()

    # add hardcoded conconstraints

    atom_covalence = [4, 2, 1, 1] 

    if hard_coded:
        pass
    else: 
        # we take d_plus [i] minus A[i,i] becuase when dplus = 0 then A[ii] is also 0
        # when dplus is at least 2 then we need to substract 1. With this we save making a new var. 
        pass
        m.addConstrs((d_plus[i] - A[i,i] == gp.quicksum((j - 4) * x[i, j] for j in range(4,9)) for i in range(n)), 
            name = 'neighbour_connection')


    m.update()
    return m, t


def make_GCN_milp(m: gp.Model, bt_procedures, x, n, F, d_max,
                  rel_data=rel_data, state_dict=state_dict,
                  bilinear=False, norm_term=None, multi_layer=False, constrained_fingerprint = False):

    GCN_output_len = 32

    h_prime = m.addVars(n, GCN_output_len, lb=0,
                        vtype=GRB.CONTINUOUS, name="h_prime")
    
    h_prime_vars = [[h_prime[i, j]
                     for j in range(GCN_output_len)] for i in range(n)]
    m.update()

    # GCN layer builder
    m, bounds = build_GCN_milp_and_run_bt(  model = m, input_vars = x, output_vars = h_prime_vars, 
                                            state_dict = state_dict, n = n, d_max = d_max, 
                                            c_id='', bt_procedures = bt_procedures, bilinear = bilinear, 
                                            norm_term = norm_term, multilayer = multi_layer)
        
    m.update()

    h = m.addVars(GCN_output_len, vtype=GRB.CONTINUOUS, name='h')
    m.update()
    m.addConstrs(((((gp.quicksum(h_prime[i, j] for i in range(n))) == h[j])) for j in range(
        GCN_output_len)), name='max_pool')

    j = 0
    for (lb, ub) in (bounds[-1]):
        h[j].setAttr(gpy.GRB.Attr.LB, (d_max + 1)*(lb))
        h[j].setAttr(gpy.GRB.Attr.UB, (d_max + 1)*(ub))
        j += 1
    
    if constrained_fingerprint:
        m.addConstr(gp.quicksum(h[i] for i in range(32)) <= 20)

    m.update()

    return h_prime, h, m, bounds


def make_ANN_milp(m: gp.Model, bt_procedures, n, h, rel_data=rel_data, state_dict=state_dict):
    ann_bt_procedures = [bt_lrr_ann]

    y = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name='y')

    m.update()
    hi = [h[i] for i in range(len(h))]

    yi = [y]

    m = build_milp_and_run_bt(m, hi, yi, state_dict,
                              bt_procedures=ann_bt_procedures)
    m.update()

    return m, y


def make_model_optimize_and_save(bt_procedures=[bt_lrr_gcn],
                                 ann_bt_procedures=[bt_lrr_ann],
                                 n=find_mol_of_length, F=F, d_max=4,
                                 rel_data=rel_data, state_dict=state_dict, bilinear=False):
    
    
    #create model and set time limit
    m = gp.Model("GNN")
    m.Params.timeLimit = time_lim
    m.update()

    # create input constriants
    m, A, x = make_input_constraints(m, rel_data, n, F, mol_dict)
    m.update()

    # make the norm layers
    m, t_vars = make_norm_layers(m, d_max, F, n, A, x)
    m.update()

    # make the GCN layers
    h_prime, h, m, bounds = make_GCN_milp(  m, bt_procedures, x, n, F, d_max, 
                                            bilinear=bilinear, norm_term=t_vars, 
                                            multi_layer=model_type_multi_layer_bool)
    m.update()

    # make the ANN layers
    m, y = make_ANN_milp(m, ann_bt_procedures, n, h, rel_data, state_dict)
    m.update()

    # set objective and make bilinear
    objective = m.getVarByName('y')
    m.setObjective(objective, GRB.MAXIMIZE)
    
    if bilinear == True:
        m.Params.NonConvex = 2
    m.update()
    
    # save the linear programming formulation
    file_name = f'{model_name}_mol_len_{find_mol_of_length}'
    m.write(f'GCN_lps/{file_name}.lp')


    # optimise and get objective
    m.optimize()
    obj = m.getObjective()
    val = obj.getValue()
    
    # save solution, relevant information and temrinal output (as backup)
    m.write(f'GCN_solutions/{file_name}.sol')

    with open(f'GCN_results', mode = 'a+', newline='') as result_file:
        wr = csv.writer(result_file, quoting=csv.QUOTE_ALL)
        wr.writerow([path, gnn, layers, neurons, find_mol_of_length, m.Runtime, m.MIPGap, val])


    return val


def main():
    make_model_optimize_and_save(bilinear=True)


if __name__ == "__main__":
    main()
