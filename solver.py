# %%

import csv
import os
import sys
from docplex.mp.progress import *
from docplex.mp.progress import SolutionRecorder
import docplex.mp.model as cpx
import networkx as nx
import time
import cvxpy as cp
import numpy as np
import random
from config import get_config
from absl import flags
FLAGS = flags.FLAGS

# %%
class Solver:
    def __init__(self):
        pass


    def CPLEX_maximizing_delta_max_min(self,wk_idx,network,config,generation_number,chromosome_id,excluded_flows_set,save_flag):

        countr = 0

        list_of_flows = network.each_wk_k_u_weight[wk_idx][0].keys()
        list_of_flows = list(list_of_flows)
        list_of_flows.sort()
        #print("# flows in weight data structure",len(list_of_flows))
        list_of_flows2 = network.each_wk_each_k_each_user_pair_id_paths[wk_idx][0].keys()
        list_of_flows2 = list(list_of_flows2)
        list_of_flows2.sort()
        #print("# flows in path data structure",len(list_of_flows2))


        list_of_flows3 = network.each_wk_each_k_user_pair_ids[wk_idx][0]
        list_of_flows3 = list(list_of_flows3)
        list_of_flows3.sort()

        #print("done!")

#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="x_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] if u not in excluded_flows_set for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u] if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0}

        delta_var  = opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="delta")


        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(wk_idx,k,u,edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                # * 1

                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] if u not in excluded_flows_set
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u] if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))

        # constraint for each flow 
        # print("excluded_flows_set",excluded_flows_set)
        # print("network.each_flow_could_have_rate",network.each_flow_could_have_rate)
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0:
                    # print("we have these paths for flow u",len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]))
                    if network.fractional_max_min_flag:
                        opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])/network.each_flow_could_have_rate[u]
                        >= delta_var , ctname="each_flow_rate_{0}_{1}".format(k,u))
                        # print("we used fractional max min")
                    else:
                        opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
                        >= delta_var , ctname="each_flow_rate_{0}_{1}".format(k,u))
            

        # objective = opt_model.sum(delta_var)


        # for maximization
        opt_model.maximize(delta_var)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        # objective_value = 0
        # try:
        #     if opt_model.solution:
        #         objective_value =opt_model.solution.get_objective_value()
        # except ValueError:
        #     print(ValueError)
        objective_value =opt_model.solution.get_objective_value()
        # print("our objective dealt value is ",objective_value,network.max_edge_capacity)
        # if objective_value>0:
        sum_rates = 0
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0:
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                        sum_rates = sum_rates +x_vars[k,p].solution_value
        
        
        
        # if saving_flow_rates_flag:
            # we store results
        # if objective_value>0:
        if save_flag:
            purification_schemes_string = ""
            for pur_scheme in network.purification.allowed_purification_schemes:
                if purification_schemes_string:
                    purification_schemes_string = purification_schemes_string+","+pur_scheme
                else:
                    purification_schemes_string = pur_scheme
            list_of_flows = []
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                    if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0 :
                        w=network.each_wk_k_u_weight[wk_idx][k][u] 
                        flow_sum = 0
                        if u not in list_of_flows:
                            list_of_flows.append(u)
                        for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                            flow_sum +=x_vars[k,p].solution_value
                        with open(network.max_min_flow_rate_file_path, 'a') as newFile:                                
                            newFileWriter = csv.writer(newFile)
                            newFileWriter.writerow([generation_number,u,flow_sum,
                                                    network.fractional_max_min_flag,config.get_sum_rate_falg,
                                                    config.equal_weight_flag,chromosome_id,
                                                    network.purification.two_qubit_gate_fidelity,k,
                                                    network.running_path_selection_scheme,
                                                    
                                                    network.purification.measurement_fidelity,
                                                    network.alpha_value,x_vars[k,p].solution_value,
                                                    wk_idx,objective_value,
                                                    purification_schemes_string,network.number_of_flows])
                    # print("flow %s from %s  has sum_rate %s could_have rate %s and fraction %s  "%(u,len(list_of_flows),flow_sum,network.each_flow_could_have_rate[u],flow_sum/network.each_flow_could_have_rate[u]))
        
        return objective_value

    
    
    def cvxpy_maximizing_delta_max_min(self,wk_idx,network,config,generation_number,chromosome_id,excluded_flows_set,save_flag):
        # Define variables
        x_vars = {
            (k, p): cp.Variable(name=f"x_{k}_{p}", nonneg=True)
            for k in network.each_wk_organizations[wk_idx]
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
            if u not in excluded_flows_set
            for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
            if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]) > 0
        }
        delta_var = cp.Variable(name="delta", nonneg=True)
    
        # Objective: Maximize delta_var
        objective = cp.Maximize(delta_var)
    
        # Constraints
        constraints = []
    
        # Edge capacity constraints
        for edge in network.set_E:
            constraints.append(
                cp.sum([
                        x_vars[k, p]
                        * network.purification.get_required_purification_EPR_pairs(
                        wk_idx, k, u, edge, p,
                        network.purification.get_each_wk_k_u_threshold(wk_idx, k, u)
                    )
                    for k in network.each_wk_organizations[wk_idx]
                    for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                    if u not in excluded_flows_set
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                    if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]) > 0
                    if network.check_path_include_edge(edge, p)
                ])
                <= network.each_edge_capacity[edge]
            )
    
        # Flow constraints
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]) > 0:
                    if network.fractional_max_min_flag:
                        constraints.append(
                            cp.sum([x_vars[k, p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]]) /
                            network.each_flow_could_have_rate[u]
                            >= delta_var
                        )
                    else:
                        constraints.append(
                            cp.sum([x_vars[k, p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]])
                            >= delta_var
                        )
    
        # Solve the problem
        problem = cp.Problem(objective, constraints)
        result = problem.solve()
    
        # Extract solution
        if problem.status == cp.OPTIMAL:
            x_values = {key: var.value for key, var in x_vars.items()}
            delta_value = delta_var.value
            return delta_value
        else:
            return 0


    
    def CPLEX_max_could_have_rate(self,wk_idx,network,config,genetic_alg,generation_number,
                             run_id,chromosome_id,scheme,max_min_fairness_flag,target_flow_ids):
        
        zero_path_flag = False
#         print("**************** start **************** ")
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                user_pair = network.each_wk_k_id_pair[wk_idx][k][u]
#                 print("wk %s k %s flow %s using %s paths "%(wk_idx,k,u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])))
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                     print("using path %s "%(p))
                    path_basic_fidelity = network.purification.each_path_basic_fidelity[p]
                    F_th = network.purification.each_path_flow_target_fidelity[p]
                    path_dis_strategy_edge_Fth = network.purification.each_path_edge_level_Fth[p]
                    path_edges = network.set_of_paths[p]
                    
                    
                    edge_fidelity = network.purification.each_edge_fidelity[path_edges[0]]
                    p_lenght = network.purification.get_path_actual_length(path_edges,network.set_of_virtual_links)
                
                    new_target = round((3*(4/3*F_th-1/3)**(1/p_lenght)+1)/4,3)
                    edge_g_value = network.purification.each_path_edge_level_g[p]
                    end_g_value = network.purification.each_path_end_level_g[p]
                   
                        
                    path_fidelity_after_edge_level_dis = network.purification.compute_e2e_fidelity(path_dis_strategy_edge_Fth,path_edges,network.set_of_virtual_links)
                if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])==0:# there is no path to deliver any
                    zero_path_flag = True
#         print("**************** end **************** ")
#         time.sleep(10)
        if zero_path_flag and network.optimization_problem_with_minimum_rate:
#             print("zero path constraint")
            return 0
#         else:
#             print("not zero path constraint")
        countr = 0

        
        
#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="w_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]}

       
        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(wk_idx,k,u,edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))
                

        # constraint for each flow 
        if network.optimization_problem_with_minimum_rate:
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    flow_weight = network.each_wk_k_u_weight[wk_idx][k][u]
#                     print("we are setting constraint min flow for k %s flow %s "%(k,u))

                    #for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
    #                 print("for org %s u %s we have %s paths user pair %s "%(k,u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),network.each_id_pair[u]))
                    flow_pair = network.each_wk_k_id_pair[wk_idx][k][u]
                    each_flow_min_rate = network.min_flow_rate
                    if config.get_flows_minimum_rate_flag:
                        each_flow_min_rate = network.each_flow_minimum_rate[wk_idx,k,flow_pair]
#                     print("for flow %s pair %s in wk %s k %s we set minimum rate %s "%(u,flow_pair,wk_idx,k,each_flow_min_rate))
                    each_flow_min_rate = 10
                    opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
                        >= each_flow_min_rate
                     , ctname="min_flow_rate_{0}_{1}".format(k,u))
        if network.optimization_problem_with_maximum_rate:
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
#                     this_flow_max_rate = network.max_flow_rate
                    
                    this_flow_max_rate = network.each_wk_k_u_max_rate_constraint[wk_idx,k,u]
#                     this_flow_max_rate = 1000
#                     print("we are setting constraint max flow for k %s flow %s "%(k,u))
                    if config.updating_max_rate_flag and this_flow_max_rate<config.up_max_rate_value:
                        this_flow_max_rate = random.randint(this_flow_max_rate,config.up_max_rate_value)
                    opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] 
                                  for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
                        <= this_flow_max_rate
                     , ctname="max_flow_rate_{0}_{1}".format(k,u))

        if max_min_fairness_flag:
            objective = opt_model.sum(x_vars[k,p]*network.q_value**(network.get_path_length(p)-1)
                              for k in network.each_wk_organizations[wk_idx]
                              for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] if u in target_flow_ids 
                              for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                              )
                
        else:
            objective = opt_model.sum(x_vars[k,p]*network.each_wk_k_weight[wk_idx][k] *network.each_wk_k_u_weight[wk_idx][k][u]*network.q_value**(network.get_path_length(p)-1)
                              for k in network.each_wk_organizations[wk_idx]
                              for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] 
                              for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                              )


        # for maximization
        opt_model.maximize(objective)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        objective_value = 0
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)

        flow_raw_rate = 0
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                if u in target_flow_ids:
                    flow_weight=network.each_wk_k_u_weight[wk_idx][k][u] 
                    user_pair = network.each_wk_k_id_pair[wk_idx][k][u]
                    user_pair_str = str(user_pair[0])+":"+str(user_pair[1])
                    k_weight = network.each_wk_k_weight[wk_idx][k]
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                        flow_raw_rate+=x_vars[k,p].solution_value
                            

                                
                                
        return objective_value,flow_raw_rate
        
    def CPLEX_maximizing_EGR(self,wk_idx,network,config,genetic_alg,generation_number,
                             run_id,chromosome_id,scheme,equal_weight_flag,excluded_flows_set,save_flag):
        
        zero_path_flag = False
#         print("**************** start **************** ")
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                user_pair = network.each_wk_k_id_pair[wk_idx][k][u]
#                 print("wk %s k %s flow %s using %s paths "%(wk_idx,k,u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])))
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                     print("using path %s "%(p))
                    path_basic_fidelity = network.purification.each_path_basic_fidelity[p]
                    F_th = network.purification.each_path_flow_target_fidelity[p]
                    path_dis_strategy_edge_Fth = network.purification.each_path_edge_level_Fth[p]
                    path_edges = network.set_of_paths[p]
                    
                    
                    edge_fidelity = network.purification.each_edge_fidelity[path_edges[0]]
                    p_lenght = network.purification.get_path_actual_length(path_edges,network.set_of_virtual_links)
                
                    new_target = round((3*(4/3*F_th-1/3)**(1/p_lenght)+1)/4,3)
                    edge_g_value = network.purification.each_path_edge_level_g[p]
                    end_g_value = network.purification.each_path_end_level_g[p]
                   
                        
                    path_fidelity_after_edge_level_dis = network.purification.compute_e2e_fidelity(path_dis_strategy_edge_Fth,path_edges,network.set_of_virtual_links)

                if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])==0:# there is no path to deliver any
                    zero_path_flag = True
#         print("**************** end **************** ")
#         time.sleep(10)
        if zero_path_flag and network.optimization_problem_with_minimum_rate:
#             print("zero path constraint")
            return 0
#         else:
#             print("not zero path constraint")
        countr = 0

        
        
#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="w_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]}

       
        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(wk_idx,k,u,edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))
                

        # constraint for each flow 
        if network.optimization_problem_with_minimum_rate:
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    flow_weight = network.each_wk_k_u_weight[wk_idx][k][u]
#                     print("we are setting constraint min flow for k %s flow %s "%(k,u))

                    #for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
    #                 print("for org %s u %s we have %s paths user pair %s "%(k,u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),network.each_id_pair[u]))
                    flow_pair = network.each_wk_k_id_pair[wk_idx][k][u]
                    each_flow_min_rate = network.min_flow_rate
                    if config.get_flows_minimum_rate_flag:
                        each_flow_min_rate = network.each_flow_minimum_rate[wk_idx,k,flow_pair]
#                     print("for flow %s pair %s in wk %s k %s we set minimum rate %s "%(u,flow_pair,wk_idx,k,each_flow_min_rate))
                    each_flow_min_rate = 10
                    opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
                        >= each_flow_min_rate
                     , ctname="min_flow_rate_{0}_{1}".format(k,u))
        if network.optimization_problem_with_maximum_rate:
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
#                     this_flow_max_rate = network.max_flow_rate
                    
                    this_flow_max_rate = network.each_wk_k_u_max_rate_constraint[wk_idx,k,u]
#                     this_flow_max_rate = 1000
#                     print("we are setting constraint max flow for k %s flow %s "%(k,u))
                    if config.updating_max_rate_flag and this_flow_max_rate<config.up_max_rate_value:
                        this_flow_max_rate = random.randint(this_flow_max_rate,config.up_max_rate_value)
                    opt_model.add_constraint(
                    opt_model.sum(x_vars[k,p] 
                                  for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
                        <= this_flow_max_rate
                     , ctname="max_flow_rate_{0}_{1}".format(k,u))


        if equal_weight_flag:
            objective = opt_model.sum(x_vars[k,p] *network.q_value**(network.get_path_length(p)-1)
                          for k in network.each_wk_organizations[wk_idx]
                          for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] 
                          for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                          )
        else:
            objective = opt_model.sum(x_vars[k,p]*network.each_wk_k_weight[wk_idx][k] *network.each_wk_k_u_weight[wk_idx][k][u]*network.q_value**(network.get_path_length(p)-1)
                          for k in network.each_wk_organizations[wk_idx]
                          for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] 
                          for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                          )


        # for maximization
        opt_model.maximize(objective)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        objective_value = 0
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)


#         opt_model.clear()
#         time.sleep(20)
        if save_flag:
            purification_schemes_string = ""
            for pur_scheme in network.purification.allowed_purification_schemes:
                if purification_schemes_string:
                    purification_schemes_string = purification_schemes_string+","+pur_scheme
                else:
                    purification_schemes_string = pur_scheme
            list_of_flows = []
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                    if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])>0 :
                        w=network.each_wk_k_u_weight[wk_idx][k][u] 
                        flow_sum = 0
                        if u not in list_of_flows:
                            list_of_flows.append(u)
                        for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                            flow_sum +=x_vars[k,p].solution_value
                        with open(network.max_min_flow_rate_file_path, 'a') as newFile:                                
                            newFileWriter = csv.writer(newFile)
                            newFileWriter.writerow([generation_number,u,flow_sum,
                                                    network.fractional_max_min_flag,config.get_sum_rate_falg,
                                                    config.equal_weight_flag,chromosome_id,
                                                    network.purification.two_qubit_gate_fidelity,k,
                                                    network.running_path_selection_scheme,
                                                    
                                                    network.purification.measurement_fidelity,
                                                    network.alpha_value,x_vars[k,p].solution_value,
                                                    wk_idx,objective_value,
                                                    purification_schemes_string,network.number_of_flows])


                                
        return objective_value

 
        
    def cvxpy_maximizing_EGR(self,wk_idx,network,config,genetic_alg,generation_number,
                             run_id,chromosome_id,scheme,equal_weight_flag,excluded_flows_set,save_flag):
        # print("we are in cvxpy_maxiizingEGR function!!")
        
        # Define variables
        x_vars = {
                (k, p): cp.Variable(name=f"x_{k}_{p}", nonneg=True)
                for k in network.each_wk_organizations[wk_idx]
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u] if u not in excluded_flows_set
            }
    
        # Objective: Maximize sum of log of sum of rates for each flow
        flow_rates = {
            (k, u): cp.sum([x_vars[k, p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]])
            if len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]) > 0 else cp.Constant(0)
            for k in network.each_wk_organizations[wk_idx]
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] if u not in excluded_flows_set 
        }
        objective = cp.Maximize(cp.sum([cp.log(flow_rates[k, u]) for (k, u) in flow_rates]))
    
        # Constraints
        constraints = []
    
        # Edge constraints
        for edge in network.set_E:
            edge_constraint = cp.sum([
                x_vars[k, p]
                * network.purification.get_required_purification_EPR_pairs(
                    wk_idx, k, u, edge, p,
                    network.purification.get_each_wk_k_u_threshold(wk_idx, k, u)
                )
                for k in network.each_wk_organizations[wk_idx]
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] if u not in excluded_flows_set 
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                if network.check_path_include_edge(edge, p) 
            ])
            constraints.append(edge_constraint <= network.each_edge_capacity[edge])
    
        # Flow constraints
        # for k in network.each_wk_organizations[wk_idx]:
        #     for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
        #         if u not in excluded_flows_set:
        #             min_rate = network.each_flow_minimum_rate.get((wk_idx, k, u), 0)
        #             max_rate = network.each_wk_k_u_max_rate_constraint.get((wk_idx, k, u), float('inf'))
    
        #             # Sum of rates over paths for the flow
        #             flow_sum = cp.sum([x_vars[k, p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]])
    
        #             # Minimum rate constraint
        #             constraints.append(flow_sum >= min_rate)
    
        #             # Maximum rate constraint
        #             constraints.append(flow_sum <= max_rate)
    
        # Solve the problem
        problem = cp.Problem(objective, constraints)
        try:
            result = problem.solve()
        except:
            try:
                result = 0
                x_values = {key: var.value for key, var in x_vars.items()}
                flow_rate_values = {key: rate.value for key, rate in flow_rates.items()}
            except:
                return 0
            #print(ValueError)
        # Extract solution
        if problem.status == cp.OPTIMAL:
            x_values = {key: var.value for key, var in x_vars.items()}
            flow_rate_values = {key: rate.value for key, rate in flow_rates.items()}

            #return result, x_values, flow_rate_values
        else:
            # print("no optimal status!")
            # time.sleep(6)
            # return None, None, None
            try:
                x_values = {key: var.value for key, var in x_vars.items()}
                flow_rate_values = {key: rate.value for key, rate in flow_rates.items()}
                result = 0
            except:
                return 0
        if save_flag:
            purification_schemes_string = ""
            for pur_scheme in network.purification.allowed_purification_schemes:
                if purification_schemes_string:
                    purification_schemes_string = purification_schemes_string+","+pur_scheme
                else:
                    purification_schemes_string = pur_scheme
            list_of_flows = []
            for flow_rate_value,flow_sum in flow_rate_values.items():
                u = flow_rate_value[0]
                # flow_sum = flow_rate_value[1]
                # print("u and flow sum rate",u,flow_sum,result)
                if u not in excluded_flows_set and len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][0][u])>0 :
                    with open(network.max_min_flow_rate_file_path, 'a') as newFile:                                
                        newFileWriter = csv.writer(newFile)
                        newFileWriter.writerow([generation_number,u,flow_sum,
                                                network.fractional_max_min_flag,config.get_sum_rate_falg,
                                                config.equal_weight_flag,chromosome_id,
                                                network.purification.two_qubit_gate_fidelity,0,
                                                network.running_path_selection_scheme,
                                                network.purification.measurement_fidelity,
                                                network.alpha_value,flow_sum,
                                                wk_idx,result,
                                                purification_schemes_string,network.number_of_flows])
        return result
    def CPLEX_maximizing_EGR_minimum_rate_constraint_gradient_based_version(self,wk_idx,network,run_step_number,chromosome_id,each_wk_k_u_path_ids):
        
        
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
#                 print("for wk %s k %s u %s we have these many users %s "%(wk_idx,k,u,(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k].keys())))
                if len(each_wk_k_u_path_ids[wk_idx][k][u])==0:
                    return 0
#                 for p in each_wk_k_u_path_ids[wk_idx][k][u]:
#                     F_threshold = network.purification.each_path_flow_target_fidelity[p]
#                     path_edges = network.set_of_paths[p]
#                     p_lenght = len(path_edges)
#                     new_target = round((3*(4/3*F_threshold-1/3)**(1/p_lenght)+1)/4,3)
#                     print("for wk %s k %s u %s p %s ( from %s) to reach Fth %s (up to %s) edge-level_g=%s end-level_g=%s "%(wk_idx,k,u,p,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),F_threshold,new_target,network.purification.each_path_edge_level_g[p],network.purification.each_path_end_level_g[p]))
        
        
        
        
#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="w_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] for p in each_wk_k_u_path_ids[wk_idx][k][u]}

#         for k in network.each_wk_organizations[wk_idx]:
#             for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
#                 for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                     flag = False
#                     for edge in network.set_E:
#                         if network.check_path_include_edge(edge,p):
#                             if not flag:
#                                 print("organization %s user %s #paths %s  p %s |p|=%s F = %s Fth=%s  g(.)=%s"%(k,
#                                                                                 u,network.num_of_paths,p,
#                                                                                 network.each_path_legth[p],
#                                                                                 network.purification.each_path_basic_fidelity[p],
#                                                                                 network.purification.get_each_wk_k_u_threshold(wk_idx,k,u),
#                     network.purification.get_required_purification_EPR_pairs(edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))))
#                                 flag = True
                          

    #     time.sleep(9)
       
        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(wk_idx,k,u,edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in each_wk_k_u_path_ids[wk_idx][k][u]
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))
                
#             else:
#                 opt_model.add_constraint(
#                     opt_model.sum(x_vars[k,p]*network.each_wk_k_u_weight[wk_idx][k][u] *
#                     network.get_required_edge_level_purification_EPR_pairs(edge,p,network.purification.each_wk_k_u_fidelity_threshold[k][u],wk_idx)
#                     for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
#                     for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
#                     if network.check_path_include_edge(edge,p))
#                      <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))
#                 print("from edge %s to F %s divide capacity by this %s"%(edge,network.each_wk_k_fidelity_threshold[0],
#                           network.get_required_edge_level_purification_EPR_pairs
#                           (edge,0,network.each_wk_k_fidelity_threshold[0],0)))

        # constraint for each flow 

        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                flow_weight = network.each_wk_k_u_weight[wk_idx][k][u]
                #for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                 print("for org %s u %s we have %s paths user pair %s "%(k,u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),network.each_id_pair[u]))
                opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] for p in each_wk_k_u_path_ids[wk_idx][k][u])
                    >= network.min_flow_rate
                 , ctname="min_flow_rate_{0}_{1}".format(k,u))

        objective = opt_model.sum(x_vars[k,p]*network.each_wk_k_weight[wk_idx][k] *network.each_wk_k_u_weight[wk_idx][k][u]*network.q_value**(network.get_path_length(p)-1)
                              for k in network.each_wk_organizations[wk_idx]
                              for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] 
                              for p in each_wk_k_u_path_ids[wk_idx][k][u]
                              )


        # for maximization
        opt_model.maximize(objective)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        objective_value = 0
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)


#         opt_model.clear()
#         time.sleep(20)
        
        if objective_value>0:
            if run_step_number  in [0,200,500,900,1000,1500,2000] and network.get_flow_rates_info:
                for k in network.each_wk_organizations[wk_idx]:
                    for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                        w=network.each_wk_k_u_weight[wk_idx][k][u] 
                        for p in each_wk_k_u_path_ids[wk_idx][k][u]:
                            F_th = network.purification.each_path_flow_target_fidelity[p]
                            path_edges = network.set_of_paths[p]
                            p_lenght = len(path_edges)
                            new_target = round((3*(4/3*F_th-1/3)**(1/p_lenght)+1)/4,3)
                            edge_g_value = network.purification.each_path_edge_level_g[p]
                            end_g_value = network.purification.each_path_end_level_g[p]
                            with open(network.path_variables_file_path, 'a') as newFile:                                
                                newFileWriter = csv.writer(newFile)
                                newFileWriter.writerow([network.running_path_selection_scheme,
                                                            run_step_number,chromosome_id,k,u,p,
                                                            network.purification.two_qubit_gate_fidelity,
                                                            network.purification.measurement_fidelity,
                                                            network.alpha_value,x_vars[k,p].solution_value,
                                                            wk_idx,objective_value,F_th,new_target,
                                                            edge_g_value,end_g_value])
    
        return objective_value
    
    
    
    def CPLEX_maximizing_EGR_old_version(self,wk_idx,network,run_step_number,chromosome_id):
#         for path_id,path_edges in network.set_of_paths.items():
#             new_edges = []
#             for edge in path_edges:
#                 if edge not in network.set_of_virtual_links:
#                     new_edges.append(edge)
#             print("path id %s length %s basic Fidelity %s "%(path_id,len(new_edges),network.purification.each_path_basic_fidelity[path_id]))
            #time.sleep(2)
        #print("scheme ",network.)
#         print("network.each_wk_k_weight",wk_idx,network.each_wk_k_weight)
#         countr = 0
#         for path_id,path_edges in network.set_of_paths.items():
#             capacities = []
#             virtual_link_counter = 0
#             for edge in path_edges:
#                 if edge in network.set_of_virtual_links or (edge[1],edge[0]) in network.set_of_virtual_links:
#                     virtual_link_counter+=1
#                 capacities.append(network.each_edge_capacity[edge])
#             print("path id %s length %s %s with %s virtual F %s capac. %s "%(path_id,len(path_edges),len(path_edges)-virtual_link_counter,virtual_link_counter,network.purification.each_path_basic_fidelity[path_id],min(capacities)))
#         time.sleep(15)
        countr = 0
#         for k in network.each_wk_organizations[wk_idx]:
#                 for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
# #                     if u <40:
# # #                         print("for flow %s we have %s paths "%(u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])))
# #                         #path_ids = network.each_user_pair_all_paths[u]
# # #                         print("selected paths:     ",network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
# #                         for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
# #                             print("flow %s path legnth %s"%(u,len(network.set_of_paths[p])))
# #                         print("these are all paths ",len(path_ids),path_ids)
# #                         print("flow F %s weight %s "%(network.purification.each_wk_k_u_fidelity_threshold[wk_idx][k][u],network.each_wk_k_u_weight[wk_idx][k][u]))
#                     for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                         if countr<150:
#                             Fth = network.purification.get_each_wk_k_u_threshold(wk_idx,k,u)
#                             path_f = network.purification.each_path_basic_fidelity[p]
#                             for edge in network.set_of_paths[p]:
#                                 g_value = network.purification.get_required_purification_EPR_pairs(edge,p,Fth)
#                             print("for flow %s we use path %s with F %s  to reach Fth %s we need g(.) %s"%(u,p,path_f,Fth,g_value))
#                     countr+=1
        list_of_flows = network.each_wk_k_u_weight[wk_idx][0].keys()
        list_of_flows = list(list_of_flows)
        list_of_flows.sort()
        #print("# flows in weight data structure",len(list_of_flows))
        list_of_flows2 = network.each_wk_each_k_each_user_pair_id_paths[wk_idx][0].keys()
        list_of_flows2 = list(list_of_flows2)
        list_of_flows2.sort()
        #print("# flows in path data structure",len(list_of_flows2))


        list_of_flows3 = network.each_wk_each_k_user_pair_ids[wk_idx][0]
        list_of_flows3 = list(list_of_flows3)
        list_of_flows3.sort()
        #print("# flows in each k flows",len(list_of_flows3))
#         if len(list_of_flows2)!=len(list_of_flows3) or list_of_flows3!=list_of_flows2:
#             print("from each k flows" ,list_of_flows3)
#             print("from each k flow paths ", list_of_flows2)


        try:
            if len(list_of_flows2)!=len(list_of_flows3) or list_of_flows3!=list_of_flows2:
                print("ERRRRRRRRRRRRRRRRROOOOOOOOORRRR!!!!!!!!!!")
                print("from each k flows" ,len(list_of_flows3),list_of_flows3)
                print("from each k flow paths ",len(list_of_flows2), list_of_flows2)
                import pdb
                pdb.set_trace()
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:                
                    #print("we are k %s u %s"%(k,u))
                    #print("network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]",network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
#                     print("u %s paths %s # paths %s allowed %s"%
#                     (u,network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u],len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),network.num_of_paths))
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                        w =network.each_wk_k_weight[wk_idx][k]
                        fw = network.each_wk_k_u_weight[wk_idx][k][u]
#                         print("wk %s k %s w %s user %s w %s path %s"%
#                         (wk_idx,k,network.each_wk_k_weight[wk_idx][k],u,network.each_wk_k_u_weight[wk_idx][k][u],p))
    #                     print("edges of the path",network.set_of_paths[p])
        except:
            print("********************************************************")
            print("********************************************************")
            print("********************************************************")
            print("********************************************************")
            import pdb
            pdb.set_trace()

        #print("done!")

#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="w_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]}


    #     for k in network.K:
    #         for u in network.each_k_user_pairs[k]:
    #             for p in network.each_k_u_paths[k][u]:
    #                 print("organization %s user %s #paths %s cost %s p %s path %s"%(k,u,network.num_of_paths,network.each_link_cost_metric,p,network.set_of_paths[p]))

    #     time.sleep(9)

        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))

#             else:
#                 opt_model.add_constraint(
#                     opt_model.sum(x_vars[k,p]*network.each_wk_k_u_weight[wk_idx][k][u] *
#                     network.get_required_edge_level_purification_EPR_pairs(edge,p,network.purification.each_wk_k_u_fidelity_threshold[k][u],wk_idx)
#                     for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
#                     for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
#                     if network.check_path_include_edge(edge,p))
#                      <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))
#                 print("from edge %s to F %s divide capacity by this %s"%(edge,network.each_wk_k_fidelity_threshold[0],
#                           network.get_required_edge_level_purification_EPR_pairs
#                           (edge,0,network.each_wk_k_fidelity_threshold[0],0)))

        objective = opt_model.sum(x_vars[k,p]*network.each_wk_k_weight[wk_idx][k]*network.each_wk_k_u_weight[wk_idx][k][u] *network.q_value**(network.get_path_length(p)-1)
                              for k in network.each_wk_organizations[wk_idx]
                              for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] 
                              for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                              )


        # for maximization
        opt_model.maximize(objective)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        objective_value = 0
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)


#         opt_model.clear()
#         time.sleep(20)

       
        if run_step_number  in [0,50,100,500,900,1000] and network.get_flow_rates_info:
            # we store results
            if objective_value>0:
                purification_schemes_string = ""
                for pur_scheme in network.purification.allowed_purification_schemes:
                    if purification_schemes_string:
                        purification_schemes_string = purification_schemes_string+","+pur_scheme
                    else:
                        purification_schemes_string = pur_scheme

                for k in network.each_wk_organizations[wk_idx]:
                    for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                        w=network.each_wk_k_u_weight[wk_idx][k][u] 
                        for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                             if x_vars[k,p].solution_value>1:
#                                 edge_capacities = []
#                                 for edge in network.set_of_paths[p]:
#                                     if network.check_path_include_edge(edge,p):
#                                         edge_capacity = network.each_edge_capacity[edge]
#                                         Fth = network.purification.get_each_wk_k_u_threshold(wk_idx,k,u)
#                                         path_f = network.purification.each_path_basic_fidelity[p]
#                                         edge_capacities.append(edge_capacity)
#                                         g_value = network.purification.get_required_purification_EPR_pairs(edge,p,Fth)
#                                 print("sum %s solver %s for x_vars[%s,%s] we have %s  which *w*g(.) %s = %s < %s"%(round(sum_rates,2),round(objective_value,2),k,p,x_vars[k,p].solution_value,g_value,g_value*w*x_vars[k,p].solution_value,round(min(edge_capacities),)))
                            with open(network.path_variables_file_path, 'a') as newFile:                                
                                newFileWriter = csv.writer(newFile)
                                newFileWriter.writerow([network.running_path_selection_scheme,
                                                        run_step_number,chromosome_id,k,u,p,
                                                        network.purification.two_qubit_gate_fidelity,
                                                        network.purification.measurement_fidelity,
                                                        network.alpha_value,x_vars[k,p].solution_value,
                                                        wk_idx,objective_value,
                                                        purification_schemes_string,network.number_of_flows])

#         for edge in network.set_E:
#             #if network.end_level_purification_flag:
#             sum_value = 0
#             for k in network.each_wk_organizations[wk_idx]:
#                 for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
#                     for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                         if network.check_path_include_edge(edge,p):
#                             sum_value = sum_value + x_vars[k,p].solution_value*network.each_wk_k_u_weight[wk_idx][k][u] * network.purification.get_required_purification_EPR_pairs(edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
#             if sum_value>1:
#                 print("this rate %s on edge %s is less than %s "%(sum_value,edge,network.each_edge_capacity[edge]))

        return objective_value

    

    def CPLEX_maximizing_delta(self,wk_idx,network,run_step_number,chromosome_id):
#         for path_id,path_edges in network.set_of_paths.items():
#             new_edges = []
#             for edge in path_edges:
#                 if edge not in network.set_of_virtual_links:
#                     new_edges.append(edge)
#             print("path id %s length %s basic Fidelity %s "%(path_id,len(new_edges),network.purification.each_path_basic_fidelity[path_id]))
            #time.sleep(2)
        #print("scheme ",network.)
#         print("network.each_wk_k_weight",wk_idx,network.each_wk_k_weight)
#         countr = 0
#         for path_id,path_edges in network.set_of_paths.items():
#             capacities = []
#             virtual_link_counter = 0
#             for edge in path_edges:
#                 if edge in network.set_of_virtual_links or (edge[1],edge[0]) in network.set_of_virtual_links:
#                     virtual_link_counter+=1
#                 capacities.append(network.each_edge_capacity[edge])
#             print("path id %s length %s %s with %s virtual F %s capac. %s "%(path_id,len(path_edges),len(path_edges)-virtual_link_counter,virtual_link_counter,network.purification.each_path_basic_fidelity[path_id],min(capacities)))
#         time.sleep(15)
        countr = 0
#         for k in network.each_wk_organizations[wk_idx]:
#                 for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
# #                     if u <40:
# # #                         print("for flow %s we have %s paths "%(u,len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])))
# #                         #path_ids = network.each_user_pair_all_paths[u]
# # #                         print("selected paths:     ",network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
# #                         for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
# #                             print("flow %s path legnth %s"%(u,len(network.set_of_paths[p])))
# #                         print("these are all paths ",len(path_ids),path_ids)
# #                         print("flow F %s weight %s "%(network.purification.each_wk_k_u_fidelity_threshold[wk_idx][k][u],network.each_wk_k_u_weight[wk_idx][k][u]))
#                     for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
#                         if countr<150:
#                             Fth = network.purification.get_each_wk_k_u_threshold(wk_idx,k,u)
#                             path_f = network.purification.each_path_basic_fidelity[p]
#                             for edge in network.set_of_paths[p]:
#                                 g_value = network.purification.get_required_purification_EPR_pairs(edge,p,Fth)
#                             print("for flow %s we use path %s with F %s  to reach Fth %s we need g(.) %s"%(u,p,path_f,Fth,g_value))
#                     countr+=1
        list_of_flows = network.each_wk_k_u_weight[wk_idx][0].keys()
        list_of_flows = list(list_of_flows)
        list_of_flows.sort()
        #print("# flows in weight data structure",len(list_of_flows))
        list_of_flows2 = network.each_wk_each_k_each_user_pair_id_paths[wk_idx][0].keys()
        list_of_flows2 = list(list_of_flows2)
        list_of_flows2.sort()
        #print("# flows in path data structure",len(list_of_flows2))


        list_of_flows3 = network.each_wk_each_k_user_pair_ids[wk_idx][0]
        list_of_flows3 = list(list_of_flows3)
        list_of_flows3.sort()


        try:
            if len(list_of_flows2)!=len(list_of_flows3) or list_of_flows3!=list_of_flows2:
                print("ERRRRRRRRRRRRRRRRROOOOOOOOORRRR!!!!!!!!!!")
                print("from each k flows" ,len(list_of_flows3),list_of_flows3)
                print("from each k flow paths ",len(list_of_flows2), list_of_flows2)
                import pdb
                pdb.set_trace()
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:                
                    #print("we are k %s u %s"%(k,u))
                    #print("network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]",network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])
#                     print("u %s paths %s # paths %s allowed %s"%
#                     (u,network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u],len(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]),network.num_of_paths))
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                        w =network.each_wk_k_weight[wk_idx][k]
                        fw = network.each_wk_k_u_weight[wk_idx][k][u]
#                         print("wk %s k %s w %s user %s w %s path %s"%
#                         (wk_idx,k,network.each_wk_k_weight[wk_idx][k],u,network.each_wk_k_u_weight[wk_idx][k][u],p))
    #                     print("edges of the path",network.set_of_paths[p])
        except:
            print("********************************************************")
            print("********************************************************")
            print("********************************************************")
            print("********************************************************")
            import pdb
            pdb.set_trace()

        #print("done!")

#         print("network.max_edge_capacity",network.max_edge_capacity,type(network.max_edge_capacity))
        opt_model = cpx.Model(name="inter_organization_EGR")
        x_vars  = {(k,p): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                                  name="w_{0}_{1}".format(k,p))  for k in network.each_wk_organizations[wk_idx]
                   for u in network.each_wk_each_k_user_pair_ids[wk_idx][k] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]}

        delta_var  = opt_model.continuous_var(lb=0, ub= network.num_of_paths*network.max_edge_capacity,
                                  name="delta")


        #Edge constraint
        for edge in network.set_E:
            #if network.end_level_purification_flag:
            opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] 
                * network.purification.get_required_purification_EPR_pairs(edge,p,network.purification.get_each_wk_k_u_threshold(wk_idx,k,u))
                for k in network.each_wk_organizations[wk_idx] for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]
                for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]
                if network.check_path_include_edge(edge,p))
                 <= network.each_edge_capacity[edge], ctname="edge_capacity_{0}".format(edge))

        # constraint for each flow 
        
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                flow_weight = network.each_wk_k_u_weight[wk_idx][k][u]
                opt_model.add_constraint(
                opt_model.sum(x_vars[k,p] for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u])*flow_weight
                    >= delta_var
                 , ctname="each_flow_rate_{0}_{1}".format(k,u))
            

        objective = opt_model.sum(delta_var)


        # for maximization
        opt_model.maximize(objective)

    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()


#         print('docplex.mp.solution',opt_model.solution)
        objective_value = 0
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)

        if objective_value>0:
            sum_rates = 0
            for k in network.each_wk_organizations[wk_idx]:
                for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                        sum_rates = sum_rates +x_vars[k,p].solution_value*network.each_wk_k_weight[wk_idx][k] * network.each_wk_k_u_weight[wk_idx][k][u]

        
        if run_step_number  in [0,50,100,500,900,1000]:
            # we store results
            if objective_value>0:
                purification_schemes_string = ""
                for pur_scheme in network.purification.allowed_purification_schemes:
                    if purification_schemes_string:
                        purification_schemes_string = purification_schemes_string+","+pur_scheme
                    else:
                        purification_schemes_string = pur_scheme
                list_of_flows = []
                for k in network.each_wk_organizations[wk_idx]:
                    for u in network.each_wk_each_k_user_pair_ids[wk_idx][k]: 
                        w=network.each_wk_k_u_weight[wk_idx][k][u] 
                        flow_sum = 0
                        if u not in list_of_flows:
                            list_of_flows.append(u)
                        for p in network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][u]:
                            flow_sum =x_vars[k,p].solution_value+flow_sum
                            with open(network.path_variables_file_path, 'a') as newFile:                                
                                newFileWriter = csv.writer(newFile)
                                newFileWriter.writerow([network.running_path_selection_scheme,
                                                        run_step_number,chromosome_id,k,u,p,
                                                        network.purification.two_qubit_gate_fidelity,
                                                        network.purification.measurement_fidelity,
                                                        network.alpha_value,x_vars[k,p].solution_value,
                                                        wk_idx,sum_rates,
                                                        purification_schemes_string,network.number_of_flows])
                        print("flow %s from %s  has rate %s "%(u,len(list_of_flows),int(flow_sum)))
        
        return objective_value
    
    def CPLEX_swap_scheduling(self,wk_idx,network):
        
        opt_model = cpx.Model(name="swap_scheduling")
        eflow_vars  = {(i,j,k,b): opt_model.continuous_var(lb=0, ub= network.max_edge_capacity,
                        name="eflow_{0}_{1}_{2}_{3}".format(i,j,k,b))
                       for i in network.nodes for j in network.nodes for k in network.nodes for b in network.nodes}

        u_vars  = {(i,j): opt_model.continuous_var(lb=0, ub= 1,
                                  name="u_{0}_{1}".format(i,j)) for i in network.nodes for j in network.nodes}
        
        # We set the capacity of pair of nodes to zero if the pair is not an edge in the network
        for i in network.nodes:
            for j in network.nodes:
                if (i,j) not in network.set_E:
                    network.each_edge_capacity[(i,j)] = 0
        
        # we check incoming and out going flow to each enode
        """we divide the edge capacity by the g(.) for the higherst threshold to count edge level purification.
        
        It is not clear how would we consider end level purification."""
        #print("this is the end level purification flag ",network.end_level_purification_flag)
        if not network.end_level_purification_flag:
            for i in network.nodes:
                for j in network.nodes:
                    #print("we are adding incoming outgoing flow constraint")
                    if i!=j:
                        opt_model.add_constraint(opt_model.sum(network.each_node_q_value[k] *
                        (eflow_vars[i,k,i,j]+eflow_vars[k,j,i,j])/2 for k in network.nodes if k not in [i,j])
                        +(network.each_edge_capacity[(i,j)] 
                         * network.check_edge_exit((i,j))) 
                        - opt_model.sum(
                        (eflow_vars[i,j,i,k]+eflow_vars[i,j,k,j]) for k in network.nodes if k not in [i,j])
                        >=0)
                        #With average round of purification required on paths
#                         opt_model.add_constraint(opt_model.sum(network.each_node_q_value[k] *
#                         (eflow_vars[i,k,i,j]+eflow_vars[k,j,i,j])/2 for k in network.nodes if k not in [i,j])
#                         +(network.each_edge_capacity[(i,j)] 
#                          * network.check_edge_exit((i,j)))
#                         /network.get_required_edge_level_purification_EPR_pairs_all_paths((i,j),i,j,wk_idx) 
#                         - opt_model.sum(
#                         (eflow_vars[i,j,i,k]+eflow_vars[i,j,k,j]) for k in network.nodes if k not in [i,j])
#                         >=0)
        for i in network.nodes:
            for j in network.nodes:
                for k in network.nodes:
                    opt_model.add_constraint(eflow_vars[i,k,i,j]==eflow_vars[k,j,i,j])
                    opt_model.add_constraint(eflow_vars[i,k,i,j]>=0)
        for k in network.each_wk_organizations[wk_idx]:
            for u in network.each_wk_each_k_user_pairs[wk_idx][k]:
                for i in network.nodes:
                    opt_model.add_constraint(eflow_vars[u[0],u[1],u[0],i]==eflow_vars[u[0],u[1],u[1],i])
                    opt_model.add_constraint(eflow_vars[u[0],u[1],u[0],i]==0)



        objective = opt_model.sum((network.each_edge_capacity[(u[0],u[1])] * network.check_edge_exit((u[0],u[1]))
                        + opt_model.sum(network.each_node_q_value[k] * 
                        (eflow_vars[u[0],k,u[0],u[1]]+eflow_vars[k,u[1],u[0],u[1]])/2 for k in network.nodes 
                            if k not in [u[0],u[1]]) )
                            *network.each_wk_k_weight[wk_idx][k] * network.each_wk_k_u_pair_weight[wk_idx][k][u]
                            for k in network.each_wk_organizations[wk_idx]
                            for u in network.each_wk_each_k_user_pairs[wk_idx][k] 
                              )


        # for maximization
        opt_model.maximize(objective)
    #     opt_model.solve()
        #opt_model.print_information()
        #try:
        opt_model.solve()
        #print('docplex.mp.solution',opt_model.solution)

        objective_value = -1
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except ValueError:
            print(ValueError)

        #print("EGR is ",objective_value)
        return objective_value

# %%
