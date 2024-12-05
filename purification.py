#!/usr/bin/env python
# coding: utf-8

# In[1]:


import time


# In[ ]:


class Purification:
    def __init__(self):
        self.each_path_basic_fidelity ={}
        self.each_wk_k_fidelity_threshold = {}
        
        self.global_each_basic_fidelity_target_fidelity_required_EPRs = {}
        self.all_basic_fidelity_target_thresholds = []
        self.each_n_f_purification_result = {}
        self.each_edge_target_fidelity = {}
        self.each_wk_k_u_fidelity_threshold = {}
        self.each_edge_fidelity = {}
        #We use this to track paths that we can not reack Fth even with edge and end level distilation
        self.unable_to_distill_paths = []
        
         # these are the two-gate measurement fidelity and measurement fidleity
        self.two_qubit_gate_fidelity = 1
        self.measurement_fidelity  = 1
        # we want to compute function g only once
        self.function_g_computed  = False
        
        
        
        #we use this to track the function g value for each path 
        self.each_path_edge_level_g ={}
        self.each_path_end_level_g = {}
        
        # by default, the purificaiton scheme is end level. But we update this in network.py file to what is set in config file
        self.purification_schemes = ["end_level"]
        self.allowed_purification_schemes = ["end_level"]
        
        # we use this to set the edge-level fidelity threshold for different distilation strategies
        self.set_of_edge_level_Fth = []
        self.each_path_edge_level_Fth = {}
        # tracking the required distilation at edge e to reach target fidleity Fth
        self.oracle_for_edge_target_fidelity = {}
        # tracking the required distilation at path p to reach target fidleity Fth
        self.oracle_for_target_fidelity = {}
        self.each_path_version_numbers= {}
        self.each_path_id_purificaiton_scheme = {}
    
        self.each_path_flow_target_fidelity = {}
        # we use this to set the function g of paths that can not be distilled
        self.upper_bound_on_distilation = 1
        
        # we track the virtual links in order to not involve them in computing end-toend fidelity
        self.virtual_links = []
        
        # we use this for tracking what target fidleity we need to reach at each edge to satify flow fidelity threshold
        
        self.each_path_target_fidelity_for_edge_dislication = {}
    
    def get_next_fidelity_and_succ_prob_BBPSSW(self,F):
        succ_prob = (F+((1-F)/3))**2 + (2*(1-F)/3)**2
        output_fidelity = (F**2 + ((1-F)/3)**2)/succ_prob

        return output_fidelity, succ_prob

    def get_next_fidelity_and_succ_prob_DEJMPS(self,F1,F2,F3,F4):
        succ_prob = (F1+F2)**2 + (F3+F4)**2
        output_fidelity1 = (F1**2 + F2**2)/succ_prob
        output_fidelity2 = (2*F3*F4)/succ_prob
        output_fidelity3 = (F3**2 + F4**2)/succ_prob
        output_fidelity4 = (2*F1*F2)/succ_prob

        return output_fidelity1, output_fidelity2, output_fidelity3, output_fidelity4, succ_prob

    def get_avg_epr_pairs_BBPSSW(self,F_init,F_target):
        F_curr = F_init
        n_avg = 1.0
        while(F_curr < F_target):
            F_curr,succ_prob = get_next_fidelity_and_succ_prob_BBPSSW(F_curr)
            n_avg = n_avg*(2/succ_prob)
        return  n_avg

    def get_avg_epr_pairs_DEJMPS(self,F_init,F_target):
        F_curr = F_init
        F2 = F3 = F4 = (1-F_curr)/3
        n_avg = 1.0
        while(F_curr < F_target):
            F_curr,F2, F3, F4, succ_prob = self.get_next_fidelity_and_succ_prob_DEJMPS(F_curr, F2, F3, F4)
            n_avg = n_avg*(2/succ_prob)
            
        return  n_avg
    def set_required_EPR_pairs_for_each_path_each_fidelity_threshold(self,wk_idx,path):
        targets = []
        for k,u_target_fidelity in self.each_wk_k_u_fidelity_threshold[wk_idx].items():
            for u,target_fidelity in u_target_fidelity.items():
                if target_fidelity not in targets:
                    targets.append(target_fidelity)
        targets.append(0.6)
        targets.append(0.55)
        targets.append(0.5)
        targets.append(0.0)
        targets.sort()
        counter = 0
        
        path_basic_fidelity = self.each_path_basic_fidelity[path]
        if path_basic_fidelity<0.5:
            n_avg = self.upper_bound_on_distilation
            for target in targets:
                try:
                    self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg 
                except:
                    self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity]={}
                    self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg 
                try:
                    self.oracle_for_target_fidelity[path][target] = n_avg
                except:
                    self.oracle_for_target_fidelity[path] = {}
                    self.oracle_for_target_fidelity[path][target] = n_avg
        else:
            try:
                if path_basic_fidelity in self.global_each_basic_fidelity_target_fidelity_required_EPRs:
                    for target in targets:
                        n_avg = self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target]
                        try:
                            self.oracle_for_target_fidelity[path][target] = n_avg
                        except:
                            self.oracle_for_target_fidelity[path] = {}
                            self.oracle_for_target_fidelity[path][target] = n_avg
                else:
                    for target in targets:
                        n_avg = self.get_avg_epr_pairs_DEJMPS(path_basic_fidelity ,target)
                        try:
                            self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg 
                        except:
                            self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity]={}
                            self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg 
                        try:
                            self.oracle_for_target_fidelity[path][target] = n_avg
                        except:
                            self.oracle_for_target_fidelity[path] = {}
                            self.oracle_for_target_fidelity[path][target] = n_avg
            except:
                for target in targets:
                    n_avg  = self.get_avg_epr_pairs_DEJMPS(path_basic_fidelity ,target)
                    try:
                        self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg 
                    except:
                        self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity]={}
                        self.global_each_basic_fidelity_target_fidelity_required_EPRs[path_basic_fidelity][target] =n_avg                     
                    try:
                        self.oracle_for_target_fidelity[path][target] = n_avg
                    except:
                        self.oracle_for_target_fidelity[path] = {}
                        self.oracle_for_target_fidelity[path][target] = n_avg
            #print("this is important ",self.topology_name,self.oracle_for_target_fidelity,self.each_path_basic_fidelity)
    def check_perfect_edge_reach_this_target(self,path_edges,new_target):
        
        F_product = (4*1-1)/3 
        for edge in path_edges[1:]:
            if edge not in self.virtual_links and (edge[1],edge[0]) not in self.virtual_links:
                F_product  = F_product*(4*1-1)/3
        N = len(path_edges)+1
        p1 = self.two_qubit_gate_fidelity
        p2 = self.measurement_fidelity
        F_final = (1/4)+(3/4)*((p1*(((4*(p2**2))-1)/3))**(N-1))*(F_product)
        if new_target<= round(F_final,3):
            return True
        else:
            return False
        
    def set_required_EPR_pairs_for_distilation(self,wk_idx,network):
        """this function computes the g value for edge level and end level distilation for each path given dis. strategy"""
        if self.function_g_computed:
            pass
        else:
            covered_path_counter = 0
#             for k,user_pair_ids in network.each_wk_each_k_user_pair_ids[wk_idx].items():
#                 for user_pair_id in user_pair_ids:
#                     print("these are the user pairs %s for org %s "%(network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k],k))
#                     path_ids = network.each_wk_each_k_each_user_pair_id_paths[wk_idx][k][user_pair_id]
            for path_id,path_edges in network.set_of_paths.items():
                #print("computing g function for path %s out of %s"%(covered_path_counter,len(list(network.set_of_paths.keys()))))
                covered_path_counter+=1
                # we get the flow-specific fidelity threshold
                F_threshold = self.each_path_flow_target_fidelity[path_id]
#                 p_lenght = len(path_edges)
                p_lenght = self.get_path_actual_length(path_edges,network.set_of_virtual_links)
                F_threshold_when_only_edge_level = round((3*(4/3*F_threshold-1/3)**(1/p_lenght)+1)/4,3)
                #since all edges have the same fideity, the first edge fidelity gives us the whole edges fidelity
                edge_fidelity = self.each_edge_fidelity[path_edges[0]]

                # e2e fidleity of path without any edge level distilation
                path_basic_fidelity = self.compute_e2e_fidelity(edge_fidelity,path_edges,network.set_of_virtual_links)

                
                #distilation strategy (the one associated to this path) indicates what is the edge-level fidelity threshold
                path_dis_strategy_edge_Fth = self.each_path_edge_level_Fth[path_id]
                """here is the hack. We distill even thogh we do not need it. 
                The reason is that the distilation strategy is blind to the fidelity threshold of flows. 
                If distilation strategy was not blind, we will uncomment the next command"""
#                 path_edge_Fth = min(path_dis_strategy_edge_Fth,F_threshold_when_only_edge_level)
                path_edge_Fth = path_dis_strategy_edge_Fth
                if edge_fidelity>=path_edge_Fth:#we do not need sitilation
                    edge_level_g = 1
                    path_fidelity_after_edge_level_dis = path_basic_fidelity
                else:# we need distilation to bring the edge fidelity to the edge-level fidelity threshold
                    edge_level_g = self.get_avg_epr_pairs_DEJMPS(edge_fidelity ,path_edge_Fth)
                    path_fidelity_after_edge_level_dis = self.compute_e2e_fidelity(path_edge_Fth,path_edges,network.set_of_virtual_links)


                if path_fidelity_after_edge_level_dis>= F_threshold:# we do not need end-level distilation
                    end_level_g = 1
                else:# we need distilation to reach the Fth
                    if path_fidelity_after_edge_level_dis<0.5:
                        end_level_g = self.upper_bound_on_distilation
                        self.unable_to_distill_paths.append(path_id)
                        
                    else:
                        end_level_g =self.get_avg_epr_pairs_DEJMPS(path_fidelity_after_edge_level_dis ,F_threshold)
                self.each_path_edge_level_g[path_id] = edge_level_g
                self.each_path_end_level_g[path_id] = end_level_g
#                 path_g_function_value = edge_level_g+end_level_g# this is wrong!
                path_g_function_value = edge_level_g*end_level_g
                try:
                    self.oracle_for_target_fidelity[path_id][F_threshold] = path_g_function_value
                except:
                    self.oracle_for_target_fidelity[path_id] = {}
                    self.oracle_for_target_fidelity[path_id][F_threshold] = path_g_function_value
#                 print("for path id %s we have Fth %s and g is %s "%(path_id,F_threshold,(edge_level_g+end_level_g)))
#                 if self.each_path_id_purificaiton_scheme[path_id]=="end_level":
#                     self.set_required_EPR_pairs_for_each_path_each_fidelity_threshold(wk_idx,path_id)
#                 if self.each_path_id_purificaiton_scheme[path_id]=="edge_level":
#                     self.set_required_edge_level_purification_EPR_pairs(wk_idx,network,path_id,path_edges)

            self.function_g_computed = True
                
    def check_path_distilability(self,path_id):
        """this function checks if this path has been set as a path that can not be distilled"""
        if path_id in self.unable_to_distill_paths:
            return False
        else:
            return True
        
    
    def set_required_edge_level_purification_EPR_pairs(self,wk_idx,network,path_id,path_edges):
        """this function, for each path, and for each edge of that path, 
        computes the function g to reach the flow fidelity of the flow that uses that path"""
        F_threshold = self.each_path_flow_target_fidelity[path_id]
        p_lenght = len(path_edges)
        new_target = round((3*(4/3*F_threshold-1/3)**(1/p_lenght)+1)/4,3)
        
        end_level_g = self.each_path_edge_level_g[path_id][new_target]
        edge_level_g = self.each_path_edge_level_g[path_id][new_target]
        self.oracle_for_edge_target_fidelity[edge][new_target] =end_level_g
        #print("going to compute edge level destilation ********")
        
        self.each_path_target_fidelity_for_edge_dislication[path_id] = new_target
        if self.check_perfect_edge_reach_this_target(path_edges,new_target):
            for edge in path_edges:
                edge_basic_fidelity = self.each_edge_fidelity[edge]
                try:
                    if new_target in self.oracle_for_edge_target_fidelity[edge]:
                        pass
                    else:
                        n_avg = self.get_avg_epr_pairs_DEJMPS(edge_basic_fidelity ,new_target)
                        try:
                            self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
                        except:
                            self.oracle_for_edge_target_fidelity[edge] = {}
                            self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
                except:

                    new_target = round((3*(4/3*F_threshold-1/3)**(1/p_lenght)+1)/4,3)
                    edge_basic_fidelity = self.each_edge_fidelity[edge]
                    n_avg = self.get_avg_epr_pairs_DEJMPS(edge_basic_fidelity ,new_target)
                    try:
                        self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
                    except:
                        self.oracle_for_edge_target_fidelity[edge] ={}
                        self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
        else:#this means we can not reach this target fidelity even if we purify all edges to 1!
            #so we set the value of g to a value higher than its capacity!
            n_avg = self.upper_bound_on_distilation
            for edge in path_edges:
                try:
                    self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
                except:
                    self.oracle_for_edge_target_fidelity[edge] ={}
                    self.oracle_for_edge_target_fidelity[edge][new_target] = n_avg
    def get_required_purification_EPR_pairs(self,wk_idx,k,flow,edge,p,threshold):
#         print("we get g value for wk %s org %s flow %s path %s edge %s threshold %s "%(wk_idx,k,flow,p,edge,threshold))
#         print("self.oracle_for_target_fidelity",self.oracle_for_target_fidelity)
#         print("self.oracle_for_target_fidelity[p]",self.oracle_for_target_fidelity[p])
        try:
            return self.oracle_for_target_fidelity[p][threshold]
        except:
            print("self.oracle_for_target_fidelity[p] ",self.oracle_for_target_fidelity[p])
#         if self.each_path_id_purificaiton_scheme[p]=="end_level":
#             return self.oracle_for_target_fidelity[p][threshold]
#         elif self.each_path_id_purificaiton_scheme[p]=="edge_level":
#             new_target = round(self.each_path_target_fidelity_for_edge_dislication[p],3)
#             return self.oracle_for_edge_target_fidelity[edge][new_target]
#         else:
#             print("we have not set any purificaiton scheme for this path!!!!!!")
#             import pdb
#             pdb.set_trace()
#             return False
    def get_each_wk_k_u_threshold(self,wk_idx,k,u):
        return self.each_wk_k_u_fidelity_threshold[wk_idx][k][u]
    
    
    def get_path_actual_length(self,path_edges,virtual_links):
        path_length = 0
        for edge in path_edges:
            if edge not in virtual_links and (edge[1],edge[0]) not in virtual_links:
                path_length+=1
        if path_length==0:
            path_length=1
        return path_length
    def compute_e2e_fidelity(self,edge_F,path_edges,virtual_links):
        if path_edges:
            F_product = (4*edge_F-1)/3 
            for edge in path_edges[1:]:
                if edge not in virtual_links and (edge[1],edge[0]) not in virtual_links:
                    F_product  = F_product*(((4*edge_F)-1)/3)
        else:
            print("Error")
            return 1.0
        N = len(path_edges)
        p1 = self.two_qubit_gate_fidelity
        p2 = self.measurement_fidelity
        F_final = (1/4)+(3/4)*((p1*(((4*(p2**2))-1)/3))**(N-1))*(F_product)
        
        return round(F_final,3)
    
    def get_fidelity(self,path_edges,virtual_links):
        if path_edges:
            F_product = (4*self.each_edge_fidelity[path_edges[0]]-1)/3 
            for edge in path_edges[1:]:
                if edge not in virtual_links and (edge[1],edge[0]) not in virtual_links:
                    F_product  = F_product*(((4*self.each_edge_fidelity[edge])-1)/3)

        else:
            print("Error")
            return 1.0
        N = len(path_edges)
        p1 = self.two_qubit_gate_fidelity
        p2 = self.measurement_fidelity
        F_final = (1/4)+(3/4)*((p1*(((4*(p2**2))-1)/3))**(N-1))*(F_product)
        
        return round(F_final,3)
    def get_required_edge_level_purification_EPR_pairs_all_paths(self,edge,source,destination,wk_idx):
        #print("checking pair %s in edges %s"%(edge,self.set_E))
        if edge not in self.set_E:
            return 1
        else:
            longest_p_lenght = self.find_longest_path(source,destination)
            #print("the length of the longest path is ",len(longest_p_lenght))
            #print("longest_p_lenght",longest_p_lenght)
            return 1
            try:
                if new_target in self.each_edge_target_fidelity[edge]:
                    return self.each_edge_target_fidelity[edge][new_target]
                else:
                    n_avg = self.get_avg_epr_pairs_DEJMPS(edge_basic_fidelity ,new_target)
                    try:
                        self.each_edge_target_fidelity[edge][new_target] = n_avg
                    except:
                        self.each_edge_target_fidelity[edge] = {}
                        self.each_edge_target_fidelity[edge][new_target] = n_avg
                    return n_avg
            except:
                if longest_p_lenght==0:
                    new_target = self.each_edge_fidelity[edge]
                else:
                    new_target = (3*(4/3*max_F_threshold-1/3)**(1/longest_p_lenght)+1)/4

                edge_basic_fidelity = self.each_edge_fidelity[edge]
                n_avg = self.get_avg_epr_pairs_DEJMPS(edge_basic_fidelity ,new_target)
                try:
                    self.each_edge_target_fidelity[edge][new_target] = n_avg
                except:
                    self.each_edge_target_fidelity[edge] ={}
                    self.each_edge_target_fidelity[edge][new_target] = n_avg
            return n_avg    


# In[2]:


# def get_next_fidelity_and_succ_prob_DEJMPS(F1,F2,F3,F4):
#     succ_prob = (F1+F2)**2 + (F3+F4)**2
#     output_fidelity1 = (F1**2 + F2**2)/succ_prob
#     output_fidelity2 = (2*F3*F4)/succ_prob
#     output_fidelity3 = (F3**2 + F4**2)/succ_prob
#     output_fidelity4 = (2*F1*F2)/succ_prob

#     return output_fidelity1, output_fidelity2, output_fidelity3, output_fidelity4, succ_prob


# def get_avg_epr_pairs_DEJMPS(F_init,F_target):
#     F_curr = F_init
#     F2 = F3 = F4 = (1-F_curr)/3
#     n_avg = 1.0
#     while(F_curr < F_target):
#         F_curr,F2, F3, F4, succ_prob = get_next_fidelity_and_succ_prob_DEJMPS(F_curr, F2, F3, F4)
#         n_avg = n_avg*(2/succ_prob)

#     return  n_avg


# In[7]:


# get_avg_epr_pairs_DEJMPS(0.65,0.85)


# In[ ]:




