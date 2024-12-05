#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
from itertools import islice
import random
from itertools import groupby
import time
import math as mt
import csv
import os
import random
import pdb

from os import listdir

from os.path import isfile, join


# In[2]:


class Genetic_algorithm:
    def __init__(self,config):
        self.crossover_p = 1
        self.mutation_p = 1
        self.selection_p = 1
        self.elit_pop_size =10 
        self.multi_point_crossover_value = config.multi_point_crossover_value
        self.multi_point_mutation_value = config.multi_point_mutation_value
        self.number_of_chromosomes = 10# you can set how many chromosome you want to have in each chromosome in config file
        self.max_runs_of_genetic_algorithm = int(config.max_runs_of_genetic_algorithm)# how many generations we want to run genetic algorithm
        self.each_path_replaceable_paths = {}
        self.elit_pop_sizes = config.elit_pop_sizes
        self.cross_over_values=config.cross_over_values
        self.mutation_op_values = config.mutation_op_values
        self.population_sizes = config.population_sizes
        self.random_initial_population = config.genetic_algorithm_random_initial_population
        self.genetic_algorithm_initial_population = config.genetic_algorithm_initial_population
        self.only_elit = config.ga_only_elit_to_generate_population
        self.dynamic_policy_flag = config.dynamic_policy_flag
        
        # to always move the best chromosome to the next population
        self.top_chromosome = []
        self.top_fitness_values = []
        # we use this to check if we are allowed to directly move the indiviual to the next generation
        self.moving_indivitual_to_next_gen = config.moving_indivituals_to_next_gen
        
        # we use the static values in static policy settign for ga
        
        self.static_crossover_p = config.static_crossover_p
        self.static_mutation_p = config.static_mutation_p
        self.static_elit_pop_size = config.static_elit_pop_size
        self.static_multi_point_crossover_value = config.static_multi_point_crossover_value
        self.static_multi_point_mutation_value = config.static_multi_point_mutation_value
            
        
        self.each_fitness_chromosomes = {}
        self.chromosomes = [[1,2,3],[4,5,6]]# an example of chromosomes. Each list in this list is one chromosome
        #each chromosome is a list of paths ids. For example if we have three user pairs and
        # each user can have at most one path, then in the first chromosom [1,2,3] \n
        # path id 1 is for the first user, 2 is for the second and 3 is for third.
        
    def crossover_op(self,chromosome1,chromosome2):
        """this function applies crossover operator to the given chromosomes"""
        new_chromosome1 = []
        new_chromosome2 = []
        points = []
    
        while (len(points)< self.multi_point_crossover_value):
            random_value = random.randint(0,len(chromosome1)-1)
            if random_value not in points:
                points.append(random_value)
        points.sort()
        start_point = 0
        points_indexes = []
        random_value_crossover = random.uniform(0,1)
        if random_value_crossover <=self.crossover_p:
            flag = True
            for random_point in points:
                if flag:
                    flag = False
                    for i in range(start_point,random_point):
                        new_chromosome2.append(chromosome1[i])
                    for i in range(start_point,random_point):
                        new_chromosome1.append(chromosome2[i])
                else:
                    flag = True
                    for i in range(start_point,random_point):
                        new_chromosome2.append(chromosome2[i])
                    for i in range(start_point,random_point):
                        new_chromosome1.append(chromosome1[i])
                start_point = random_point


            for j in range(max(points),len(chromosome2)):
                    new_chromosome1.append(chromosome2[j])
            for j in range(max(points),len(chromosome2)):
                    new_chromosome2.append(chromosome1[j]) 
        else:
            new_chromosome1= chromosome1
            new_chromosome2 = chromosome2
        assert len(new_chromosome1)==len(new_chromosome2),"Two chromosomes don't have the same length!"
        
        return new_chromosome1,new_chromosome2
    
    def mutation_op(self,chromosome):
        """this function applies mutation operator to the given chromosome"""
        
        points = []
        tries = 0
        while(tries< self.multi_point_mutation_value):
            tries+=1
            random_value1 = random.uniform(0,1)
            if random_value1 <=self.mutation_p:
                random_value = random.randint(0,len(chromosome)-1)
                if random_value not in points:
                    points.append(random_value)
            
        points.sort()
        for random_point in points:
            
            possible_values = self.each_path_replaceable_paths[chromosome[random_point]]
            new_value = chromosome[random_point]
            if len(possible_values)>1:
                while(new_value==chromosome[random_point]):
                    new_value = possible_values[random.randint(0,len(possible_values)-1)]
            chromosome[random_point]= new_value
        return chromosome
    def selection_random_chromosomes(self,number):
        assert number in [1,2],"Either ask one or two chromosomes!"
        if number ==2:
            if len(self.elit_population)==0:
                chromosome1= self.chromosomes[random.randint(0,len(self.chromosomes)-1)]
                chromosome2 = self.chromosomes[random.randint(0,len(self.chromosomes)-1)]
                return chromosome1,chromosome2
            elif len(self.elit_population)==1:
                chromosome1= self.elit_population[random.randint(0,len(self.elit_population)-1)]
                chromosome2 = chromosome1
                return chromosome1,chromosome2
            else:
                chromosome1= self.elit_population[random.randint(0,len(self.elit_population)-1)]
                chromosome2 = chromosome1
                #while(chromosome2==chromosome1):
                chromosome2= self.elit_population[random.randint(0,len(self.elit_population)-1)]
               
                return chromosome1,chromosome2
        elif number ==1:
            random_value = random.uniform(0,1)
            if len(self.elit_population)==0:
                chromosome1= self.chromosomes[random.randint(0,len(self.chromosomes)-1)]
            else:
                chromosome1= self.elit_population[random.randint(0,len(self.elit_population)-1)]
            return chromosome1
        
    def select_elit_population(self):
        self.elit_population = []
        # Get Top N fitness values from Records
        N = int(int(self.elit_pop_size/100*len(self.chromosomes))+1)
        self.top_fitness_values = sorted(list(self.each_fitness_chromosomes.keys()),reverse=True)
        all_fitness_values = list(self.each_fitness_chromosomes.keys())
#         print("self.each_fitness_chromosomes")
        self.top_chromosome = self.each_fitness_chromosomes[self.top_fitness_values[0]][0]
        
        for fitness in self.top_fitness_values:
            for chromosome in self.each_fitness_chromosomes[fitness]:
                if len(self.elit_population)<=N and chromosome not in self.elit_population:
                    self.elit_population.append(chromosome)
        for fitness in self.top_fitness_values:
            for chromosome in self.each_fitness_chromosomes[fitness]:
                if len(self.elit_population)<=N:
                    self.elit_population.append(chromosome)
        if not self.only_elit:
            while(len(self.elit_population)<self.number_of_chromosomes):
                random_fitnes = all_fitness_values[random.randint(0,len(all_fitness_values)-1)]
                random_chromosome = self.each_fitness_chromosomes[random_fitnes][0]
                counter = 0
                while(random_chromosome in self.elit_population and counter<400):
                    random_fitnes = all_fitness_values[random.randint(0,len(all_fitness_values)-1)]
                    random_chromosome = self.each_fitness_chromosomes[random_fitnes][0]
                    counter+=1
                self.elit_population.append(random_chromosome)
    def population_gen_op(self):
        self.select_elit_population()
    
        new_population = []
        new_chromosome = []
        for item in self.top_chromosome:
            new_chromosome.append(item)
        new_population.append(new_chromosome)
        
        while(len(new_population)<len(self.chromosomes)):
            new_elected_chromosomes = []
            chromosome1,chromosome2 = self.selection_random_chromosomes(2)
            if self.moving_indivitual_to_next_gen:
                new_chromosome = []
                for item in chromosome1:
                    new_chromosome.append(item)
                new_elected_chromosomes.append(new_chromosome)
                new_chromosome = []
                for item in chromosome2:
                    new_chromosome.append(item)
                new_elected_chromosomes.append(new_chromosome)
            chromosome1,chromosome2 = self.crossover_op(chromosome1,chromosome2)
           
            new_chromosome = []
            for item in chromosome1:
                new_chromosome.append(item)
            new_elected_chromosomes.append(new_chromosome)
            new_chromosome = []
            for item in chromosome2:
                new_chromosome.append(item)
            new_elected_chromosomes.append(new_chromosome)

            chromosome3 = self.selection_random_chromosomes(1)
            if self.moving_indivitual_to_next_gen:
                new_chromosome = []
                for item in chromosome3:
                    new_chromosome.append(item)
                new_elected_chromosomes.append(new_chromosome)
                
            chromosome3 = self.mutation_op(chromosome3)
            new_chromosome = []
            for item in chromosome3:
                new_chromosome.append(item)
            new_elected_chromosomes.append(new_chromosome)
            if new_elected_chromosomes:
                for chromosome in new_elected_chromosomes:
                    new_population.append(chromosome)
        self.chromosomes = new_population
#         print("****** start printing new generation ******* "+"\n \n")
#         for chromosome in self.chromosomes:
#             print("new chromosome ",chromosome)
#         print("****** done printing new generation  ******* "+"\n \n")
        
    def update_operation_probabilities(self,runs_of_genetic_algorithm,config):
        """this function update the selection, crossover, and mutation probabilities to
        trade-of between exploration and exploitation"""
        if not self.dynamic_policy_flag:
            self.crossover_p = self.static_crossover_p
            self.mutation_p = self.static_mutation_p
            self.elit_pop_size = self.static_elit_pop_size
            self.multi_point_crossover_value = self.static_multi_point_crossover_value
            self.multi_point_mutation_value = self.static_multi_point_mutation_value
            
        if self.dynamic_policy_flag and runs_of_genetic_algorithm % config.ga_elit_pop_update_step == config.ga_elit_pop_update_step - 1:            
            self.elit_pop_size =max(self.elit_pop_size-3,20) 
            
        if self.dynamic_policy_flag and runs_of_genetic_algorithm % config.ga_crossover_mutation_update_step == config.ga_crossover_mutation_update_step - 1:            

            self.crossover_p = max(self.crossover_p-0.05,0.5)
            self.mutation_p = max(self.mutation_p-0.05,0.2)
        if self.dynamic_policy_flag and runs_of_genetic_algorithm % config.ga_crossover_mutation_multi_point_update_step == config.ga_crossover_mutation_multi_point_update_step - 1:
            self.multi_point_crossover_value = max(self.multi_point_crossover_value-1,20)
            self.multi_point_mutation_value = max(self.multi_point_mutation_value-1,20)
        
            #             self.selection_p = min(self.selection_p+0.05,0.4)

#             self.mutation_p = max(self.crossover_p-0.03,0.1)
            #print(" ********* Hyper parameters are updated *************** ")
        
    def generate_chromosomes(self,wk_idx,round_of_running,network,config):
        """this function generates a population from finite number of chromosomes"""
        network.ordered_user_pairs = []# we use this list to identify the paths of each user pair in the chromosome later
        # we use the data structure filled in function get_each_user_all_paths() for this
#         network.valid_flows = []
        for k in network.each_wk_organizations[wk_idx]:
            for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                user_pair = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                if user_pair_id not in network.valid_flows and len(network.each_user_pair_id_all_paths[user_pair_id])>0:
                    network.valid_flows.append(user_pair_id)
                network.each_user_organization[user_pair_id] = k
        def get_path_ids(network,user_pair_id,saved_path_edges):
            list_of_path_ids = []
            for path_id in network.each_pair_id_paths[user_pair_id]:
                edges = network.set_of_paths[path_id]
                if edges in saved_path_edges:
                    list_of_path_ids.append(path_id)
            return list_of_path_ids
        def check_path_exist(user_pair_id,saved_path_edges,network):
            for path_id,edges in network.set_of_paths.items():
                if saved_path_edges==edges:
                    if path_id in network.each_pair_id_paths[user_pair_id]:
                        return True,path_id
            return False,-1
        file_scanned_flag = False
        each_k_u_scanned_paths ={}
        for k in network.each_wk_organizations[wk_idx]:# Since we have only one organization, the value of k is always zero
#            for user_pair_id in network.valid_flows:
            for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                user_pair = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                selected_path_ids = []
                all_unique_path_ids = network.each_user_pair_id_all_paths[user_pair_id]
                path_ids = []
                if self.genetic_algorithm_initial_population=="RL":
                    print("******* initializing from RL result *********. ")
                    with open(network.toplogy_wk_rl_for_initialization_ga_result_file, "r") as f:
                        reader = csv.reader( (line.replace('\0','') for line in f) )
                        for line in reader:#0:topology_name,1:wk_idx,2:num_of_paths,
        #                             3:fidelity_threshold_range,4:q_value,
        #                             5:each_link_cost_metric,6:number_of_flows,7:number_of_training_wks,
        #                                 8:egr,9:epoch_number,
        #                                 10:cut_off_for_path_searching,11:k,12:u,13:p
                            #print("line is ",line)
                            if line and len(line)>5:
                                corrupted_line_flag = True
                                try:
                                    topology_name = line[0]
                                    wk_idx = float(line[1])
                                    num_of_paths = int(line[2])
                                    F = round(float(line[3]),3)
                                    q_value = round(float(line[4]),3)
                                    num_flows =int(line[6]) 
                                    number_of_training_wks = int(line[7])
                                    egr = float(line[8])/1000
                                    epoch_number = int(line[9])
                                    k = int(line[11])
                                    flow_id = int(line[12])
                                    path_id = int(line[13])
                                    corrupted_line_flag = False
                                except ValueError:
                                    print(ValueError)
        #                             corrupted_line_flag = True
                                if not corrupted_line_flag:
                                    if (number_of_training_wks == network.number_of_training_wks and 
                                        user_pair_id == flow_id and 
                                        epoch_number==network.genetic_algorithm_initial_population_rl_epoch_number):
#                                         print("we use the result! that is good")
#                                         print("we trained rl with %s wks want to use the result of %s wks RL. this the for %s epochs we want to use %s epochs "%(number_of_training_wks,network.number_of_training_wks,epoch_number,network.genetic_algorithm_initial_population_rl_epoch_number))
                                        try:
                                            network.each_scheme_each_user_pair_id_paths["RL"][flow_id].append(path_id)
                                        except:
                                            try:
                                                network.each_scheme_each_user_pair_id_paths["RL"][flow_id]= [path_id]
                                            except:
                                                try:
                                                    network.each_scheme_each_user_pair_id_paths["RL"]= {}
                                                    network.each_scheme_each_user_pair_id_paths["RL"][flow_id]=[path_id]
                                                except ValueError:
                                                    print("Error",ValueError)
#                                 else:
#                                     print("not a valid line! ",line)
                    for k in network.each_wk_organizations[wk_idx]:
                        for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                            if user_pair_id not in network.each_scheme_each_user_pair_id_paths["RL"]:
                                network.each_scheme_each_user_pair_id_paths["RL"][user_pair_id] = network.each_scheme_each_user_pair_id_paths["Hop"][user_pair_id]
                    path_ids = network.each_scheme_each_user_pair_id_paths["RL"][user_pair_id]
                elif self.genetic_algorithm_initial_population=="EGR":
                    path_ids = network.each_scheme_each_user_pair_id_paths["Hop"][user_pair_id]
                elif self.genetic_algorithm_initial_population=="Hop":
                    path_ids = network.each_scheme_each_user_pair_id_paths["Hop"][user_pair_id]
                if self.dynamic_policy_flag and network.candidate_paths_size_for_genetic_alg == max(config.candidate_path_sizes_for_genetic_alg):
                    path_ids = network.each_user_pair_id_all_paths[user_pair_id]
#                     print("we will use all possible paths for this setup")
                    #time.sleep(2)
                else:
#                     print("we will use the paths used for the upper level for this setup current %s using of %s "%(network.candidate_paths_size_for_genetic_alg,config.each_cut_of_upper_value[network.candidate_paths_size_for_genetic_alg]))
                    #time.sleep(3)
                    """we read and use the paths that have been used for
                    the higer value of cand_paths_size_for_genetic_alg"""
                    we_got_paths_from_file_from_other_setup = False
                    if not file_scanned_flag:
                        
                        print(" for user pair %s we are on #run %s #flow %s |path_cand| %s need to have results for |cand|%s "%(user_pair_id,
                                                                    round_of_running,
                                                                    network.number_of_flows,
                                                                    network.candidate_paths_size_for_genetic_alg,
                                           (config.each_cut_of_upper_value[network.candidate_paths_size_for_genetic_alg])))
#                         time.sleep(10)
                        line_counter = 0

                        with open(network.toplogy_wk_paths_for_each_candidate_size, "r") as f:
#                             number_of_lines = len(f.readlines())
                            reader = csv.reader( (line.replace('\0','') for line in f) )
                            for line in reader:
                                successful_line = True
#                                 try:                    
                                saved_network= (line[0])
                                saved_wk_idx = int(line[1])
                                saved_round_of_running = int(line[2])
                                saved_k = int(line[3])
                                saved_user_pair_id = int(line[4])
                                saved_user_pair = (line[5]).split(":")
                                saved_user_pair = (int(saved_user_pair[0]),int(saved_user_pair[1]))
                                saved_path_id = int(line[6])
                                saved_edges = (line[7])
                                saved_edges = saved_edges.split(",")
                                saved_path_edges =[]
                                for edge in saved_edges:
                                    edge = edge.split(":")
                                    saved_path_edges.append((int(edge[0]),int(edge[1])))
                                #path_existing_flag,current_path_id = check_path_exist(saved_path_edges,network)
                                saved_num_org = int(line[8])
                                saved_num_of_paths = int(line[9])
                                saved_cut_off_for_path_searching =int(line[10])
                                saved_num_of_flows=int(line[12])
                                saved_alpha_value=float(line[13])
                                saved_p2 = float(line[14])
                                saved_eta = float(line[15])
                                saved_len_set_of_edge_level_Fth = int(line[16])
                                successful_line = True
#                                 except ValueError:
#                                     print(ValueError)
                                if successful_line:
                                    if (saved_network==network.topology_name and 
                                        saved_wk_idx==wk_idx and 
                                        saved_round_of_running==round_of_running and 
                                        saved_num_org==config.num_of_organizations and 
                                        saved_num_of_flows ==network.number_of_flows and 
                                        saved_alpha_value == network.alpha_value and 
                                        saved_p2== network.purification.two_qubit_gate_fidelity and 
                                        saved_eta == network.purification.measurement_fidelity 
                                        and 
                                        saved_cut_off_for_path_searching== config.each_cut_of_upper_value[network.candidate_paths_size_for_genetic_alg]):
                                        try:
                                            each_k_u_scanned_paths[saved_k,saved_user_pair].append(saved_path_edges)
                                        except:
                                            each_k_u_scanned_paths[saved_k,saved_user_pair]=[saved_path_edges]
#                                     print("we read %s lines out of %s lines "%(line_counter,123))
#                                     line_counter+=1
#                                      
#                                     print(saved_network ,network.topology_name, saved_wk_idx,wk_idx , 
#                                           saved_round_of_running,round_of_running , saved_k,k , 
#                                         saved_user_pair , user_pair, 
# #                                         saved_user_pair_id == user_pair_id and 
#                                         saved_path_id ,saved_edges, all_unique_path_ids ,                                             
#                                         saved_num_org,config.num_of_organizations ,
#                                         saved_num_of_paths,network.num_of_paths , 
#                                         saved_num_of_flows,network.number_of_flows , 
#                                         saved_alpha_value , network.alpha_value ,
#                                         saved_p2,network.purification.two_qubit_gate_fidelity , 
#                                         saved_eta,network.purification.measurement_fidelity,
#                                         saved_cut_off_for_path_searching, config.each_cut_of_upper_value[network.candidate_paths_size_for_genetic_alg])
#                                     if (
#                                         saved_network==network.topology_name and saved_wk_idx==wk_idx and 
#                                         saved_round_of_running==round_of_running and saved_k==k and 
#                                         saved_user_pair == user_pair and 
# #                                         saved_user_pair_id == user_pair_id and 
#                                            path_existing_flag and
#                                         saved_num_org==config.num_of_organizations and 
#                                         saved_num_of_paths==network.num_of_paths and 
#                                         saved_num_of_flows ==network.number_of_flows and 
#                                         saved_alpha_value == network.alpha_value and 
#                                         saved_p2== network.purification.two_qubit_gate_fidelity and 
#                                         saved_eta == network.purification.measurement_fidelity 
#                                         and saved_cut_off_for_path_searching== config.each_cut_of_upper_value[network.candidate_paths_size_for_genetic_alg]):


#                                         we_got_paths_from_file_from_other_setup  =True
# #                                         if current_path_id in all_unique_path_ids:
# #                                             print("we are adding this path %s for user pair %s "%(current_path_id,user_pair_id))
# #                                         else:
# #                                             print("we do not add this path %s for user pair %s "%(current_path_id,user_pair_id))
#                                         if current_path_id not in path_ids:
#                                             path_ids.append(current_path_id)

                        print("finished scanning the file")
                        file_scanned_flag = True
                    else:
                        print("we have already scanned the file!")
#                     print("each_k_u_scanned_paths ",each_k_u_scanned_paths.keys())
                    path_ids =get_path_ids(network,user_pair_id,each_k_u_scanned_paths[k,user_pair])
                    
#                 print("we have %s paths for user pair %s got maybe from file "%(len(path_ids),user_pair_id))
                random.shuffle(path_ids)
                random.shuffle(all_unique_path_ids)
                permitted_path_ids = []
                permitted_unique_path_ids = []
                for path_id in path_ids:
                    if len(permitted_unique_path_ids)<network.candidate_paths_size_for_genetic_alg:
                        permitted_unique_path_ids.append(path_id)
                for path_id in all_unique_path_ids:
                    if (path_id not in permitted_unique_path_ids and 
                    len(permitted_unique_path_ids)<network.candidate_paths_size_for_genetic_alg):
                        permitted_unique_path_ids.append(path_id)
                for path_id in permitted_unique_path_ids:
                    if path_id not in permitted_path_ids:
                        permitted_path_ids.append(path_id)
                    for path_version in network.purification.each_path_version_numbers[path_id]:
                        if path_version not in permitted_path_ids:
                            permitted_path_ids.append(path_version)
                            
                # we want to know what possible values are for this path to be replaced in mutation
                for p_id1 in permitted_path_ids: 
                    possible_values = []
                    for p_id2 in permitted_path_ids:
                        if p_id1 !=p_id2:
                            possible_values.append(p_id2)
                    self.each_path_replaceable_paths[p_id1] = possible_values
                # these are the paths allowed to be used by this pair
                network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id]=permitted_unique_path_ids
#                 print("for flow id %s we have %s candidate paths "%(user_pair_id,len(permitted_unique_path_ids)))
                
                if self.dynamic_policy_flag and config.saving_used_paths_flag:
                    """we save these paths to use them in the setup for 5 and 10 paths """
                    
#                     print("let save the results at ",network.toplogy_wk_paths_for_each_candidate_size)
                    with open(network.toplogy_wk_paths_for_each_candidate_size, 'a') as newFile:                                
                        newFileWriter = csv.writer(newFile)
                        for p in permitted_unique_path_ids:
                            edges_string = ""
                            edges = network.set_of_paths[p]
                            for edge in edges:
                                if edges_string:
                                    edges_string = edges_string+","+str(edge[0])+":"+str(edge[1])
                                else:
                                    edges_string = str(edge[0])+":"+str(edge[1])
                            user_pair_string = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                            user_pair_string = str(user_pair_string[0])+":"+str(user_pair_string[1])
                            newFileWriter.writerow([network.topology_name,wk_idx,round_of_running,k,user_pair_id
                            ,user_pair_string,p,edges_string,
                            config.num_of_organizations,network.num_of_paths,
                            config.candidate_paths_size_for_genetic_alg,
                            network.each_link_cost_metric,network.number_of_flows,
                            network.alpha_value,
                            network.purification.two_qubit_gate_fidelity,
                            network.purification.measurement_fidelity,
                           len(config.set_of_edge_level_Fth)])


        self.chromosomes = []
#         print("for work load %s we have these flows %s "%(wk_idx,network.valid_flows))
        added_solutions = 0
        
        # We use the paths computed by EGR scheme to initiate the genetic algorth population
        if  self.genetic_algorithm_initial_population in ["EGR","Hop","EGRSquare","RL"]:
            for i in range(int(self.random_initial_population*self.number_of_chromosomes/100)):
                chromosome = []# this is a new chromosome
                for k in network.each_wk_organizations[wk_idx]:# Since we have only one organization, the value of k is always zero
#                         for user_pair_id in network.valid_flows:
                    for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                        path_counter = 0
                        user_pair = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                        path_ids = network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id]
                        selected_path_ids = []
                        random.shuffle(path_ids)
#                         print("for flow %s checking a path %s in permitted paths %s "%(user_pair_id,len(path_ids),len(network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id])))
                        #time.sleep(1)
                        for path_id in path_ids:
                            if path_id in network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id]:
                                selected_path_ids.append(path_id)
                                for path_version in network.purification.each_path_version_numbers[path_id]:
                                    if path_version not in selected_path_ids:
                                        selected_path_ids.append(path_version)
                        path_ids=selected_path_ids
                        for path_id in path_ids:
                            if path_counter <network.num_of_paths:# check how many paths have been selected for each user pair
                                chromosome.append(path_id)
                                path_counter+=1
                        for i in range(path_counter,network.num_of_paths):
#                             print("path_ids to select one random for flow in scheme ",path_ids,user_pair_id,self.genetic_algorithm_initial_population)
#                             time.sleep(0.7)
                            path_id = path_ids[random.randint(0,len(path_ids)-1)]# get one random path
                            chromosome.append(path_id)
                print("len of new added solution ",len(chromosome))
                self.chromosomes.append(chromosome)
                added_solutions+=1
        
        if config.set_paths_from_file:
            print("getting path info ")
            chromosome=[]
            each_flow_paths ={}
            each_path_Fth = {}
            each_flow_Fth = {}
            each_flow_weight = {}
            each_path_edges = {}
            with open(config.previous_run_path_info_file, "r") as f:
                reader = csv.reader( (line.replace('\0','') for line in f) )
                for line in reader:
                    topology_name = (line[0])
                    wk_idx = int(line[1])
                    num_of_organizations = int(line[2])
                    generation_number = int(line[3])
                    k = int(line[4])
                    k_weight = float(line[5])
                    flow_pair = line[6]
                    flow_pair = flow_pair.replace("(","")
                    flow_pair = flow_pair.replace(")","")
                    flow_pair = flow_pair.split(",")
                    flow_pair_src = flow_pair[0]
                    flow_pair_dst = flow_pair[1]
                    flow = (int(flow_pair_src),int(flow_pair_dst))
                    u = int(line[7])
                    flow_weight = float(line[8])
                    p_lenght = int(line[9])
                    x_vars_solution_value = float(line[10])
                    flow_Fth = float(line[11])
                    path_edge_Fth = float(line[12])
                    scheme = line[38]
                    p_id = int(line[42])
                    if (wk_idx in [1] 
                        and generation_number in [0] and 
                        "Genetic" in scheme):
#                             print("we have wk %s k %s u %s flow %s path_edge_Fth %s "%(
#                                             wk_idx,k,u,flow,path_edge_Fth))
                            #time.sleep(0.1)
                            try:
                                each_flow_paths[wk_idx,k,flow].append(p_id)
                            except:
                                each_flow_paths[wk_idx,k,flow]= [p_id]
                            each_path_Fth[p_id] = path_edge_Fth
                            each_flow_Fth[k,flow] = flow_Fth
                            each_flow_weight[k,flow] = flow_weight
                            path_edges = line[44]
                            path_edges = path_edges.split(":")
                            list_of_path_edges = []
                            for edge in path_edges:
                                edge = edge.split("|")
                                list_of_path_edges.append((int(edge[0]),int(edge[1])))
                                
                            each_path_edges[wk_idx,k,p_id] = list_of_path_edges
                            
                            
            for k in network.each_wk_organizations[wk_idx]:
                for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    
                    user_pair = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                    flow_weight=network.each_wk_k_u_weight[wk_idx][k][user_pair_id] 
                    flow_Fth = network.purification.each_wk_k_u_fidelity_threshold[wk_idx][k][user_pair_id]
                    from_file_Fth = each_flow_Fth[k,user_pair]
                    from_file_weight = each_flow_weight[k,user_pair]
                    
#                     print("flows Fth is %s from file it is %s weight %s from file %s "%(flow_Fth,from_file_Fth
#                                                                       ,flow_weight,from_file_weight))
#                     time.sleep(0.1)
                    
                    path_ids = each_flow_paths[wk_idx,k,user_pair]
#                     print("for k %s flow %s we have these path IDs %s "%(k,user_pair,path_ids))
                    selected_path_ids = []
                    extracted_paths = []
                    this_round_path_ids = network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id]
                   
                    for path_id in this_round_path_ids:
                        for path_version in network.purification.each_path_version_numbers[path_id]:
                            for extracted_path_id in path_ids:    
                                path_edges = network.set_of_paths[path_version]
                                path_edges_from_file  =each_path_edges[wk_idx,k,extracted_path_id]
                                if path_edges ==path_edges_from_file:
                                    if each_path_Fth[extracted_path_id]==network.purification.each_path_edge_level_Fth[path_version]:
                                        if path_version not in extracted_paths:
                                            extracted_paths.append(path_version)
                            if path_version not in selected_path_ids:
                                selected_path_ids.append(path_version)
                    path_ids = extracted_paths
#                     for p in path_ids:
#                         if p not in selected_path_ids:
#                             print("this is an error since p %s is not in selected paths %s min %s max %s "%(p,len(selected_path_ids),min(selected_path_ids),max(selected_path_ids)))
#                             time.sleep(2)


                    path_counter = 0
#                     random.shuffle(path_ids)
                    for p in path_ids:
                        chromosome.append(p)
                        path_counter+=1
                    while(path_counter <network.num_of_paths):
                        print("we are adding new path!!!!!!!!!!")
                        if len(path_ids)>0:
                            path_id = path_ids[random.randint(0,len(path_ids)-1)]# get one random path
                        else:
                            path_id==-1
                        chromosome.append(path_id)
                        path_counter+=1
                    
                    
                    
            print("from file len of new added solution ",len(chromosome))
#             print("best chromosome ",chromosome)
#             time.sleep(10)
            chromosome=[18267, 18232, 18236, 18296, 18354, 18294, 18471, 18457, 18462, 18568, 18525, 18574, 18586, 18598, 18613, 18682, 18701, 18696, 18854, 18797, 18795, 18946, 18952, 18954, 19045, 19007, 19040, 19133, 19104, 19110, 19270, 19230, 19220, 19336, 19319, 19324, 19455, 19472, 19484, 19615, 19645, 19612, 19775, 19715, 19752, 19920, 19913, 19921, 20000, 20051, 20010, 20131, 20131, 20134, 20215, 20342, 20234, 20439, 20451, 20451, 20486, 20565, 20549, 20671, 20759, 20712, 20768, 20788, 20844, 20925, 20901, 20862, 21064, 21069, 20964, 21099, 21131, 21080, 21228, 21207, 21170, 21288, 21289, 21245, 21336, 21379, 21331, 21401, 21451, 21413, 21669, 21496, 21474, 21752, 21720, 21843, 21922, 21896, 21906, 22011, 22023, 22003, 22127, 22166, 22130, 22192, 22271, 22271, 22323, 22323, 22390, 22563, 22609, 22569, 22663, 22796, 22722, 22821, 22810, 22813, 22902, 22928, 22925, 23035, 23014, 22971, 23093, 23126, 23136, 23321, 23197, 23215, 23422, 23411, 23456, 23566, 23542, 23522, 23616, 23732, 23726, 23749, 23858, 23863, 24091, 23984, 24009, 24158, 24132, 24196, 24228, 24250, 24224, 24379, 24374, 24389, 24496, 24435, 24496, 24603, 24594, 24590, 24759, 24679, 24727, 24799, 24818, 24840, 24895, 24850, 25006, 25186, 25079, 25191, 25209, 25280, 25244, 25422, 25366, 25414, 25503, 25442, 25475, 25579, 25568, 25556, 25670, 25635, 25707, 25742, 25759, 25729, 25923, 25807, 25966, 26051, 26062, 26015, 26192, 26169, 26228, 26373, 26270, 26300, 26500, 26482, 26506, 26652, 26545, 26602, 26753, 26710, 26780, 26888, 26820, 26884, 26979, 26902, 26918, 27101, 27029, 27020, 27194, 27165, 27263, 27364, 27326, 27442, 27500, 27521, 27496, 27630, 27628, 27553, 27648, 27648, 27760, 27797, 27855, 27852, 27952, 27890, 27939, 28105, 28212, 28209, 28306, 28384, 28233, 28464, 28438, 28459, 28569, 28555, 28561, 28673, 28635, 28654, 28797, 28882, 28809, 28906, 28919, 29016, 29198, 29117, 29229, 29262, 29309, 29247, 29397, 29346, 29387, 29477, 29434, 29472, 29544, 29547, 29558, 29646, 29600, 29630, 29769, 29789, 29740, 30060, 29988, 29927, 30151, 30129, 30144, 30334, 30254, 30308, 30354, 30426, 30366, 30588, 30458, 30508, 30792, 30815, 30738, 30822, 30853, 30833, 30996, 30943, 30919, 31122, 31100, 31207, 31295, 31302, 31265, 31386, 31353, 31357, 31475, 31484, 31459, 31602, 31510, 31564, 31816, 31733, 31791, 31924, 31888, 31864, 32014, 32009, 32014, 32060, 32033, 32077, 32175, 32152, 32140, 32204, 32203, 32209, 32347, 32345, 32301, 32495, 32469, 32465, 32504, 32529, 32555, 32638, 32701, 32676, 32755, 32751, 32747, 32859, 32841, 32824, 33041, 32929, 33023, 33191, 33186, 33140, 33285, 33311, 33292, 33328, 33332, 33394, 33489, 33484, 33450, 33642, 33614, 33602, 33684, 33668, 33663, 33936, 33747, 33794, 34011, 33980, 33979, 34155, 34158, 34149, 34246, 34237, 34226, 34347, 34333, 34339, 34478, 34555, 34553, 34663, 34593, 34650, 34720, 34740, 34731, 34961, 34883, 34856, 35061, 35029, 35068, 35102, 35083, 35180, 35218, 35218, 35255, 35433, 35458, 35393, 35518, 35570, 35540, 35786, 35760, 35808, 35858, 35855, 35922, 35964, 36004, 35964, 36036, 36088, 36076, 36282, 36199, 36168, 36409, 36481, 36436, 36590, 36545, 36582, 36690, 36665, 36648, 36764, 36762, 36799]
            self.chromosomes.append(chromosome)
        
            path_indx = 0
            for k in network.each_wk_organizations[wk_idx]:
                for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    path_ids = chromosome[path_indx:path_indx+network.num_of_paths]
                    path_ids = list(set(path_ids))
                    for p_id1 in path_ids: 
                        possible_values = []
                        for p_id2 in path_ids:
                            if p_id1 !=p_id2:
                                possible_values.append(p_id2)
                        for p in possible_values:
                            try:
                                self.each_path_replaceable_paths[p_id1].append(p)
                            except:
                                self.each_path_replaceable_paths[p_id1]= [p]
                    
                    path_indx = path_indx+network.num_of_paths
        
            
        print("we initialized %s chromosomes "%(len(self.chromosomes)))
        print("and we will add %s more "%(self.number_of_chromosomes-len(self.chromosomes)))
#         time.sleep(1)
        for i in range(self.number_of_chromosomes-int(len(self.chromosomes))):
            chromosome = []# this is a new chromosome
            for k in network.each_wk_organizations[wk_idx]:# Since we have only one organization, the value of k is always zero
#                 for user_pair_id in network.valid_flows:
                for user_pair_id in network.each_wk_each_k_user_pair_ids[wk_idx][k]:
                    path_counter = 0
                    user_pair = network.each_wk_k_id_pair[wk_idx][k][user_pair_id]
                    while(path_counter <network.num_of_paths):# check how many paths have been selected for each user pair
                        path_ids = network.each_wk_k_user_pair_id_permitted_unique_paths[wk_idx,k,user_pair_id]
                        selected_path_ids = []
                        random.shuffle(path_ids)
                        for path_id in path_ids:
                            if len(selected_path_ids)<=network.candidate_paths_size_for_genetic_alg:
                                selected_path_ids.append(path_id)
                                for path_version in network.purification.each_path_version_numbers[path_id]:
                                    if path_version not in selected_path_ids:
                                        selected_path_ids.append(path_version)
                        path_ids=selected_path_ids
                        if len(path_ids)>0:
                            path_id = path_ids[random.randint(0,len(path_ids)-1)]# get one random path
                        else:
                            path_id==-1
                        chromosome.append(path_id)
                        path_counter+=1
            print("one solution was added %s left %s "%(added_solutions,(self.number_of_chromosomes-added_solutions)))
            self.chromosomes.append(chromosome)
            added_solutions+=1
        
                                
        print("we have %s chromosomes "%(len(self.chromosomes)))
#         time.sleep(5)


# In[ ]:





# In[ ]:




