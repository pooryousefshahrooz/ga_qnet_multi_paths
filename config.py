# This is the configuration file for the gradient-based scheme. 
class NetworkConfig(object):
    
  """this class would be used in the reinforcement learning implementation"""  
  version = 'RL_v26'
  project_name = 'machine1_B_20_128X128_LRDR096_Ent01v26'
  testing_results = "results/juniper_path_selection_evaluationsurfnet26.csv"
  training_testing_switching_file = "tf_ckpts/training_testing_switching_file_surfnetv26.txt"
  tf_ckpts = "tf_ckpts/tf_ckpts_surfnetv26"
  method = 'actor_critic'
#   method = 'pure_policy'
  
  model_type = 'Conv'
  scale = 100

  max_step = 80 * scale
  
  initial_learning_rate = 0.0001
  learning_rate_decay_rate = 0.96
  learning_rate_decay_step_multiplier = 5
  learning_rate_decay_step = learning_rate_decay_step_multiplier * scale
  moving_average_decay = 0.9999
  entropy_weight = 0.1
  rl_batch_size = 20
  
 
  save_step = 20 
  max_to_keep = 1000

  Conv2D_out = 128
  Dense_out = 128
#   Conv2D_out = 255
#   Dense_out = 255
#   Conv2D_out = 224
#   Dense_out = 224
  
  optimizer = 'RMSprop'
#   optimizer = 'Adam'
    
  logit_clipping = 10       #10 or 0, = 0 means logit clipping is disabled



  test_traffic_file = 'WK2'
  max_moves = 30            #percentage
  # For pure policy
  baseline = 'avg'          #avg, best
#   baseline = 'best'          #avg, best
  training = False
class Config(NetworkConfig):
  """ this class includes all the experiments setup"""
  list_of_topologies = ["ModifiedSurfnet"]#"ModifiedSurfnet_Eng_7hops_flowsv2","ModifiedSurfnet_Eng_7hops_flowsv2_wider_flows","ModifiedSurfnet_Eng_7hops_flowsv2_wider_flowsFths","ModifiedSurfnetv2.txt"
  multiplexing_value = 1
  work_load_file = 'WK'
  schemes = ["EGR","Hop","EGRSquare","Genetic","RL"]# Set the schemes that you want to evamluate. If you want to evaluate only genetic algorithm, set only Genetic keyword in the list
  schemes = ["Genetic"]  
  set_of_edge_level_Fth= [0.8,0.85,0.90,0.992]#4 distilation strategies
  distilation_scheme_for_huristic = "random"
  num_of_organizations = 10
  min_num_of_paths = 3
  num_of_paths = 3
  candidate_paths_size_for_genetic_alg = 15
  candidate_path_sizes_for_genetic_alg = [3,5,8,10,15]
  each_cut_of_upper_value = {3:8,5:10,8:15,10:15,15:15}
  cut_off_for_path_searching = 5# We use this cut off when we search for the set of paths that a user pair can use
  toplogy_wk_paths_for_each_candidate_size = "results/qVPN_paths_for_each_candidate_size_final.csv"
  toplogy_wk_rl_for_initialization_ga_result_file = "results/gradient_scheme_for_ga_initialization_results_machine1.csv"
  saving_used_paths_flag = False
  q_values = [1]# Different values of q (swap success probability) that we want to evaluate in the experiment
  alpha_values = [0.2]
  two_qubit_gate_fidelity_set = [1.0]
  measurement_fidelity_set =    [0.978]
  toplogy_wk_setup_feasibility_checking = "results/feasibility_checking_engineering_SURFnet_machine30.csv"
  toplogy_wk_scheme_fairness_result = "results/toplogy_wk_scheme_fairness_resultv3.csv"#toplogy_wk_scheme_fairness_resultv3
  max_min_flow_rate_file_path = "results/max_min_flow_rate_file_path_final_non_linear_node3.csv"
  set_equal_Fth = True
  number_of_work_loads = 1
  number_of_flow_set = [50]# the number of flows that we want to maximize EGR for them
  min_flow_rate = 10
  get_flow_rates_info = False
  get_could_have_rate_flag = False
  fractional_max_min_flag = False
  get_sum_rate_falg = True
  equal_weight_flag = True
  non_linear_objective = True
   
#   flow_path_values_results = "results/QVPN_flow_values_machine15_without_rate_constraint.csv"
  flow_path_values_results = "results/ICDCS_final_result_with_rate_constraint_node3.csv"
  """parameters of genetic algorithm"""
  population_sizes = [50]
  max_runs_of_genetic_algorithm = 100
  elit_pop_sizes = [100]
  selection_op_values =[1.0] # values we want to test genetic algorithm
  cross_over_values =[1.0] # values we want to check crossover
  mutation_op_values =[1.0] # values we want to check mutation
  runs_of_genetic_algorithm = 1
  multi_point_mutation_value = 50
  multi_point_crossover_value = 50
  genetic_algorithm_initial_population = "Hop"
  genetic_algorithm_initial_population_rl_epoch_number = 2019
  genetic_algorithm_random_initial_population = 10# only this percent of the inial population is random
  ga_elit_pop_update_step = 10
  ga_selection_prob_update_step = 10
  ga_crossover_mutation_update_step = 10
  ga_crossover_mutation_multi_point_update_step = 10
  ga_only_elit_to_generate_population = True
  dynamic_policy_flag = True
  moving_indivituals_to_next_gen = True
  static_crossover_p = 0.5
  static_mutation_p = 0.2
  static_elit_pop_size = 50
  static_multi_point_crossover_value = 50
  static_multi_point_mutation_value =50
        
  each_flow_minimum_rate_value_file  ="results/each_flow_min_rate_from_min_10_rate_max1k_ratev5.csv"
  previous_run_path_info_file = "results/QVPN_flow_values_experimenting_fairness_with_min_rate_10_with_max_ratev5.csv"
  set_paths_from_file = False
  get_flows_minimum_rate_flag = False
  optimization_problem_with_minimum_rate= False
  min_flow_rates = [10]
  optimization_problem_with_maximum_rate= False
  up_max_rate_value = 1000
  updating_max_rate_flag = False
  max_flow_rate = 10
    
  number_of_training_wks = 1
  workloads_to_test = 1
  wkidx_min_value = 0
  wkidx_max_value = 1
  workloads_with_feasible_solution = [0]
  save_rl_results_for_initialization = False

  toplogy_wk_scheme_result_file = "results/qVPN_final_results_machine33_for_distillation_old_topology.csv"
#   toplogy_wk_scheme_result_file = "results/qVPN_final_results_machine7_for_fairness.csv"
    

def get_config(FLAGS):
  config = Config

  for k, v in FLAGS.__flags.items():
    if hasattr(config, k):
      setattr(config, k, v.value)

  return config
