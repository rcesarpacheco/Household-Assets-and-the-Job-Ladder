clear;
clc;
wage_type = "log"; %either log or level
integral_type = "sum"; % either trapezoidal or sum
solve_for_dist = true;
name_file_firm_data  = "Data/firms_out.csv";
share_j2j_flows_out_total_flows_data = 0.27;
load_parameters
options = optimoptions('fmincon','Display','iter');
lb = 0;
fmincon(@(x) solve_model_calibration(x,param,integral_type),lambda_0,[],[],[],[],lb,[],[],options);


