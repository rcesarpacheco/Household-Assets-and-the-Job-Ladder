clear;
clc;
wage_type = "log"; %either log or level
integral_type = "sum"; % either trapezoidal or sum
solve_for_dist = true;
name_file_firm_data  = "Data/firms_out.csv";

load_parameters
options = optimoptions('fmincon','Display','iter');
lb = zeros(1,N_f);
Aeq = ones(1,N_f);
beq = 1;
fmincon(@(x) solve_model_calibration(x,param,integral_type),1/N_f*ones(1,N_f),[],[],Aeq,beq,lb,[],[],options);


