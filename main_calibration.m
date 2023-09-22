clear;
clc;
load_parameters

ineq_matrix_constraints = mean(param.firm_rank)-param.firm_rank; % Imposing probabilities must be positive
options = optimoptions('fmincon','Display','iter');

fmincon(@(x) solve_model_calibration(x,param),0,ineq_matrix_constraints,1/N_f*ones(N_f,1),[],[],[],[],[],options);


