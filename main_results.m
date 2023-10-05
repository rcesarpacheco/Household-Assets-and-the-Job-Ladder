clear;
clc;
wage_type = "log"; %either log or level
integral_type = "sum"; % either trapezoidal or sum
solve_for_dist = true;
lambda_1 = 0.1207;
share_j2j_flows_out_total_flows_data = 0.27;
name_file_firm_data  = "Data/firms_out.csv";
load_parameters
solve_model
frac_unemployed
prob_accepting_offer = sum((U<V_new_job_offer).*repmat(g_u,1,N_f).*repmat(p,N_a,1),'all')/sum(g_u,'all')
%% unconditional wage distribution
mass_each_firm = squeeze(sum(sum(g_e,1),2))*da*dw;
g_e_w = squeeze(sum(g_e,1)*da)./repmat(mass_each_firm',N_w,1); % this assumes w is in levels
plot(w,g_e_w(:,2))
% export to csv

export_g = [w',g_e_w];
writematrix(export_g,'Data/wage_densities_model.csv')
%% indifference curves

idx_current_wealth = 10;
n_means = 1000;
n_sigma = 1000;

v_offer_current_wealth = V_new_job_offer(idx_current_wealth,:);
min_w_mean = min(w_log_mean_vec);
max_w_mean = max(w_log_mean_vec);
min_w_sigma = min(sigma_log_vec);
max_w_sigma = max(sigma_log_vec);
grid_means = linspace(min_w_mean,max_w_mean,n_means);
grid_sigma = linspace(min_w_sigma,max_w_sigma,n_sigma);


[means_q,sigma_q] = meshgrid(grid_means,grid_sigma);
v_interpolate = griddata(w_log_mean_vec,sigma_log_vec,v_offer_current_wealth,means_q,sigma_q,'natural');
mesh(means_q,sigma_q,v_interpolate)
hold on
plot3(w_log_mean_vec,sigma_log_vec,v_offer_current_wealth,"*k")
xlabel('Mean')
ylabel('Sigma')
hold off

figure(2)
hold on
contour(means_q,sigma_q,v_interpolate,10)
plot(w_log_mean_vec,sigma_log_vec,'.')

figure(3)
plot3(w_log_mean_vec,sigma_log_vec,v_offer_current_wealth,"*k")

%%
export_v = [means_q(:),sigma_q(:),v_interpolate(:)];
writematrix(export_v,'indifference_curve.csv')