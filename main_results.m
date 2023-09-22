clear;
clc;

wage_type = "log";
solve_for_dist = true;
name_file_firm_data  = "firms_out.csv";
slope = 0;
load_parameters
solve_model
frac_unemployed
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