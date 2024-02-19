clear;
clc;

wage_type = "log";
solve_for_dist = false;
name_file_firm_data  = "Data/firms_out_fake_feb_2024.csv";
slope = 0;
load_parameters
solve_model
save('Workspace\wkspace_sigma_ex2.mat')    
%% indifference curves

idx_current_wealth = 1000;
n_means = 500;
n_sigma = 500;

v_offer_current_wealth = V_new_job_offer(idx_current_wealth,:);
min_w_mean = min(w_mean_vec);
max_w_mean = max(w_mean_vec);
min_w_sigma = min(sigma_vec);
max_w_sigma = max(sigma_vec);
grid_means = linspace(min_w_mean,max_w_mean,n_means);
grid_sigma = linspace(min_w_sigma,max_w_sigma,n_sigma);


[means_q,sigma_q] = meshgrid(grid_means,grid_sigma);
v_interpolate = griddata(w_mean_vec,sigma_vec,v_offer_current_wealth,means_q,sigma_q,'natural');
mesh(means_q,sigma_q,v_interpolate)
hold on
plot3(w_mean_vec,sigma_vec,v_offer_current_wealth,"*k")
xlabel('Mean')
ylabel('Sigma')
hold off

figure(2)
hold on
contour(means_q,sigma_q,v_interpolate,10)
plot(w_mean_vec,sigma_vec,'.')

figure(3)
plot3(w_mean_vec,sigma_vec,v_offer_current_wealth,"*k")

%%
export_v = [means_q(:),sigma_q(:),v_interpolate(:)];
writematrix(export_v,'Results/indifference_curve_rich.csv')

%% export sigma exercise
export_v = [sigma_vec,V_new_job_offer'];

writematrix(export_v,'Results/sigma_exercise_output_matlab2.csv')


