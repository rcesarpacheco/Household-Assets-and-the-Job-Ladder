%Finite difference approximation of the partial derivatives
% employed
constant = 1/N_f*(1-slope*sum(firm_rank));
p = slope*firm_rank'+constant;
Vaf = zeros(N_a,N_w,N_f); % forward difference            
Vab = zeros(N_a,N_w,N_f); % backward difference

%unemployed
Uaf = zeros(N_a,1); % forward difference      
Uab = zeros(N_a,1); % backward difference
c_employed = zeros(N_a,N_w,N_f); % consumption decision when employed

V0 = (w + r.*aa).^(1-ga)/(1-ga)/rho;
U0 = (inc_unemployed + r*a).^(1-ga)/(1-ga)/rho;
V = V0;
exp_value_job_offer_employed = zeros(N_a,N_w,N_f);
U = U0;
V_new_job_offer = zeros(N_a,N_f);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
% Main loop Value Function Iteration %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err = 2*tol_vf;
while err>tol_vf
    V_old = V;
    U_old = U;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                    %
    %   Update Value Function Employed   %
    %                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % forward difference
    Vaf(1:N_a-1,:,:) = (V(2:N_a,:,:)-V(1:N_a-1,:,:))/da;
    Vaf(N_a,:,:) = repmat((w + r.*amax).^(-ga),1,1,N_f); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    Vab(2:N_a,:,:) = (V(2:N_a,:,:)-V(1:N_a-1,:,:))/da;
    Vab(1,:,:) = repmat((w + r.*amin).^(-ga),1,1,N_f);  %state constraint boundary condition

    %consumption and savings with forward difference
    cf = Vaf.^(-1/ga);
    sf = ww + r.*aa - cf;
    %consumption and savings with backward difference
    cb = Vab.^(-1/ga);
    sb = ww + r.*aa - cb;
    %consumption and derivative of value function at steady state. Savings
    %is equal to 0
    c0 = ww + r.*aa;
    Va0 = c0.^(-ga);

    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
  
    % value of receiving a offer from each firm, for each possible wealth
    for i=1:N_f
        V_new_job_offer(:,i)=V(:,idx_wages_mean(i),i);    
    end
    % calculate expected value when a job offer arrives. try to vectorize
    % this later to increase speed
    for i=1:N_a
        for j=1:N_w
            for k=1:N_f
                exp_value_job_offer_employed(i,j,k)=sum(max(repmat(V(i,j,k),1,N_f),V_new_job_offer(i,:)).*p);
            end
        end
    end
    
    Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
    c_employed = Va_Upwind.^(-1/ga);
    s_employed = ww+r*aa-c_employed;
    u_employed = c_employed.^(1-ga)/(1-ga)-disutility_working;

    %CONSTRUCT MATRIX A
    for i=1:N_f
        X = - min(sb(:,:,i),0)/da;
        Y = - max(sf(:,:,i),0)/da + min(sb(:,:,i),0)/da-lambda_1-firing_rate;
        Z = max(sf(:,:,i),0)/da;
        
        updiag=[0]; %This is needed because of the peculiarity of spdiags.
        for j=1:N_w
            updiag=[updiag;Z(1:N_a-1,j);0];
        end
        
        centdiag=reshape(Y,N_a*N_w,1);
        
        lowdiag=X(2:N_a,1);
        for j=2:N_w
            lowdiag=[lowdiag;0;X(2:N_a,j)];
        end
        Atilde_block = spdiags(centdiag,0,N_a*N_w,N_a*N_w)+spdiags([updiag;0],1,N_a*N_w,N_a*N_w)+ ...
        spdiags([lowdiag;0],-1,N_a*N_w,N_a*N_w);
        if i==1
            Atilde = [Atilde_block,repmat(sparse(N_a*N_w,N_a*N_w),1,N_f-i)];
    
        else
            Atilde = [Atilde;repmat(sparse(N_a*N_w,N_a*N_w),1,i-1),Atilde_block,repmat(sparse(N_a*N_w,N_a*N_w),1,N_f-i)];
      end
    end
       
    A = Atilde + Abar;

    u_stacked = reshape(u_employed,N_a*N_w*N_f,1);
    V_stacked = reshape(V,N_a*N_w*N_f,1);
    exp_value_job_offer_employed_stacked = reshape(exp_value_job_offer_employed,N_a*N_w*N_f,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                    %
    %  Update Value Function Unemployed  %
    %                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % forward difference
    Uaf(1:N_a-1) = (U(2:N_a)-U(1:N_a-1))/da;
    Uaf(N_a) = (inc_unemployed + r.*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    Uab(2:N_a) = (U(2:N_a)-U(1:N_a-1))/da;
    Uab(1) = (inc_unemployed + r.*amin).^(-ga);  %state constraint boundary condition

    %consumption and savings with forward difference
    cf = Uaf.^(-1/ga);
    sf = inc_unemployed + r.*a - cf;
    %consumption and savings with backward difference
    cb = Uab.^(-1/ga);
    sb = inc_unemployed + r.*a - cb;
    %consumption and derivative of value function at steady state
    c0 = inc_unemployed + r.*a;
    Ua0 = c0.^(-ga);

    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
  
    Ua_Upwind = Uaf.*If + Uab.*Ib + Ua0.*I0; 

    c_unemployed = Ua_Upwind.^(-1/ga);
    s_unemployed = inc_unemployed + r.*a - c_unemployed;
    u_unemployed = c_unemployed.^(1-ga)/(1-ga);
    
    %CONSTRUCT MATRIX D for unemployment
    X = - min(sb,0)/da;
    Y = - max(sf,0)/da + min(sb,0)/da-lambda_0;
    Z = max(sf,0)/da;
    
    updiag=[0;Z(1:N_a-1)]; %need to put a fucking 0 because the spdiags is the most stupid function in matlab
    centdiag=Y;
    lowdiag = [X(2:N_a);0]; 

    D=spdiags(centdiag,0,N_a,N_a)+spdiags(updiag,1,N_a,N_a)+spdiags(lowdiag,-1,N_a,N_a);
    
    % Expected value when job offer arrives and unemployed
    exp_value_job_offer_unemployed = sum(max(repmat(U,1,N_f),V_new_job_offer).*p,2);
    
    b_u = U/Delta+u_unemployed+lambda_0*exp_value_job_offer_unemployed;

    M = (1/Delta + rho)*speye(N_a*N_w*N_f+N_a) - [A,B;C,D];
    b = [u_stacked + V_stacked/Delta+lambda_1*exp_value_job_offer_employed_stacked;b_u];
    
    sol_system = M\b; %SOLVE SYSTEM OF EQUATIONS
    V_stacked = sol_system(1:N_a*N_w*N_f);
    U  = sol_system(N_a*N_w*N_f+1:end);
    V = reshape(V_stacked,N_a,N_w,N_f);
    dist_V = max(abs(V-V_old),[],"all");
    dist_U =  max(abs(U-U_old));
    disp(dist_V);
    disp(dist_U);
    err = max(dist_V,dist_U);
end
% update V_new_job_offer last time
for i=1:N_f
    V_new_job_offer(:,i)=V(:,idx_wages_mean(i),i);    
end

if solve_for_dist
    %% solve for stationary distribution
   
    % A_aux: additional matrix to A that captures flows from Employment to
    % Employment
    % to employment.
    
    % D_Aux: aditional matrix to D that captures flows from Unemployment to
    % Unemployment

    % C_aux: additional matrix to C that captures flows from Unemployment
    % to employment.

    % D_aux
    I_offer_empl_better_unp = V_new_job_offer>repmat(U,1,N_f);
    D_aux = spdiags(lambda_0*sum((1-I_offer_empl_better_unp).*repmat(p,N_a,1),2),0,N_a,N_a);
    elements_C = I_offer_empl_better_unp.*repmat(p,N_a,1)*lambda_0;
    %% C_aux
    row = repelem(1:N_a,N_f);
    sz = [N_a,N_w,N_f];
    
    w_array_idx = repmat(idx_wages_mean',1,N_a);
    f_array_idx = repmat(1:N_f,1,N_a);
    
    columns = sub2ind(sz,row,w_array_idx,f_array_idx);
    C_aux = sparse(row,columns,reshape(elements_C',N_a*N_f,1),N_a,N_a*N_w*N_f);
    
    %% A_aux
    diag_A_aux=zeros(N_a*N_w*N_f,1);
    elements_A_Aux_2 = zeros(N_a*N_w*N_f,N_f);
    idx=1;
    for k=1:N_f
      for j=1:N_w
          for i=1:N_a
            ind_current_job_better_than_offer = V(i,j,k)>=V_new_job_offer(i,:);
            diag_A_aux(idx)=lambda_1*sum(ind_current_job_better_than_offer.*p);
            elements_A_Aux_2(idx,:) = lambda_1*(1-ind_current_job_better_than_offer).*p;
            idx=idx+1;
          end
      end
    end      
    %%
    row_to_firm = repmat(1:N_a,1,N_f*N_w);
    row_to_firm = repelem(row_to_firm,N_f);
    w_array_idx = repmat(idx_wages_mean',1,N_a*N_w*N_f);
    f_array_idx = repmat(1:N_f,1,N_a*N_w*N_f);
    
    columns = sub2ind(sz,row_to_firm,w_array_idx,f_array_idx);
    row = repelem(1:N_a*N_f*N_w,N_f);
    A_aux_2 = sparse(row,columns,reshape(elements_A_Aux_2',N_a*N_f*N_w*N_f,1),N_a*N_w*N_f,N_a*N_w*N_f);
    
    %%
    
    A_aux = A_aux_2+spdiags(diag_A_aux,0,N_a*N_w*N_f,N_a*N_w*N_f);
    
    A_dist = [A+A_aux,B;C+C_aux,D+D_aux];
    % V_iteration = [reshape(V,N_a*N_f*N_w,1);U];
    % V_dist = ([u_stacked;u_unemployed]+A_dist*V_iteration)/rho;
    % test = max(abs(V_iteration-V_dist)) % if this number is big, i think something is wrong. If it's small, it still might be wrong.
    
    
    
    AT = A_dist';
    b = zeros(N_a*N_w*N_f+N_a,1);
    
    % %need to fix one value, otherwise matrix is singular
    i_fix = 1;
    b(i_fix)=.1;
    row = [zeros(1,i_fix-1),1,zeros(1,N_a*N_w*N_f+N_a-i_fix)];
    AT(i_fix,:) = row;
    % b(1)=1;
    % row = [da*dw*ones(1,N_a*N_w*N_f),ones(1,N_a)*da];
    % AT(1,:) = row;
    
    %Solve linear system
    disp('solving g')
    gg = AT\b;
    disp('done solving g')
    g_u = gg(N_a*N_w*N_f+1:end);
    g_e = reshape(gg(1:N_a*N_w*N_f),N_a,N_w,N_f);
    sum_g = sum(g_e,"all")*da*dw+sum(g_u)*da; %normalizing so that the integral be equal to one
    g_u=g_u/sum_g;
    g_e=g_e/sum_g;

    %% fraction of people employed in each firm
    frac_employed = sum(g_e,'all')*da*dw;
    frac_employed_per_firm = squeeze(sum(sum(g_e,1),2)*da*dw)/frac_employed;
    frac_unemployed = sum(g_u)*da;

    loss = sum((frac_employed_per_firm-firm_size_share).^2);
end




% %% checking cross
% 
% diff_v_offers = zeros(N_a,N_f*(N_f-1)/2);
% same_sign = zeros(N_f*(N_f-1)/2,1);
% idx=1;
% for i=1:N_f
%     for j=(i+1):N_f
%         current_column=V_new_job_offer(:,i)-V_new_job_offer(:,j);
%         diff_v_offers(:,idx) = current_column;
%         same_sign(idx) = ~any(diff(sign(current_column(current_column~=0))));
%         idx=idx+1;
%     end
% end
% 
% %% checking preferences over firms for each asset level
% 
% preferences = zeros(N_a,N_f);
% 
% for i=1:N_a
%     current_asset = V_new_job_offer(i,:);
%     [~,preferences(i,:)] = sort(current_asset,'descend');
% end
% 
% 
