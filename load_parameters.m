%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%             Parameters             %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ga = 2;       % CRRA utility with parameter gamma
rho = 0.05;   % discount rate
firing_rate = 0.175; % poisson rate of getting fired
lambda_0 = 0.51; % poisson rate of job offers when unemployed.
lambda_1 = 0.43; % poisson rate of job offers when employed.
disutility_working = 0.0;
r = 0.03; %interest rate
inc_unemployed =0.1; %income when unemployed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%          Grid Parameters           %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assets
amin = 0;    % borrowing constraint
amax = 100;  % range a
N_a=500;    % number of a grid points

%simulation parameters
tol_vf = 10^(-6); %criterion HJB loop
Delta = 100;   %delta in HJB algorithm

% FIRMS - sigmas and mu for each wage path.

% Import NLSY data from stata that we will use to calibrate model.
% Firms ranked from lowest to highest in dataset
firm_data = readtable(name_file_firm_data);
w_mean_vec = table2array(firm_data(:,"mean_wage"));
sigma_vec = table2array(firm_data(:,"sigma_ou_matlab"));
sigma2_vec = sigma_vec.^2;
firm_size_share = table2array(firm_data(:,'firm_share'));
firm_rank = table2array(firm_data(:,"firm_rank")); %Firms ranked from lowest to highest
thetas = table2array(firm_data(:,"theta")); % thetas for OU process

% wage
N_w=100;         % number of w grid points 
wmin = min(table2array(firm_data(:,'lower_w')));     % Range wages
wmax = max(table2array(firm_data(:,'upper_w')));     % Range wages
% each index of the sigma2_vec and mean_vec is a different firm.
N_f = length(w_mean_vec);

%%
%%--------------------------------------------------
%VARIABLES 
a = linspace(amin,amax,N_a)';  %wealth grid
w = linspace(wmin,wmax,N_w);   % wage grid
da = (amax-amin)/(N_a-1);      
dw = (wmax-wmin)/(N_w-1);
dw2 = dw^2;
aa = repmat(a,1,N_w,N_f);
ww = repmat(w,N_a,1,N_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%           Create Abar              %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Abar should be N_a*N_w*N_f

idx_wages_mean = zeros(N_f,1); % vector of the indexes of the closest element in the wage grid from the wmean
for i=1:N_f
    w_mean = w_mean_vec(i);
    sig2 = sigma2_vec(i);
    if wage_type == "log"
        mu = (thetas(i)*(w_mean - log(w))+sig2/2).*w;        %DRIFT (FROM ITO'S LEMMA)
        s2 = sig2.*w.^2;        %VARIANCE (FROM ITO'S LEMMA)
    else
        mu = thetas(i)*(w_mean - w);        %DRIFT (FROM ITO'S LEMMA)
        s2 = sig2.*ones(1,N_w);        %VARIANCE (FROM ITO'S LEMMA)
    end
      
    chi =  - min(mu,0)/dw + s2/(2*dw2);
    yy =  min(mu,0)/dw - max(mu,0)/dw - s2/dw2;
    zeta = max(mu,0)/dw + s2/(2*dw2);
    
    %This will be the upperdiagonal of the matrix Abar
    updiag=zeros(N_a,1); %This is necessary because of the peculiar way spdiags is defined.
    for j=1:N_w
        updiag=[updiag;repmat(zeta(j),N_a,1)];
    end
    
    %This will be the center diagonal of the matrix Abar
    centdiag=repmat(chi(1)+yy(1),N_a,1);
    for j=2:N_w-1
        centdiag=[centdiag;repmat(yy(j),N_a,1)];
    end
    centdiag=[centdiag;repmat(yy(N_w)+zeta(N_w),N_a,1)];
    
    %This will be the lower diagonal of the matrix Abar
    lowdiag=repmat(chi(2),N_a,1);
    for j=3:N_w
        lowdiag=[lowdiag;repmat(chi(j),N_a,1)];
    end
    
    %Add up the upper, center, and lower diagonal into a sparse matrix
    Abar_block=spdiags(centdiag,0,N_a*N_w,N_a*N_w)+spdiags(lowdiag,-N_a,N_a*N_w,N_a*N_w)+spdiags(updiag,N_a,N_a*N_w,N_a*N_w);
    if i==1
        Abar = [Abar_block,repmat(sparse(N_a*N_w,N_a*N_w),1,N_f-i)];
    
    else
        Abar = [Abar;repmat(sparse(N_a*N_w,N_a*N_w),1,i-1),Abar_block,repmat(sparse(N_a*N_w,N_a*N_w),1,N_f-i)];
    end
    % find closest wage in the grid to the mean_wages. Not used here, will
    % be used later in the VF iteration!
    if wage_type == "log"
        [~,idx_wages_mean(i)] = min(abs(w-exp(w_mean_vec(i))));
    else
         [~,idx_wages_mean(i)] = min(abs(w-w_mean_vec(i)));
    end
end


% rho*V = u(c)-disutility_working + expected valued given new job offer + [A,B]x[V,U]
% rho*U = u(c) + expected valued given new job offer + [C,D]x[V,U]
% A = Abar + Atilde

% A is N_a*N_w*N_f X N_a*N_w*N_f
% B is N_a*N_w*N_f X N_a
% C is N_a X N_a*N_w*N_f
% D is N_a X N_a

C = sparse(N_a,N_a*N_w*N_f);
% B = firing_rate*sparse(repmat(eye(N_a),N_w*N_f,1));

B = firing_rate*sparse(1:N_a*N_w*N_f,repmat(1:N_a,1,N_w*N_f),1);
% save all the variables in a structure 
param.Abar=Abar;
param.B=B;
param.C=C;
param.Delta=Delta;
param.N_a=N_a;
param.N_f=N_f;
param.N_w=N_w;
param.a=a;
param.aa=aa;
param.amax=amax;
param.amin=amin;
param.da=da;
param.disutility_working=disutility_working;
param.dw=dw;
param.firing_rate=firing_rate;
param.firm_rank=firm_rank;
param.firm_size_share=firm_size_share;
param.ga=ga;
param.idx_wages_mean=idx_wages_mean;
param.inc_unemployed=inc_unemployed;
param.lambda_0=lambda_0;
param.lambda_1=lambda_1;
param.r=r;
param.rho=rho;
param.tol_vf=tol_vf;
param.w=w;
param.ww=ww;