clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% data loading
load('lung_mutational_data.mat') 
lung_EGFR = lung_data.X5039;
lung_KRAS = lung_data.X82;
lung_TP53 = lung_data.X1;

lung_data = 100*[lung_EGFR, lung_KRAS, lung_TP53];

obs_data_EGFR = lung_data(1:100,1);
obs_data_KRAS = lung_data(1:100,2);
obs_data_TP53 = lung_data(1:100,3);

EGFR_Signal = obs_data_EGFR';
KRAS_Signal = obs_data_KRAS';
TP53_Signal = obs_data_TP53';

%% Multivariate ARMA_1_1 of type: A(t)x(1:t) = B(t)u(1:t) ---------------
iTi = 2000; 
corr = zeros(1,iTi);
m_u = zeros(iTi,1); 
nu_square = 1.78; 
H = 0.58; % 0.7/// 0.9 % Hurst parameter

%% Computation of the correlation matrix using the Browmian Law ---------------
for t=1:iTi 
    corr(t) = 0.5.*(((abs(t+1)).^(2.*H))-2*((abs(t)).^(2.*H))+((abs(t-1)).^(2.*H)));
end
corr_mat = toeplitz(corr);
sigma = nu_square*corr_mat; 
state_x = mvnrnd(m_u,sigma); 

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% parameters ----------------------
N = 2000; % number of Particles 
iT = length(EGFR_Signal);
muMult = zeros(1,iT);
statVecPart_ini = zeros(iT,N);
statVecPart = zeros(iT,N);
errorRatio = zeros(iT,N);
statexPart = zeros(iT,N);
statexPartArma = zeros(iT,N); %%£
invCovMat = zeros(iT,N);
variance = zeros(iT,N);
invCovMatArma = zeros(iT,N);
obsVec = zeros(1,iT);
ObsMat = zeros(1,iT);
gamFunc = zeros(1,iT);
inv_covMat = zeros(iT,N);
vESS = zeros(1,iT);
vESSArma = zeros(1,iT);
vEps = randn(iT,1);
sisrPartState = zeros(1,iT);
EGFR_Model = zeros(1,iT);
KRAS_Model = zeros(1,iT);
TP53_Model = zeros(1,iT);
sisrPartStateArma = zeros(1,iT);
Weights_pred = ones(1,N)*(1/N); %initial weight dist
logWeight0 = log(Weights_pred);
MSE = zeros(1,iT);
ArmaCoef = zeros(1,iT);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Initialization ----------------------
mu=0; 
sigm=0.05; 
MU=mu*ones(1,N); 
Cxx= (sigm^2)*diag(ones(N,1));
R = chol(Cxx); %Cholesky of Covariance Matrix
StateVec = repmat(MU,iT,1) + randn(iT,N)*R; %Generating a Multivariate Gaussian Distribution with given mean vector and Covariance Matrix Cxx

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Estimation of the distribution of arma coefficient ----------------------
for j = 1:iT
    obsVec(1,j) = EGFR_Signal(1,j);
    ObsMat(1,j) = obsVec(1,j)*obsVec(1,j)'; 
    coef = 0.61; 
    gamFunc(1,j) = (coef.^-1).*ObsMat(1,j);
    priorArmaCoef = exp(sum(diag(gamFunc(1,j))));
    ArmaCoef(j) = exp(sum(diag(gamFunc(1,j)))); 
end
mean_arma_coef = mean(ArmaCoef);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Advancing Particles ----------------------
for t = 1:iT
%% Calculate the prior density ----------------------   
    for i=1:N
        statexPart(:,i) = priorArmaCoef.*state_x(:,i);  
    end   

%% Calculate importance weights ----------------------   
    obsVec(1,t) = EGFR_Signal(1,t); 
    variance(t,:) = exp(statexPart(t,:));
    invCovMat(t,:) = 1./abs(variance(t,:)); 
    logWeight_1 = 0.0015; 
    logWeight_2 = -0.5*iT*log(2*pi);

    logWeight_3= -0.5*log(abs(variance(t,:))); 
    logWeight_4 = -0.5*obsVec(1,t)'*invCovMat(t,:)*obsVec(1,t);
    logWeight_5 = sum((obsVec(t) - statexPart(t))'*ArmaCoef(1,t)*(obsVec(t) - statexPart(t)));

    logWeight_T = logWeight_1 + logWeight_2 + logWeight_3 + logWeight_4 + logWeight_5;
    logWeight = logWeight0 + logWeight_T;
    
%% Stabilize the importance weight ----------------------
    dMaxWeight = max(logWeight);
	vWeights = exp(logWeight - dMaxWeight);
    vWeights = vWeights./sum(vWeights);
    sisrPartState(t) = sum(vWeights.*statexPart(t,:));

%% Compute the ESS ----------------------
    n_thr = 0.5*N;
    vESS(t) = 1/sum(vWeights.^2);
    
%% Resample and reset the weights ----------------------
    if(vESS(t)<=n_thr)
        vIndex = systematicR_matlab(1:N,vWeights');
        statexPart = statexPart(:,vIndex);
        vWeights = ones(1,N)*(1/N);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% statement
for l = 1:iT
    if sisrPartState(l) < 0
        EGFR_Model(l) = 0; 
    else
       EGFR_Model(l) = sisrPartState(l);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Mean Square error
for k = 1:iT
    MSE(1,k) = sqrt(mean(abs(EGFR_Model(1,k) - EGFR_Signal(1,k)).^2));
    err(1,k) = immse(EGFR_Model(1,k),EGFR_Signal(1,k));
end

%%
% sisrPartState = sisrPartState + .1.*ones(1,100);
figure
plot(EGFR_Signal,'b')
hold on
plot(EGFR_Model,'r')
xlabel('Time serie mutation')
ylabel('Magnitude')


