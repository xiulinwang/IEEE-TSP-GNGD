function [X, AA, X_W, P, St, info] =  f_EXP3_DATA(N, mix_opt, cond, src, T, R, nse, rou, snr)
%--------------------------------------------------------------------------
% This function synthesize the noisy multi-set data for Experiment 3 in the
% reference below.  
%--------------------------------------------------------------------------
%**************************************************************************
%% Parameter book
% INPUTS:
% ----- Mixing matrix parameters -----
% N: number of observations
% mix_opt: options for mixing matrix
%          1 - randomly generated mixing matrices
%          2 - mixing matrices with a given condition number
% cond   : condition number for mixing matrices

% -----------  Source parameters -----------
% src: types of sources
%      1 - BPSK signals
%      2 - BPSK signals with short time non-stationarity
%      3 - Xilin Li's simulated signals (reference 2)
% T:   number of snapshots
% R:   number of datasets

% -----------  Noise parameters -----------
% nse: types of noises
%      1 - colored stationary gaussian noise with covariance rou
%      2 - colored non-stationary gaussian noise with covariance rou
% rou: spatial colorness between adjacent components
% snr: signal-to-noise ratio in terms of dB

% OUTPUTS:
% X   : array observations
% AA  : mixing matrices
% X_W : pre-whitened mixtures
% P   : de-whitening matrix
% St  : source signals
% info: information about the simulation settings
%% References:
%  1. X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%     Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%  2. http://202.118.75.4/gong/GNJD.html
%  3. X.-L. Li, T. Adal?, M. Anderson, "Joint blind source separation by generalized joint diagonalization of cumulant matrices" 
%     Signal Process., vol. 91, no. 10, pp. 2314-2322, Oct. 2011.
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.

%**************************************************************************
%--------------------------------------------------------------------------
%******************************************************************


%% Constants
Pn=0.01;                                                                                
Ps=Pn*10^(snr/10);
J = sqrt(-1);
%% Select the type of sources
St = zeros(N,T,R);
switch src
    case 1
        %% BPSK signals
        symb = [1+J,-1-J];
        St = randsrc(N,T*R,symb);
        St = reshape(St,[N,T,R]);
        for n = 1 : N
            temp  = squeeze(St(n,:,:)); % T x R
            B = randn(R,R) + J*randn(R,R);
            temp = temp * B; 
            St(n,:,:) = temp;            
            for k = 1 : R    
                St(n,:,k) = St(n,:,k) - mean(St(n,:,k));
                St(n,:,k) = St(n,:,k)/sqrt(var(St(n,:,k)));
            end
        end
    case 2
        %% BPSK signals with short-time nonstationarity
        blk_snap = T/10;  % Number of samples in each block
        overlap = 0.5;   % Overlapping rate
        n_blk = floor((T - blk_snap * overlap)/(blk_snap*(1-overlap))); % Number of blocks
        frame_signal = zeros(N*n_blk,T,R);
        St = zeros(N,T,R);
        symb = [1+J,-1-J];
        for ii = 1:n_blk
            prob = randn(N,1,R);
            prob = repmat(prob,[1,blk_snap]);
            temp_cur = randsrc(N,blk_snap*R,symb);
            temp_cur = reshape(temp_cur,[N,blk_snap,R]);
            temp_cur = temp_cur.*prob/sqrt(2);
            frame_signal(ii:n_blk:end,(blk_snap*(1-overlap)*(ii-1)+1):(blk_snap*(1-overlap)*(ii-1)+blk_snap),:) = temp_cur;            
        end 
        ind = find(frame_signal);
        temp_1 = frame_signal;
        temp_1(ind) = 1;
        temp_1 = sum(temp_1,1)/N;
        for jj = 1:N
            St(jj,:,:) = sum(frame_signal((jj-1)*n_blk+1:jj*n_blk,:,:),1)./temp_1;            
        end        
        for n = 1 : N
            temp  = squeeze(St(n,:,:)); % T x R
            B = randn(R,R) + J*randn(R,R);
            temp = temp * B; 
            St(n,:,:) = temp;            
            for k = 1 : R    
                St(n,:,k) = St(n,:,k) - mean(St(n,:,k));
                St(n,:,k) = St(n,:,k)/sqrt(var(St(n,:,k)));
            end
        end
    case 3
        %% Li xilin        
        St = zeros(N,T,R);
        for n = 1 : N
            temp1 = randn(R,T)+J*randn(R,T);  
            temp  = zeros(R,T);
            B = randn(R,R,3)+J*randn(R,R,3);
            for p = 0 : 2
                for t = 3 : T
                    temp(:,t) = temp(:,t) + B(:,:,p+1)*temp1(:,t-p); % introduce nonwhiteness and spatial correlation
                end
            end
            for k = 1 : R
                St(n,:,k) = temp(k,:);
                St(n,:,k) = St(n,:,k) - mean(St(n,:,k));
                St(n,:,k) = St(n,:,k)/sqrt(var(St(n,:,k)));
            end
        end
end

%% generate the mixtures
switch mix_opt
    case 1
        AA = randn(N,N,R)+J*randn(N,N,R);
    case 2        
        for ii = 1:R
            leig = cond;
            seig = 1;
            eigv = seig + (leig-seig).*rand(N,1);
            eigv(1) = leig;
            eigv(end) = seig;
            eigv = diag(eigv);
            U = randn(N) + J*randn(N);
            V = randn(N) + J*randn(N);
            U = orth(U); V = orth(V);
            AA(:,:,ii) = U*eigv*V';            
        end
end
X = zeros(N,T,R);
for k = 1 : R
    X(:,:,k) = AA(:,:,k)*St(:,:,k);
end

% Remove mean
for k = 1 : R
    for n = 1 : N
        X(n,:,k) = X(n,:,k) - mean(X(n,:,k));
    end
end
Signal = X;
%% noise
switch nse    
    case 1 
        %% Colored stationary noise
        N_c = zeros(N,T,R);
        for r = 1:R
            N_c1 = crlt_noise(T,rou,N,1);
            N_c(:,:,r) = N_c1;
        end
        % introduce spatial correlation
        for n = 1 : N
            temp  = squeeze(N_c(n,:,:)); % T x R
            B = randn(R,R) + J*randn(R,R);
            temp = temp * B; 
            N_c(n,:,:) = temp;
            for k = 1 : R        
                N_c(n,:,k) = N_c(n,:,k) - mean(N_c(n,:,k));
                N_c(n,:,k) = N_c(n,:,k)/sqrt(var(N_c(n,:,k))); 
            end
        end
    case 2
        %% Colored non-stationary noise
        blk_snap = T/10;  % Number of samples in each block
        overlap = 0;      % Overlapping rate
        n_blk = floor((T - blk_snap * overlap)/(blk_snap*(1-overlap))); % Number of blocks
        frame_signal = zeros(N*n_blk,T,R);
        N_c = zeros(N,T,R);
        N_c1 = zeros(N,blk_snap,R);
        for ii = 1:n_blk
            prob = rand(N,1,R);
            prob = repmat(prob,[1,blk_snap]); % amplitude modulation coefficients
            for r = 1:R
                N_c1(:,:,r) = crlt_noise(blk_snap,rou,N,1);                
            end
            N_c1 = N_c1.*prob/sqrt(2);
            frame_signal(ii:n_blk:end,(blk_snap*(1-overlap)*(ii-1)+1):(blk_snap*(1-overlap)*(ii-1)+blk_snap),:) = N_c1;            
        end 
        ind = find(frame_signal);
        temp_1 = frame_signal;
        temp_1(ind) = 1;
        temp_1 = sum(temp_1,1)/N;
        for jj = 1:N
            N_c(jj,:,:) = sum(frame_signal((jj-1)*n_blk+1:jj*n_blk,:,:),1)./temp_1;            
        end        
        % introduce spatial correlation
        for n = 1 : N
            temp  = squeeze(N_c(n,:,:)); % T x R
            B = randn(R,R) + J*randn(R,R);
            temp = temp * B; 
            N_c(n,:,:) = temp;
            for k = 1 : R        
                N_c(n,:,k) = N_c(n,:,k) - mean(N_c(n,:,k));
                N_c(n,:,k) = N_c(n,:,k)/sqrt(var(N_c(n,:,k))); 
            end
        end
end        
%% Mixture
X = Ps*Signal + Pn*N_c;
%% whitening the mixtures
P = zeros(N,N,R);
X_W = zeros(size(X));
for k = 1 : R
    R1 = X(:,:,k)*X(:,:,k)'/T;
    P(:,:,k) = inv(sqrtm(R1));
    X_W(:,:,k) = P(:,:,k)*X(:,:,k);
end
P = permute(P,[1,3,2]); P = reshape(P,[N*R,N]);
AA = permute(AA,[1,3,2]); AA = reshape(AA,[N*R,N]);
%% Information about the generated data
% src: types of sources
%      1 - BPSK signals
%      2 - BPSK signals with short time non-stationarity
%      3 - Xilin Li's simulated signals
% T:   number of snapshots
% R:   number of datasets

% -----------  Noise parameters -----------
% nse: types of noises
%      1 - colored stationary gaussian noise with covariance rou
%      2 - colored non-stationary gaussian noise with covariance rou
mix_list = {['random'],['with preset condition number']};
src_list = {['BPSK'],['Nonstationary BPSK'],['Xilin Li']};
nse_list = {['colored stationary gaussian noise'],['colored non-stationary gaussian noise']};
info.mix_type = mix_list{mix_opt};
info.mix_dim = N; 
info.mix_cond = cond;
info.dataset = R;
info.source_type = src_list{src};
info.source_snapshots = T;
info.noise_type = nse_list{nse};
info.noise_colorness = rou;
info.noise_snr = snr;
function crlt_noise = crlt_noise(snapshot,rou,space_spot,Pn)
% this function is used to generate gaussian spatial colored noises.
% the generated spatial colored noises are spatially correlated with the correlation coefficient equaling rou
% space_spot is the number of colored noise
%randn('state',sum(100*clock));
Nt=randn(1,2*snapshot);
%randn('state',sum(100*clock));
Nt=Nt+1i*randn(1,2*snapshot); 
crlt_Nt = Nt;
for i = 2:space_spot
    temp = filter([rou,sqrt(1-rou^2)], 1 , crlt_Nt(i-1,:));
    temp1=randn(1,snapshot);
    temp1=temp1+1i*temp1;
    temp(2:2:end) = temp1;
    crlt_Nt = [crlt_Nt;temp];
end
crlt_noise = crlt_Nt(:,1:2:end);
crlt_noise = sqrt(Pn) * crlt_noise;