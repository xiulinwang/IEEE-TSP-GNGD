function [GNJD_Ten, NOJoB_Ten, GOJD_Ten, JNJD_Ten, P, sv_A, info] = f_EXP1_DATA(N,R,K)
%--------------------------------------------------------------------------
% This function synthesize the noiseless data for Experiment 1 in the reference
% below.  
%--------------------------------------------------------------------------
%**************************************************************************
%% inputs:
%     R is the number of datasets
%     N is the matrix dimensionality
%     K is the number of target matrices in each NJD problem. 
%% outputs:
%     The target tensors for competitors
%     P: Pre-whitening matrix for GOJD
%     sv_A: the true loading matrices
%     info: experiment setting
%% References:
%  1. X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%     Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%  2. http://202.118.75.4/gong/GNJD.html
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
%**************************************************************************
%--------------------------------------------------------------------------
%******************************************************************
GNJD_Ten = zeros(N,N*K,R,R);
NOJoB_Ten = GNJD_Ten;
JNJD_Ten = zeros(N*(R-1),N*K);
GOJD_Ten = zeros(N,N*K,R,R);
sv_A = randn(N*R,N)+1i*randn(N*R,N);
%% obtain the target matrices
% calculate the unitary matrix for prewhitening in GOJD
P = zeros(N*R,N);
Q = zeros(N*R,N);

for r=1:R
    [a,b,c]=svd(sv_A(N*(r-1)+1:N*r,:));
    P(N*(r-1)+1:N*r,:)=a*b;   %
    Q(N*(r-1)+1:N*r,:)=(P(N*(r-1)+1:N*r,:))\sv_A(N*(r-1)+1:N*r,:); %
end
%% GNJD target matrices
for r1 = 1:R
    for r2 = r1:R
        for m=1:K
            if r1 ~= r2
                D = diag(randn(N,1)) + i*diag(randn(N,1));  
                GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2) = sv_A((r1-1)*N+1:r1*N,:)*D*sv_A((r2-1)*N+1:r2*N,:)';                    
                GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2)/norm(GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2),'fro'); 
                GNJD_Ten(:,(m-1)*N+1:m*N,r2,r1) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2)';
                NOJoB_Ten(:,(m-1)*N+1:m*N,r1,r2) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2); 
                NOJoB_Ten(:,(m-1)*N+1:m*N,r2,r1) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2)'; 
            else %% if r1 == r2
                if R > 5
                    D = diag(randn(N,1)) + i*diag(randn(N,1));  
                    GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2) = sv_A((r1-1)*N+1:r1*N,:)*D*sv_A((r2-1)*N+1:r2*N,:)'; 
                    GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2)/norm(GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2),'fro');                     
                    NOJoB_Ten(:,(m-1)*N+1:m*N,r1,r2) = GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2); 
                end                
            end
        end
    end
end

%%  JNJD target matrices
for r=1:R-1
    for k=1:K
        JNJD_Ten(N*(r-1)+1:N*r,N*(k-1)+1:N*k)=GNJD_Ten(:,N*(k-1)+1:N*k,r,r+1);
    end
end
%% GOJD target matrices
for r1=1:R
    for r2=1:R
        for m=1:K
            GOJD_Ten(:,(m-1)*N+1:m*N,r1,r2)=(P(N*(r1-1)+1:N*r1,:))\GNJD_Ten(:,(m-1)*N+1:m*N,r1,r2)/(P(N*(r2-1)+1:N*r2,:)');            
        end
    end
end
%% Information about the generated data
% N:   Matrix dimensionality
% R:   number of datasets
% K:   number of target matrices in each NJD set
% snr: signal-to-noise ratio

info.matrix_dim = N; 
info.no_of_datasets = R;
info.no_of_matrices = K;
end