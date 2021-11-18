%% Notes
%  1. This program reproduces the results of 1 realization for GNJD, JNJD,
%     and LUCJD of Experiment 3 in the reference below. 
%  2. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 3 and thus DO NOT provide codes for their results.
%  3. The authors DO NOT gurantee reproduction of exactly identical results
%     to those in the reference below with this program. Results may vary 
%     with different realizations of independent runs, or software/hardware
%     versions/configurations
%% References:
%  1. X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%     Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%  2. X.-F. Gong, Q.-H. Lin, K. Wang, "Joint non-orthogonal joint diagonalization based on LU decomposition and Jacobi scheme" 
%     in Proc. CHINASIP¡¯2013, Beijing, China, Jul. 6-10. 2013.
%  3. http://202.118.75.4/gong/GNJD.html
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
% *************************************************************************
clear all
clc
N = 5; R = 5; T = 2000;
mix = 1; src = 2; nse = 1;
cond = []; rou = 0.8; snr = 3;
cd ..
cd DATA
[X, AA, X_W, P, St, info] =  f_EXP3_DATA(N, mix, cond, src, T, R, nse, rou, snr);
%% display information about datasets
info
%% 
tgt_opt = 0; 
com_opt = 0;
mcca_opt = 2;
switch tgt_opt
    case 0
        alfa = 0.5; beta = T/20; tau = [];
    case 1
        alfa = [];  beta = [];   tau = [1:20:400];
    case 2
        alfa = [];  beta = [];   tau = [];
end
%% Separate
cd ..
cd JBSS_FUNCS
%% GNJD
tic
[Y1,GNJD_W, GNJD_sub_off, GNJD_sum_off, GNJD_sub_diag, GNJD_sum_diag] = f_JBSS_SOS_GNJD(X,tgt_opt,com_opt,alfa,beta,tau);
time_GNJD = toc;
[GNJD_ISI,GNJD_J_ISI] = f_ISI(AA,GNJD_W);
%% JNJD
tic
[Y1,JNJD_W, JNJD_sub_off, JNJD_sum_off, JNJD_sub_diag, JNJD_sum_diag] = f_JBSS_SOS_JNJD(X,tgt_opt,alfa,beta,tau);
time_JNJD = toc;
[JNJD_ISI,JNJD_J_ISI] = f_ISI(AA,JNJD_W);
%% LUCJD
tic
NJD_W = f_JBSS_SOS_LUCJD(X,tgt_opt,alfa,beta,tau);
time_NJD = toc;
[NJD_ISI,NJD_J_ISI] = f_ISI(AA,NJD_W);
cd ..
cd EXP3_JBSS2_GJD
disp(['GNJD_J_ISI: ', num2str(GNJD_J_ISI)]);
disp(['GNJD_ISI  : ', num2str(GNJD_ISI')]);
disp(['JNJD_J_ISI: ', num2str(JNJD_J_ISI)]);
disp(['JNJD_ISI  : ', num2str(JNJD_ISI')]);
disp(['NJD_J_ISI : ', num2str(NJD_J_ISI)]);
disp(['NJD_ISI   : ', num2str(NJD_ISI')]);
%% GOJD
% tic
% [GOJD_W, sub_off, sum_off, sub_diag, sum_diag] = f_JBSS_SOS_GOJD(X_W,P,tgt_opt,com_opt,alfa,beta,tau);
% for ii = 1:R
%     GOJD_W((ii-1)*N+1:ii*N,:) = GOJD_W((ii-1)*N+1:ii*N,:)*P((ii-1)*N+1:ii*N,:);
% end
% time_GOJD = toc;
% [GOJD_ISI,GOJD_J_ISI] = f_ISI(AA,GOJD_W);
%% MCCA
% tic
% MCCA_W = f_JBSS_SOS_MCCA(X_W,mcca_opt);
% for ii = 1:R
%     MCCA_W((ii-1)*N+1:ii*N,:) = MCCA_W((ii-1)*N+1:ii*N,:)*P((ii-1)*N+1:ii*N,:);
% end
% time_MCCA = toc;
% [MCCA_ISI,MCCA_J_ISI] = f_ISI(AA,MCCA_W);
%% NOJoB
% tic
% [Y1,NOJoB_W, NOJoB_sub_off, NOJoB_sum_off, NOJoB_sub_diag, NOJoB_sum_diag] = f_JBSS_SOS_NOJoB(X,tgt_opt,com_opt,alfa,beta,tau);
% time_NOJoB = toc;
% [NOJoB_ISI,NOJoB_J_ISI] = f_ISI(AA,NOJoB_W);