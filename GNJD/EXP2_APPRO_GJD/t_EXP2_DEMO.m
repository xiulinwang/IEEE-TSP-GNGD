%% Notes
%  1. This program reproduces the results of 1 realization for GNJD and JNJD
%     of Experiment 2 in the reference below. 
%  2. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 2 and thus DO NOT provide codes for their results.
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
%**************************************************************************
clear all
clc
N = 5; R = 5; K = 20;
snr = 10;
cd ..
cd DATA
[GNJD_Ten, NOJoB_Ten, GOJD_Ten, JNJD_Ten, P, sv_A, info] = f_EXP2_DATA(N,R,K,snr);
info
cd ..
cd GJD_FUNCS
%% GNJD
tic
[GNJD_Y, GNJD_W, GNJD_SUB, GNJD_SUM] = f_GNJD(GNJD_Ten,2);  
time_GNJD = toc;
%% JNJD
tic
[JNJD_W, JNJD_Y, JNJD_SUB_off, JNJD_SUM_off, JNJD_SUB_diag, JNJD_SUM_diag] = f_JNJD(JNJD_Ten, R);
time_JNJD = toc;
%% NOJoB
% tic
% [NOJoB_W, NOJoB_Y, NOJoB_SUB_off, NOJoB_SUM_off, NOJoB_SUB_diag, NOJoB_SUM_diag] = f_NOJoB(NOJoB_Ten);
% time_NOJoB = toc;
%% GOJD
% tic
% [GOJD_W, GOJD_SUB_off, GOJD_SUM_off, GOJD_SUB_diag, GOJD_SUM_diag] = f_GOJD(GOJD_Ten);
% for ii = 1:R
%     GOJD_W((ii-1)*N+1:ii*N,:) = GOJD_W((ii-1)*N+1:ii*N,:)/P((ii-1)*N+1:ii*N,:);
% end
% time_GOJD = toc;
%% Evaluation
[GNJD_ISI,GNJD_J_ISI] = f_ISI(sv_A,GNJD_W);
[JNJD_ISI,JNJD_J_ISI] = f_ISI(sv_A,JNJD_W);
% [NOJoB_ISI,NOJoB_J_ISI] = f_ISI(sv_A,NOJoB_W);
% [GOJD_ISI,GOJD_J_ISI] = f_ISI(sv_A,GOJD_W);
cd ..
cd EXP2_APPRO_GJD
disp(['GNJD_J_ISI: ', num2str(GNJD_J_ISI)]);
disp(['GNJD_ISI  : ', num2str(GNJD_ISI')]);
disp(['JNJD_J_ISI: ', num2str(JNJD_J_ISI)]);
disp(['JNJD_ISI  : ', num2str(JNJD_ISI')]);