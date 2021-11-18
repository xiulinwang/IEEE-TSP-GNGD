%% Notes
%  1. This program reproduces the results of GNJD and JNJD of Experiment 1
%     in the reference below. 
%  2. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 1 and thus DO NOT provide codes for their results.
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
%--------------------------------------------------------------------------
clear all;
clc;
%--------------------------------------------------------------------------
N = 5; %matrix dimensionality
K = 20;%number of matrices in each NJD set
R = 3;
Mont = 5;
for iter=1:Mont
    if rem(iter,10) == 1
        fprintf('\n The %d st Monte-Carlo run: \n',iter);
    elseif rem(iter,10) == 2
        fprintf('\n The %d nd Monte-Carlo run: \n',iter);
    elseif rem(iter,10) == 3
        fprintf('\n The %d rd Monte-Carlo run: \n',iter);
    else
        fprintf('\n The %d th Monte-Carlo run: \n',iter);
    end
    
    %% initialize
    cd ..
    cd DATA
    [GNJD_Ten, NOJoB_Ten, GOJD_Ten, JNJD_Ten, P, sv_A, info] = f_EXP1_DATA(N,R,K);
    cd ..
    cd GJD_FUNCS
    %% GNJD
    [GNJD_Y, GNJD_W, GNJD_SUB_DIAG, GNJD_SUM_DIAG] = f_GNJD(GNJD_Ten,1);  
    %% JNJD
    [JNJD_W, JNJD_Y, JNJD_SUB_off, JNJD_SUM_off, JNJD_SUB_diag, JNJD_SUM_diag] = f_JNJD(JNJD_Ten, R);
    %% NOJoB
%     [NOJoB_W, NOJoB_Y, NOJoB_SUB_off, NOJoB_SUM_off, NOJoB_SUB_diag, NOJoB_SUM_diag] = f_NOJoB(NOJoB_Ten);
    %% GOJD
%     [GOJD_W, GOJD_PI, GOJD_SUB_off, GOJD_SUM_off, GOJD_SUB_diag, GOJD_SUM_diag] = f_GOJD(GOJD_Ten,sv_A,P);   
    cd ..
    cd EXP1_EXACT_GJD
    semilogy(1:length(GNJD_SUM_DIAG(1,:)),GNJD_SUM_DIAG(1,:)./GNJD_SUM_DIAG(2,:),'-b','linewidth',1.5,'MarkerSize',10);
    hold on;
    semilogy(1:length(JNJD_SUM_off),JNJD_SUM_off./JNJD_SUM_diag,'-.k','linewidth',1.5,'MarkerSize',5);
    hold on;
%     semilogy(1:length(NOJoB_SUM_off),NOJoB_SUM_off./NOJoB_SUM_diag,'--','linewidth',1.5,'color',[0 0.5 0],'MarkerSize',10);
%     hold on; 
%     semilogy(1:length(GOJD_SUM_off),GOJD_SUM_off./GOJD_SUM_diag,'-r.','linewidth',1.5,'MarkerSize',20);
%     hold on;   
    if (iter == 1)
        xlabel('Number of sweeps');
        ylabel('Overall off-norm');
        legend('GNJD','JNJD');
%         legend('GNJD','JNJD','NOJoB','GOJD');
    end
end
%% PI