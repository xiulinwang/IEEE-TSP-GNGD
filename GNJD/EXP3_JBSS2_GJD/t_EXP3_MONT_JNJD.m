%% Notes
%  1. This program reproduces the results of Monte-Carlo runs for JNJD
%     of Experiment 3 in the reference below. 
%  2. The results are stored in the folder named RESULTS in folder
%     EXP_3.
%  3. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 3 and thus DO NOT provide codes for their results.
%  4. The authors DO NOT gurantee reproduction of exactly identical results
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
%% Parameters
mont = 10;
N = 5; R = 3; T = 2000;
mix = 1; src = 2; nse = 1;
cond = []; rou = 0; snr = 2;
tgt_opt = 0;
com_opt = 0;
switch tgt_opt
    case 0
        alfa = 0.5; beta = T/20; tau = [];
    case 1
        alfa = [];  beta = [];   tau = [1:T/100:T/5];
    case 2
        alfa = [];  beta = [];   tau = [];
end

%%
% for R = 4:4:20
for snr = 6:10
    J_ISI  = zeros(mont,1);
    ISI    = zeros(mont,R);
    PI     = zeros(mont,R);    
    time_  = zeros(mont,1);
    W      = zeros(N*R,N,mont);        
    for numm = 1:mont
        clc;
        fprintf('--------------JNJD---------------\n');
        if rem(numm,10) == 1
            fprintf('\n The %d st Monte-Carlo run: \n',numm);
        elseif rem(numm,10) == 2
            fprintf('\n The %d nd Monte-Carlo run: \n',numm);
        elseif rem(numm,10) == 3
            fprintf('\n The %d rd Monte-Carlo run: \n',numm);
        else
            fprintf('\n The %d th Monte-Carlo run: \n',numm);
        end
        cd ..
        cd DATA
        [X, AA, X_W, P, St, info] =  f_EXP3_DATA(N, mix, cond, src, T, R, nse, rou, snr);
        cd ..
        cd JBSS_FUNCS
       %% Separate
        tic
        [Y1, W(:,:,numm), sub_off, sum_off, sub_diag, sum_diag] = f_JBSS_SOS_JNJD(X,tgt_opt,alfa,beta,tau);
        time_(numm) = toc;
        [ISI(numm,:),J_ISI(numm)] = f_ISI(AA,W(:,:,numm));        
        temp  = sum_off./sum_diag;
        off{numm} = temp;
        fprintf('running time£º%d s\n',time_(numm));
        fprintf('----------------------------------\n');      
    end
    cd ..
    cd EXP3_JBSS2_GJD
    cd RESULTS
    file_name = ['EXP3_SOS_JNJD_SNR=' num2str(snr) 'dB_N=' num2str(N) '_R=' num2str(R) '_T=' num2str(T) '_MIX=' num2str(mix) '_SRC=' num2str(src) '_NSE=' num2str(nse) '_ROU=' num2str(rou) '_TGT=' num2str(tgt_opt) '.mat'];
    save (file_name,'snr','N','R','T','mont','J_ISI','ISI','time_','W','off','rou','src','mix','nse','tgt_opt','alfa','beta','tau','info');  
    cd ..
end
% end
    