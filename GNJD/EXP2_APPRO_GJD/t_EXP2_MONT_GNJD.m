%% Notes
%  1. This program reproduces the results of Monte-Carlo runs for GNJD
%     of Experiment 2 in the reference below. 
%  2. The results are stored in the folder named RESULTS in folder
%     EXP_2.
%  3. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 2 and thus DO NOT provide codes for their results.
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
N = 5; R = 3; K = 20;
%%
for snr = 16:20
    J_ISI  = zeros(mont,1);
    ISI    = zeros(mont,R);
    time_  = zeros(mont,1);
    W      = zeros(N*R,N,mont);    
    for numm = 1:mont
        clc;
        fprintf('--------------GNJD---------------\n');
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
        [GNJD_Ten, ~, ~, ~, ~, AA,info] = f_EXP2_DATA(N,R,K,snr);
        cd ..
        cd GJD_FUNCS
        %% Separate
        tic
        [Y1, W(:,:,numm), P1, P2] = f_GNJD(GNJD_Ten,1);     
        %  P1/P2: Performance indices for all sweeps including 
%       P1 = [sub_off; sub_diag];
%       P2 = [sum_off; sum_diag];
        time_(numm) = toc;
        [ISI(numm,:),J_ISI(numm)] = f_ISI(AA,W(:,:,numm));        
        temp  = P2(1,:)./P2(2,:);
        off{numm} = temp;
        fprintf('running time£º%d s\n',time_(numm));
        fprintf('----------------------------------\n');        
    end
    cd ..
    cd EXP2_APPRO_GJD
    cd RESULTS
    file_name = ['EXP2_GNJD_SNR=' num2str(snr) 'dB_N=' num2str(N) '_R=' num2str(R) '_K=' num2str(K) '.mat'];
    save (file_name,'snr','N','R','K','mont','J_ISI','ISI','time_','W','off','info');  
    cd ..    
end
    