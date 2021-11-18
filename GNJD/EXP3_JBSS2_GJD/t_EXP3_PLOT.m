%% Notes
%  1. This program plots the results from Monte-Carlo runs of GNJD, JNJD,
%     and LUCJD of Experiment 3 in the reference below. 
%  2. The results are stored in the folder named RESULTS in folder
%     EXP_3.
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
clear all;
clc;
cd RESULTS;
alg_list = {'GNJD','JNJD','NJD','NOJoB','GOJD','MCCA'};
token = {'-b','-c*',':md','--g','.-k','-r.'};
no_alg = length(alg_list);
JNT_ISI = cell(1,length(alg_list));
AVE_ISI = JNT_ISI;
AVE_TIME = AVE_ISI;
lower_bound = 6; 
MEP_step = 1;
upper_bound = 10;
show = 1;    % options for plotting: 1 - plot; 2-semilogy
snr = 2; N = 5; R = 3; T = 2000; rou = 0; src = 2; nse = 1; mix = 1; 
no_meth = 3; % number of methods in alg_list that are performed, 
             % it is 3 indicating the first 3 methods in the list are plotted 
tgt = 0;
%% 
N_str = num2str(N); R_str = num2str(R); T_str = num2str(T); rou_str = num2str(rou); src_str = num2str(src); mix_str = num2str(mix); nse_str = num2str(nse); tgt_str = num2str(tgt);
SNR = num2str(snr);
for m = 1:no_meth
    temp_J_ISI  = [];
    temp_A_ISI  = [];
    temp_TIME = [];
    for snr = lower_bound:MEP_step:upper_bound
%     for R = lower_bound:MEP_step:upper_bound
        R_str = num2str(R);
        SNR = num2str(snr);        
        Path = pwd; 
        file = ['\EXP3_SOS_' alg_list{m} '_SNR=' SNR 'dB_N=' N_str '_R=' R_str '_T=' T_str '_MIX=' mix_str '_SRC=' src_str '_NSE=' nse_str '_ROU=' rou_str '_TGT=' tgt_str '.mat'];
        Path = [Path, file];        
        load(Path);
        temp_J_ISI = [temp_J_ISI,J_ISI];
        temp_A_ISI = [temp_A_ISI,mean(ISI,2)];
        temp_TIME  = [temp_TIME, time_];
    end
    JNT_ISI{m}  = mean(temp_J_ISI,1);
    AVE_ISI{m}  = mean(temp_A_ISI,1);
    AVE_TIME{m} = mean(temp_TIME ,1);    
end
cd ..
%% Displaying the results
MEP = lower_bound:MEP_step:upper_bound;
figure (1)
switch show
    case {1} 
        figure (1)
        for m = 1:no_meth
            plot(MEP,JNT_ISI{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('J-ISI');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
        figure (2)
        for m = 1:no_meth
            plot(MEP,AVE_ISI{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('ISI');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
        figure (3)
        for m = 1:no_meth
            plot(MEP,AVE_TIME{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('Running time(s)');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
    case {2}
        figure (1)
        for m = 1:no_meth
            semilogy(MEP,JNT_ISI{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('J-ISI');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
        figure (2)
        for m = 1:no_meth
            semilogy(MEP,AVE_ISI{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('ISI');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
        figure (3)
        for m = 1:no_meth
            semilogy(MEP,AVE_TIME{m},num2str(token{m}));            
            hold on;
        end
        xlabel('SNR(dB)');
        ylabel('Running time(s)');
        legend('GNJD','JNJD','NJD','NOJoB','GOJD','MCCA');
end