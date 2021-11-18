%% Notes
%  1. This program reproduces the results for GNJD and JNJD of Experiment 4
%     in the reference below. 
%  2. The authors DO NOT own the copyrights of other competing algorithms 
%     in Experiment 4 and thus DO NOT provide codes for their results.
%  3. The authors DO NOT gurantee reproduction of exactly identical results
%     to those in the reference below with this program. Results may vary 
%     with different realizations of independent runs, or software/hardware
%     versions/configurations.
%  4. The real ECG dataset is obtained from DAISY database, see references
%  5. The code for bandpass filter by fft is provided by Prof. Cong Fengyu,
%     please see the comment in the .m file.
%% References:
%  1. X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%     Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%  2. X.-F. Gong, Q.-H. Lin, K. Wang, "Joint non-orthogonal joint diagonalization based on LU decomposition and Jacobi scheme" 
%     in Proc. CHINASIP¡¯2013, Beijing, China, Jul. 6-10. 2013.
%  3. http://202.118.75.4/gong/GNJD.html
%  4. D. De Moor (Ed.), DAISY: Database for the identification of systems. Available online at: http://www.esat.kuleuven.ac.be/sista/daisy.
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
%**************************************************************************
clear all
clc
tic
%% Generate the array data
cd ..
cd DATA
cd EXP4_DATA
X = load('foetal_ECG_data.dat');%2500*9
time_step = X(:,1);   
ecg_data = X(:,2:end);
num_sensor = size(ecg_data,2); % number of sensors
%% Remove the mean
for k = 1 : num_sensor
    ecg_data(:,k) = ecg_data(:,k) - mean(ecg_data(:,k));
end
%% Illustrate the mixed signal
figure (1)
for ii = 1:num_sensor
    subplot(num_sensor,1,ii); plot(time_step,ecg_data(:,ii));
end
%% add a filter
% for k=1:num_sensor
%     data2=ecg_data(:,k); 
%     data=data2-ones(size(data2,1),1)*mean(data2);
%     M=length(time_step);
%     fs=250000;
%     freqLowCut=200;freqHighCut=30000;
%     data_fft1 = f_filterFFT(data,M,fs,freqLowCut,freqHighCut);   % 
%     data_fft = data_fft1 + ones(size(data2,1),1)*mean(data2);    % add the mean 
%     ecg_data(:,k)=data_fft;
% end
%% Generate four datasets
set_no = 5; % num of datasets
dataset = zeros(length(time_step),num_sensor-set_no+1,set_no);
for ii = 1:set_no
    dataset(:,:,ii) = ecg_data(:,ii:ii+num_sensor-set_no);
end
dataset = permute(dataset,[2,1,3]);
w_dataset = dataset;
%% whitening the mixtures
P = zeros(num_sensor-set_no+1,num_sensor-set_no+1,set_no);
for k = 1 : set_no
    R = dataset(:,:,k)*dataset(:,:,k)'/length(time_step);
    P(:,:,k) = inv(sqrtm(R));
    w_dataset(:,:,k) = P(:,:,k)*dataset(:,:,k);
end
%% Separate 
cd ..
cd ..
cd JBSS_FUNCS
tgt_opt = 0;
com_opt = 0;
alfa = 0.5;
beta = 200;
tau = [1:20:400];
mcca_opt = 2;
[Y1,GNJD_W, GNJD_sub_off, GNJD_sum_off, GNJD_sub_diag, GNJD_sum_diag] = f_JBSS_SOS_GNJD(dataset,tgt_opt,com_opt,alfa,beta,tau);
[Y1,JNJD_W, JNJD_sub_off, JNJD_sum_off, JNJD_sub_diag, JNJD_sum_diag] = f_JBSS_SOS_JNJD(dataset,tgt_opt,alfa,beta,tau);
% [GOJD_W, GOJD_sub_off, GOJD_sum_off, GOJD_sub_diag, GOJD_sum_diag] = f_JBSS_SOS_GOJD(w_dataset,P,tgt_opt,com_opt,alfa,beta,tau);
% [Y1, NOJoB_W, NOJoB_sub_off, NOJoB_sum_off, NOJoB_sub_diag, NOJoB_sum_diag] = f_JBSS_SOS_NOJoB(dataset,tgt_opt,com_opt,alfa,beta,tau);
% MCCA_W = f_JBSS_SOS_MCCA(w_dataset,mcca_opt);
GNJD_W = reshape(GNJD_W,[size(GNJD_W,2),set_no,size(GNJD_W,2)]); GNJD_W = permute(GNJD_W,[1,3,2]);
JNJD_W = reshape(JNJD_W,[size(JNJD_W,2),set_no,size(JNJD_W,2)]); JNJD_W = permute(JNJD_W,[1,3,2]);
% NOJoB_W = reshape(NOJoB_W,[size(NOJoB_W,2),set_no,size(NOJoB_W,2)]); NOJoB_W = permute(NOJoB_W,[1,3,2]);
% GOJD_W = reshape(GOJD_W,[size(GOJD_W,2),set_no,size(GOJD_W,2)]); GOJD_W = permute(GOJD_W,[1,3,2]);
% MCCA_W = reshape(MCCA_W,[size(MCCA_W,2),set_no,size(MCCA_W,2)]); MCCA_W = permute(MCCA_W,[1,3,2]);
% recover the whitening 
% for k = 1:set_no
%     GOJD_W(:,:,k) = GOJD_W(:,:,k) * P(:,:,k);
%     MCCA_W(:,:,k) = MCCA_W(:,:,k) * P(:,:,k);
% end
src_est = zeros(size(dataset));
GNJD_src_est = src_est;
JNJD_src_est = src_est;
% GOJD_src_est = src_est;
% NOJOB_src_est = src_est;
% MCCA_src_est = src_est;
for ii = 1:set_no
    GNJD_src_est(:,:,ii) = GNJD_W(:,:,ii)*dataset(:,:,ii);
    JNJD_src_est(:,:,ii) = JNJD_W(:,:,ii)*dataset(:,:,ii);
%     GOJD_src_est(:,:,ii) = GOJD_W(:,:,ii)*dataset(:,:,ii);
%     NOJOB_src_est(:,:,ii) = NOJoB_W(:,:,ii)*dataset(:,:,ii);
%     MCCA_src_est(:,:,ii) = MCCA_W(:,:,ii)*dataset(:,:,ii);
    for jj = 1:size(GNJD_src_est,1)
        GNJD_src_est(jj,:,ii) = GNJD_src_est(jj,:,ii) - mean(GNJD_src_est(jj,:,ii));
        GNJD_src_est(jj,:,ii) = GNJD_src_est(jj,:,ii)/sqrt(var(GNJD_src_est(jj,:,ii)));
        JNJD_src_est(jj,:,ii) = JNJD_src_est(jj,:,ii) - mean(JNJD_src_est(jj,:,ii));
        JNJD_src_est(jj,:,ii) = JNJD_src_est(jj,:,ii)/sqrt(var(JNJD_src_est(jj,:,ii)));
%         GOJD_src_est(jj,:,ii) = GOJD_src_est(jj,:,ii) - mean(GOJD_src_est(jj,:,ii));
%         GOJD_src_est(jj,:,ii) = GOJD_src_est(jj,:,ii)/sqrt(var(GOJD_src_est(jj,:,ii)));
%         NOJOB_src_est(jj,:,ii) = NOJOB_src_est(jj,:,ii) - mean(NOJOB_src_est(jj,:,ii));
%         NOJOB_src_est(jj,:,ii) = NOJOB_src_est(jj,:,ii)/sqrt(var(NOJOB_src_est(jj,:,ii)));
%         MCCA_src_est(jj,:,ii) = MCCA_src_est(jj,:,ii) - mean(MCCA_src_est(jj,:,ii));
%         MCCA_src_est(jj,:,ii) = MCCA_src_est(jj,:,ii)/sqrt(var(MCCA_src_est(jj,:,ii)));
    end
end
toc
cd ..
cd EXP4_ECG_GJD
%%
figure (2)
GNJD_ind = [3,1,4,2];
for ii = 1:set_no    
    for jj = 1:num_sensor-set_no+1
        subplot(set_no,num_sensor-set_no+1,(ii-1)*(num_sensor-set_no+1)+jj);
        plot(time_step,GNJD_src_est(GNJD_ind(jj),:,ii));
        axis([0,10,-10,10]);
    end
end
JNJD_ind = [2,1,4,3];
figure (3)
for ii = 1:set_no    
    for jj = 1:num_sensor-set_no+1
        subplot(set_no,num_sensor-set_no+1,(ii-1)*(num_sensor-set_no+1)+jj);
        plot(time_step,JNJD_src_est(JNJD_ind(jj),:,ii));
        axis([0,10,-10,10]);
    end
end
% GOJD_ind = [1,3,4,2];
% figure (4)
% for ii = 1:set_no    
%     for jj = 1:num_sensor-set_no+1
%         subplot(set_no,num_sensor-set_no+1,(ii-1)*(num_sensor-set_no+1)+jj);
%         plot(time_step,GOJD_src_est(GOJD_ind(jj),:,ii));
%         axis([0,10,-10,10]);
%     end
% end
% NOJOB_ind = [2,1,3,4];
% figure (5)
% for ii = 1:set_no    
%     for jj = 1:num_sensor-set_no+1
%         subplot(set_no,num_sensor-set_no+1,(ii-1)*(num_sensor-set_no+1)+jj);
%         plot(time_step,NOJOB_src_est(NOJOB_ind(jj),:,ii));
%         axis([0,10,-10,10]);
%     end
% end
% MCCA_ind = [4,3,2,1];
% figure (6)
% for ii = 1:set_no    
%     for jj = 1:num_sensor-set_no+1
%         subplot(set_no,num_sensor-set_no+1,(ii-1)*(num_sensor-set_no+1)+jj);
%         plot(time_step,MCCA_src_est(MCCA_ind(jj),:,ii));
%         axis([0,10,-10,10]);
%     end
% end