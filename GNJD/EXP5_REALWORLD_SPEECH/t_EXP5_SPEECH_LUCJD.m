%% Notes
%  1. This program reproduces the results for LUCJD of Experiment 5 in reference 3. 
%  2. The authors DO NOT own the copyrights of some other competing algorithms 
%     in Experiment 5 and thus DO NOT provide codes for their results.
%  3. The authors DO NOT gurantee reproduction of exactly identical results
%     to those in the reference below with this program. Results may vary 
%     with different realizations of independent runs, or software/hardware
%     versions/configurations.
%  4. The real speech mixture is obtained from SISEC2010 benchmarks, the speech
%     sources are obtained from Prof. Sawada's website see references 5 and
%     6 below.
%  5. The metrics for the separation, namely SIR, SDR and SAR, are obtained from
%     Prof. E. Vincent's website. See reference 7.
%  6. The softwares for performing windowed short time fourier transform
%     and inverse short time fourier transform are obtained from Prof. E.
%     Vincent's website. See reference 7.
%% References:
%  1. X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%     Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%  2. X.-F. Gong, Q.-H. Lin, K. Wang, "Joint non-orthogonal joint diagonalization based on LU decomposition and Jacobi scheme" 
%     in Proc. CHINASIP'2013, Beijing, China, Jul. 6-10. 2013.
%  3. K. Wang, X. -F. Gong, Q. -H. Lin, "Complex non-orthogonal joint diagonalization based on LU and LQ decompositions" 
%     in Proc. LVA/ICA'2012, Tel Aviv, Israel, Mar. 12-15. 2012.
%  4. http://202.118.75.4/gong/GNJD.html
%  5. http://sisec2010.wiki.irisa.fr/tiki-index.php?page=Determined+and+over-determined+speech+and+music+mixtures+speech+and+music+mixtures.
%  6. http://www.kecl.ntt.co.jp/icl/signal/sawada/demo/bss2to4/index.html
%  7. http://www.irisa.fr/metiss/members/evincent/software
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
% -------------------------------------------------------------------------
clear all;
close all;
clc

if ~exist('setting')                  % settings  0: Room 4; 
    setting = 1;                      %           1: Room 5 
end

if ~exist('pplot')                    %  Plot sources/s and its estimates/y
    pplot = 1;                        %  0: No
end                                   %  1: Yes

disp('Speech separation with LUCJD');
%% Room 4 or Room 5
switch setting
    case 0
       %% Room 4
        fftlen = 2048;                % fft length   
        blklen = 9;                   % block length 
        t_dur = 7;                    % speech duration is set to 7 seconds as the sources only have 7 seconds duration
        cd ..
        cd DATA
        cd EXP5_DATA
        x = 'room4_2sources_mix.wav'; % 
        [x,fs] = audioread(x);        % fs = 16000         
        X = x(1:16000*t_dur,1:2);     % speech mixtures trancated into 7 seconds
        slen = size(X,1); 
        s1 = 's1-1.wav';              % source files
        s2 = 's2-1.wav';  
        [s1,fs] = audioread(s1);      % fs = 8000, duration = 7s, samples = 56000
        [s2,fs] = audioread(s2);
        s = [s1,s2];
        s = s(1:8000*t_dur,:);        
    case 1
       %% Room 5
        fftlen = 2048;                % fft length   
        blklen = 4;                   % block length 
        t_dur = 7;                    % speech duration is set to 7s as the sources only have 7s duration
        cd ..
        cd DATA
        cd EXP5_DATA
        x = 'room5_2sources_mix.wav';
        [x,fs] = audioread(x);        
        X = x(1:16000*t_dur,1:2);
        slen = size(X,1); 
        s1 = 's1-1.wav';              % fs = 8000, duration = 7s, samples = 56000
        s2 = 's2-1.wav';  
        [s1,fs] = audioread(s1);
        [s2,fs] = audioread(s2);
        s = [s1,s2];
        if t_dur <= 7
            s = s(1:8000*t_dur,:);
        end
end;
%% STFT to time-frequency domain
cd ..
cd ..
cd EXP5_REALWORLD_SPEECH
nam = ['fft:' num2str(fftlen)];
disp(nam);  
% % % fft
Xfft = f_stft_multi(X',fftlen);  % nbin x nfram x nchan
[nbin,nfram,nsrc] = size(Xfft);
Xfft = permute(Xfft,[3,2,1]);
% source fft
s3 = resample(s,16000,8000);
Sfft = f_stft_multi(s3',fftlen);   % nbin x nfram x nchan
[nbin,nfram,nsrc] = size(Sfft);
Sfft = permute(Sfft,[3,2,1]);
%% Separation
Yfft = zeros(nsrc,nfram,nbin);
Wfft = zeros(nsrc,nsrc,nbin);
cd ..
cd GJD_FUNCS
tic
for fbin=2:nbin      
    %% LUCJD
    X1 = squeeze(Xfft(:,:,fbin));        
    X_len = size(X1,2);
    if fbin >= 3
        Gf = Wfft(:,:,fbin-1);
        Gf = diag(sqrt(sum(abs(Gf).^2,2)/2))\Gf;        
        X1 = Gf * X1;        
    end
    M = nsrc;
    overlap = (blklen-1)/blklen;%overlapping coefficient
    n_blk = round((size(X1,2)-blklen)/((1-overlap)*blklen)) + 1;%
    ibloc = round((0:n_blk-1)*(X_len-blklen)/(n_blk-1));%                 
    ccov_ten = zeros(nsrc,nsrc*n_blk);
    for ii = 1:n_blk
        temp_1 = X1(:,ibloc(ii)+1:ibloc(ii)+blklen);             
        ccov_ten(:,(ii-1)*nsrc+1:ii*nsrc) = temp_1*temp_1';
    end
    %%  Compression to accelerate
    ccov_ten = ccov_ten/blklen;
    ccov_ten = reshape(ccov_ten,[nsrc*nsrc,n_blk]);
    [U,S,V] = svd(ccov_ten,'econ');
    ccov_ten = reshape(U*S,[nsrc,nsrc,nsrc^2]);
    ccov_ten = reshape(ccov_ten,[nsrc,nsrc^3]);    
    [Y,G] = f_LUCJD(ccov_ten,1);    
    Wfft(:,:,fbin) = G;
    if fbin >= 3
        Wfft(:,:,fbin) = Wfft(:,:,fbin) * Gf;
    end    
    Yfft(:,:,fbin) = Wfft(:,:,fbin) * Xfft(:,:,fbin);       
end
toc
%% Permutation alignment via sequential amplitude correlation for every 2 adjacent frequency bins  
order_amplitude=repmat([1;2],1,nbin);
for i=2:nbin-1  
    pA1=corrcoef(abs(Yfft(1,:,i)),abs(Yfft(1,:,i+1)));
    pA2=corrcoef(abs(Yfft(2,:,i)),abs(Yfft(2,:,i+1)));
    pB1=corrcoef(abs(Yfft(1,:,i)),abs(Yfft(2,:,i+1)));
    pB2=corrcoef(abs(Yfft(2,:,i)),abs(Yfft(1,:,i+1)));
    A1=pA1(1,2);A2=pA2(1,2);B1=pB1(1,2);B2=pB2(1,2);
    if (A1<B1)&&(A2<B2)
        Yfft(:,:,i+1) = Yfft([2,1],:,i+1);
        Wfft(:,:,i+1) = Wfft([2,1],:,i+1);
    end
end
%% Minimal distortion principle for scaling de-ambiguity
Yfft_WMDP=zeros(nsrc,nfram,nbin);
W_scale=zeros(nsrc,nsrc,nbin);
for f=2:nbin  
     W_scale(:,:,f)=diag(diag(inv(Wfft(:,:,f))))*Wfft(:,:,f);   %scale correction with MDP
     Yfft_WMDP(:,:,f)=W_scale(:,:,f)*Xfft(:,:,f);            
end

%% Back transform the estimates to the temporal domain
disp( 'Back transform the estimates to the temporal domain');
Yfft=Yfft_WMDP; % with scale correction
S=zeros(nbin,nfram,nsrc);
for f=1:nbin,
    Sf=zeros(nsrc,nfram);
    Sf=Yfft(:,:,f);
    S(f,:,:)=reshape(Sf.',1,nfram,nsrc);
end
cd ..
cd EXP5_REALWORLD_SPEECH
Ye = f_istft_multi(S,slen);

%% smoothing
Y_smooth=zeros(2,slen);
Y_smooth(1,:)=smooth(Ye(1,:),5);  % the smooth length has impact on the final results
Y_smooth(2,:)=smooth(Ye(2,:),5);
if t_dur <= 7
    Y_smooth = Y_smooth(:,1:2:16000*t_dur);
end
if t_dur > 7
    Y_smooth = Y_smooth(:,1:2:16000*7);
end
if pplot == 1
    figure;
    subplot(2,2,1); plot(s1); axis tight; xlabel('samples');ylabel('s1');
    subplot(2,2,2); plot(s2); axis tight; xlabel('samples');ylabel('s2');
    subplot(2,2,3); plot(Y_smooth(2,:)); axis tight; xlabel('samples');ylabel('y1');
    subplot(2,2,4); plot(Y_smooth(1,:)); axis tight; xlabel('samples');ylabel('y2');
end
cd RESULTS
audiowrite('LUCJD_y1.wav',4*Y_smooth(1,:),fs);
audiowrite('LUCJD_y2.wav',4*Y_smooth(2,:),fs);
%% computing SIR
cd ..
[SDR_Y,SIRs_Y,SAR_Y,perm] = f_bss_eval_sources(Y_smooth,s');
cd ..
disp(['SIR_Y=', num2str(mean(SIRs_Y)),'    SDR_Y=',num2str(mean(SDR_Y)), '   SAR_Y=',num2str(mean(SAR_Y))]);
