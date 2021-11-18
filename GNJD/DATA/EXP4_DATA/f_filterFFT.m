%%% This code was written by Dr. Fengyu Cong in April 2012
%%% Department of Mathematical Information Technology, University of Jyväskyl?, Finland
%%% Address: PL35(Agora), 40014, Jyväskyl?, Finland
%%% Tel.: +358-40-8053255 (Mobile)
%%% E-mails: fengyu.cong@jyu.fi, fengyucong@gmail.com
%%% Homepage: http://users.jyu.fi/~fecong

%%% Using this code please cite the following article:
%%% Fengyu Cong, Yixiang Huang, Igor Kalyakin, Hong Li, Tiina Huttunen-Scott, Heikki Lyytinen, Tapani Ristaniemi, 
%%% Frequency Response based Wavelet Decomposition to Extract Children's Mismatch Negativity Elicited by Uninterrupted Sound, 
%%% Journal of Medical and Biological Engineering, DOI: 10.5405/jmbe.908 (2012, article in press).

function y=f_filterFFT(x,M,fs,freqLowCut,freqHighCut)

freqBin = fs/M;%% frequency represented by one frequency bin 

%% parameters for FFT-filter
binLowCut = ceil(freqLowCut/freqBin);%% the frequency bin corresponds to the low frequency 
binHighCut = ceil(freqHighCut/freqBin); %% the frequency bin corresponds to the high frequency 
binLowCut2 = (M/2+1)+(M/2+1-binLowCut);%% the frequency bin corresponds to the low frequency 
binHighCut2 = (M/2+1)+(M/2+1-binHighCut); %% the frequency bin corresponds to the high frequency 
%% FFT-filter 
%% M-point DFT 
%%%%% signal is length - N 
X = fft(x,M); 
Z = zeros(M,size(x,2)); 
Z(binLowCut:binHighCut,:) = X(binLowCut:binHighCut,:); 
Z(binHighCut2:binLowCut2,:) = X(binHighCut2:binLowCut2,:); 
z = ifft(Z,M); 
y = z(1:size(x,1),:); 