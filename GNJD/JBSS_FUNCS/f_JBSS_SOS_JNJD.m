function [Y1, W, sub_off, sum_off, sub_diag, sum_diag] = f_JBSS_SOS_JNJD(X,tgt_opt,alfa,beta,tau)
%% This function do Joint Blind Source Separation 
% with Joint-Nonorthogonal-Joint diagonalization(JNJD) and second-order statistics 
%% The inputs:
% X: The observations, N x T x K, where N is the number of channels, 
%                                       T is the number of temporal samples
%                                       K is the number of data sets

% tgt_opt: options for constructing target matrices
%          0 - time-varying covariance matrices, in this case parameters
%              alfa and beta are valid:
%              alfa: the overlapping rate; beta: block length
%          1 - time-shifted covariance matrices, in this case parameter tau
%              is valid: tau: a vector containing time shifts     
%          2 - covariance matrix
%% The outputs:
% Y: N x N*n_blk x K x K, The unloaded approximatedly diagonalized matrices arranged as:
%     Y(1:N,(r-1)*N+1:r*N,k1,k2) = A_r1*Ck_r1_r2*A_r2';
% W: The demixing matrices concatenated as A = [A1;A2;...;AK];
% sub_off/sub_diag/sum_off/sum_diag: Performance indices for all sweeps 
%% References:
%  X.-F. Gong, Q.-H. Lin, K. Wang, "Joint non-orthogonal joint diagonalization based on LU decomposition and Jacobi scheme" 
%  in Proc. CHINASIP’2013, Beijing, China, Jul. 6-10. 2013.
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome. 
%**************************************************************************
%--------------------------------------------------------------------------
%% Calculate the target matrices
switch tgt_opt
    case 0
        overlap = alfa;    % Overlapping rate
        blk_snap = beta;   % Number of samples in each block
        n_blk = round((T-blk_snap)/((1-overlap)*blk_snap)) + 1;% Number of blocks
        ibloc = round((0:n_blk-1)*(T-blk_snap)/(n_blk-1)); % Indices of each block
        tgt_ten = zeros(N*(R-1),N*n_blk); % No compression
        for ii = 1:R-1        
            for kk = 1:n_blk            
                temp_i = X(:,ibloc(kk)+1:ibloc(kk)+blk_snap,ii);                
                temp_j = X(:,ibloc(kk)+1:ibloc(kk)+blk_snap,ii+1);
                tgt_ten(N*(ii-1)+1:N*ii,N*(kk-1)+1:N*kk) = temp_i*temp_j'/blk_snap; 
            end
        end
    case 1
        n_blk = length(tau);
        tgt_ten = zeros(N*(R-1),N*n_blk); % No compression
        for ii = 1:R-1
            for kk = 1:n_blk            
                temp_i = X(:,1:T-tau(kk),ii);
                temp_j = X(:,tau(kk)+1:T,ii+1);                
                tgt_ten(N*(ii-1)+1:N*ii,N*(kk-1)+1:N*kk) = temp_i*temp_j'/(T-tau(kk));
            end
        end
    case 2
        n_blk = 1;
        for ii = 1:R-1
            temp_i = X(:,:,ii);
            temp_j = X(:,:,ii+1);                
            tgt_ten(N*(ii-1)+1:N*ii,:) = temp_i*temp_j'/T;
        end
end

%%
[W, Y1, sub_off, sum_off, sub_diag, sum_diag] = JNJD_H(tgt_ten, R);
W = permute(W,[1,3,2]); W = reshape(W,[N*R,N]);
function [A, Y1, sub_off, sum_off, sub_diag, sum_diag] = JNJD_H(C, R)
%% This function do Joint-Nonorthogonal-Joint diagonalization of R-1 sets
% Generalized-Nonorthogonal-Joint diagonalization and second-order statistics
%**************************************************************************************
% The inputs:
% C: the target matrix
% R: the number of datasets
% sv_A: the loading matrix, N x N x R
% tgt_opt: options for constructing target matrices
%          0 - time-varying covariance matrices, in this case parameters
%              alfa and beta are valid:
%              alfa: the overlapping rate; beta: block length
%          1 - time-shifted covariance matrices, in this case parameter tau
%              is valid: tau: a vector containing time shifts   
% The outputs:
% W: N x N x R, the demixing matrices concatenated as A = [A1;A2;...;AR]; 
%%
[m,r] = size(C);
n = m/(R-1); N = r/n; % n: dimentionality of target matrices; N: number of target matrices in each set of JD problem
                      % R: number of JD problems
ERR = 1e-10;           % Stopping threshold
RBALANCE = 1;         %Normalization of the results every RBALANCE times
ITER = 100;           %Upper bound of iteration number
A = eye(n);           %joint diagonalizer
A = repmat(A,[R,1]);  %复制子矩阵A
err = ERR*n+1;        %error
k = 0;                %Number of sweeps
err_old = 0;
err_new = zeros(R,1);
%LU分解
while err>ERR && k<ITER
      k = k + 1;
      L = eye(n); L = repmat(L,[R,1]);
      U = eye(n); U = repmat(U,[R,1]);   
      %% compute off-norm
      M = C;
      sum_off(k) = 0;
      sum_diag(k) = 0;
      for r1 = 1 : R-1
          Mtx = M((r1-1)*n+1:r1*n,:);
          [c1,c2,c] = offdiag(Mtx);
          sum_off(k) = sum_off(k) + c1;
          sum_diag(k)= sum_diag(k)+ c2; 
          sub_off(k,r1) = c1;
          sub_diag(k,r1) = c2;
      end
      %% U
      for i = 1:n-1 % l in the document
         for j = i+1:n % k in the document
            a = zeros(R,1); 
            cindex = 1:r;
            cindex(i:n:r) = [];
            W = reshape(C,[n,R-1,n,N]); W = permute(W,[3,2,1,4]); W = reshape(W,[n*(R-1),n*N]);
            % U            
            a(1) = -C(i,cindex)*C(j,cindex)'/(C(j,cindex)*C(j,cindex)');
            if abs(a(1))>1, a(1) = sign(a(1))*1;end;
            C(i,:)=a(1)*C(j,:)+C(i,:);
            I=i:n:r; J=j:n:r;
            for ii = 2:R-1
                C_row_ij = C((ii-1)*n+i,cindex) * C((ii-1)*n+j,cindex)';   
                C_row_jj = C((ii-1)*n+j,cindex) * C((ii-1)*n+j,cindex)';
                C_col_ji = W((ii-2)*n+j,cindex) * W((ii-2)*n+i,cindex)'; %@改@  C_col_ij = W((ii-2)*n+i,cindex) * W((ii-2)*n+j,cindex)';
                C_col_jj = W((ii-2)*n+j,cindex) * W((ii-2)*n+j,cindex)';
                a(ii) = -(C_row_ij + C_col_ji)./(C_row_jj + C_col_jj);   %@改@  C_col_ij
                if abs(a(ii))>1, a(ii) = sign(a(ii))*1;end;
                C((ii-1)*n+i,:)=a(ii)*C((ii-1)*n+j,:)+C((ii-1)*n+i,:);   %左乘当前目标矩阵
                C((ii-2)*n+1:(ii-1)*n,I)=a(ii)'*C((ii-2)*n+1:(ii-1)*n,J)+C((ii-2)*n+1:(ii-1)*n,I);  %@改@   右乘前一目标矩阵
                t2 = reshape(C((ii-2)*n+1:(ii-1)*n,:),[n,n,N]);   t2 = permute(t2,[2,1,3]); W2((ii-2)*n+1:(ii-1)*n,:) = reshape(t2,[n,n*N]);
            end
            a(R) = -(W(n*(R-2)+j,cindex)*W(n*(R-2)+i,cindex)')./(W(n*(R-2)+j,cindex)*W(n*(R-2)+j,cindex)');%@改@
            if abs(a(R))>1, a(R) = sign(a(R))*1;end; %?                            
            C((R-2)*n+1:(R-1)*n,I)=a(R)'*C((R-2)*n+1:(R-1)*n,J)+C((R-2)*n+1:(R-1)*n,I); %@改@
            U(i:n:R*n,:) = repmat(a,[1,n]).*U(j:n:R*n,:) + U(i:n:R*n,:);            
         end
      end     
      
      %% L
      for i = 1:n-1
         for j = i+1:n
            a = zeros(R,1); 
            cindex = 1:r;
            cindex(j:n:r) = [];
            W = reshape(C,[n,R-1,n,N]); W = permute(W,[3,2,1,4]); W = reshape(W,[n*(R-1),n*N]);            
            a(1) = -(C(j,cindex)*C(i,cindex)')/(C(i,cindex)*C(i,cindex)');
            if abs(a(1))>1, a(1) = sign(a(1))*1;end;            
            C(j,:)=a(1)*C(i,:)+C(j,:); 
            I=i:n:r; J=j:n:r;
            for ii = 2:R-1                
                C_row_ji = C((ii-1)*n+j,cindex) * C((ii-1)*n+i,cindex)';   
                C_row_ii = C((ii-1)*n+i,cindex) * C((ii-1)*n+i,cindex)';
                C_col_ij = W((ii-2)*n+i,cindex) * W((ii-2)*n+j,cindex)'; %@改@  C_col_ji = W((ii-2)*n+j,cindex) * W((ii-2)*n+i,cindex)';  
                C_col_ii = W((ii-2)*n+i,cindex) * W((ii-2)*n+i,cindex)';
                a(ii) = -(C_row_ji + C_col_ij)./(C_row_ii + C_col_ii);   %@改@  C_col_ji
                if abs(a(ii))>1, a(ii) = sign(a(ii))*1;end;
                C((ii-1)*n+j,:)=a(ii)*C((ii-1)*n+i,:)+C((ii-1)*n+j,:); 
                C((ii-2)*n+1:(ii-1)*n,J)=a(ii)'*C((ii-2)*n+1:(ii-1)*n,I)+C((ii-2)*n+1:(ii-1)*n,J);%@改@   右乘前一目标矩阵
                t2 = reshape(C((ii-2)*n+1:(ii-1)*n,:),[n,n,N]);   t2 = permute(t2,[2,1,3]); W((ii-2)*n+1:(ii-1)*n,:) = reshape(t2,[n,n*N]);
            end
            a(R) = -(W(n*(R-2)+i,cindex)*W(n*(R-2)+j,cindex)')/(W(n*(R-2)+i,cindex)*W(n*(R-2)+i,cindex)');%@改@ 
            if abs(a(R))>1, a(R) = sign(a(R))*1;end; %?            
            C((R-2)*n+1:(R-1)*n,J)=a(R)'*C((R-2)*n+1:(R-1)*n,I)+C((R-2)*n+1:(R-1)*n,J);   %@改@          
            L(j:n:R*n,:) = repmat(a,[1,n]).*L(i:n:R*n,:) + L(j:n:R*n,:);            
         end
      end
      for ii = 1:R           
          A((ii-1)*n+1:ii*n,:) = L((ii-1)*n+1:ii*n,:)*U((ii-1)*n+1:ii*n,:)*A((ii-1)*n+1:ii*n,:); %
          err_new(ii) = norm(L((ii-1)*n+1:ii*n,:)*U((ii-1)*n+1:ii*n,:)-eye(n,n),'fro')^2;
      end
      err = abs(err_old-max(err_new));
      err_old = max(err_new);
      
      % Normalization 
      if rem(k,RBALANCE)==0
          d = zeros(R,1);          
          for ii = 2:R
              if ii == 2
                  d(1) = norm(A(1:n,:),'fro');
                  A(1:n,:) = A(1:n,:)/d(1);                  
              end
              d(ii) = norm(A((ii-1)*n+1:ii*n,:),'fro');
              A((ii-1)*n+1:ii*n,:) = A((ii-1)*n+1:ii*n,:)/d(ii);              
              for t = 1:N                  
                  C((ii-2)*n+1:(ii-1)*n,(t-1)*n+1:t*n) = d(ii-1)*C((ii-2)*n+1:(ii-1)*n,(t-1)*n+1:t*n)*d(ii);                   
              end
          end               
      end      
end
Y1 = C;
function [off_d,diag_d,c] = offdiag(M)
% This function calculates the off_norm of a set of matrices
% Input:  M = [M1,M2,...,MK]
% Output: off_d: the absolute off_norm
%            _d: the absolute diag_norm
%             c: relative off_norm w.r.t. diag_norm 
[I,J] = size(M);
K = J/I; 
M = reshape(M,[I,I,K]);
c = 0; c1 = 0;
for ii = 1:K
    mtx = M(:,:,ii);
    mtx1 = abs(mtx - diag(diag(mtx))).^2;
    c = c + sum(mtx1(:));%非对角元模平方之和
    c1 = c1 + norm(diag(diag(mtx)))^2;%对角元模平方之和    
end
off_d = c; diag_d = c1; c = c/c1;