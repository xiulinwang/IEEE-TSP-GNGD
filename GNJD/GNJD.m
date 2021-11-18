function [Y, A, P1, P2] = GNJD(tgt_ten, opt)
%--------------------------------------------------------------------------
% Generalized Nonorthogonal Joint Diagonalization(GNJD)
% with LU Decomposition and Successive Rotations 
%--------------------------------------------------------------------------
%**************************************************************************
%% inputs:
%  tgt_ten:  The input target matrices are arranged as a fourth-order tensor:
%     Ck_r1_r2 = {Ar1*Sk_r1_r2*Ar2'}; k = 1:K, r1 = 1:R, r2 = r1:R
%     N x NK x R x R
%     R is the number of datasets
%     N is the matrix dimensionality
%     K is the number of target matrices in each NJD problem. 
%  opt: options for parallelization
%       0 - without parallelization
%       1 - columnwise parallelization
%       2 - diagonalwise parallelization
%% outputs:
%  Y:  unloaded target matrices
%  A:  demixing matrix concatenated along column dimension: [A1;A2;...;AR];
%  P1/P2: Performance indices for all sweeps including 
%       P1 = [sub_off; sub_diag];
%       P2 = [sum_off; sum_diag];
%% References:
%  X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%  Pre-print available at arXiv:1312.0712v2. This paper is accepted IEEE Transactions on Signal Processing;
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
%**************************************************************************
%--------------------------------------------------------------------------
clc
if nargin==0   
    help f_GNJD
    disp('Executing test program.')
    DEMO_GNJD;
    return ;
end
%--------------------------------------------------------------------------
% Generalized Nonorthogonal Joint Diagonalization(GNJD)
%--------------------------------------------------------------------------
[N,r,R,R] = size(tgt_ten); % N by NK by R by R
if R < 6  %  remove the intra-statistics of target matrices for small R 
     for r = 1:R
         tgt_ten(:,:,r,r) = zeros(N,r);
     end
end
[Y, A, P1, P2] = f_GNJD(tgt_ten,opt);  

end
function DEMO_GNJD
    % construct the target matrices
    N = 6;             % N: dimentionality of target matrices; 
    K = 20;            % K: number of target matrices in each set of JD problem
    R = 3;             % R: number of multiple datasets
    tgt_ten = zeros(N,N*K,R,R);      % target matrices
    sv_A = randn(N*R,N)+1i*randn(N*R,N); % mixing matrices
    for r1 = 1:R
        for r2 = 1:R
            for m=1:K
                D = diag(randn(N,1));        
                tgt_ten(:,(m-1)*N+1:m*N,r1,r2) = sv_A((r1-1)*N+1:r1*N,:)*D*sv_A((r2-1)*N+1:r2*N,:)';
                tgt_ten(:,(m-1)*N+1:m*N,r1,r2) = tgt_ten(:,(m-1)*N+1:m*N,r1,r2)/norm(tgt_ten(:,(m-1)*N+1:m*N,r1,r2),'fro');  
            end
        end
    end
    opt = 2;
    % Generalized Nonorthogonal Joint Diagonalization(GNJD)
    if R < 6  %  remove the intra-statistics of target matrices for small R 
         for r = 1:R
             tgt_ten(:,:,r,r) = zeros(N,N*K);
         end
    end
    [Y, A, P1, P2] = f_GNJD(tgt_ten,opt); 
    [~,J_ISI] = f_ISI(sv_A,A);   % evaluate the performances 
    semilogy(P2(1,:)./P2(2,:));% off-norm
    xlabel('Number of Iterations')
    legend('GNJD');
    title('ORON'); % overall off-norm
    disp('J_ISI evaluates both accuracy and permutation');
    J_ISI
end
function [Y, A, P1, P2] = f_GNJD(C, opt)
%--------------------------------------------------------------------------
% Generalized Nonorthogonal Joint Diagonalization(GNJD)
% with LU Decomposition and Successive Rotations 
%--------------------------------------------------------------------------
%**************************************************************************
%% inputs:
%  tgt_ten:  The input target matrices are arranged as a fourth-order tensor:
%     Ck_r1_r2 = {Ar1*Sk_r1_r2*Ar2'}; k = 1:K, r1 = 1:R, r2 = r1:R
%     N x NK x R x R
%     R is the number of datasets
%     N is the matrix dimensionality
%     K is the number of target matrices in each NJD problem. 
%  opt: options for parallelization
%       0 - without parallelization
%       1 - columnwise parallelization
%% outputs:
%  Y:  unloaded target matrices
%  A:  demixing matrix concatenated along column dimension: [A1;A2;...;AR];
%  P1/P2: Performance indices for all sweeps including 
%       P1 = [sub_off; sub_diag];
%       P2 = [sum_off; sum_diag];
%% References:
%  X.-F. Gong, X.-L. Wang, Q.-H. Lin "Generalized Non-orthogonal Joint Diagonalization with LU Decomposition and Successive Rotations"
%  Pre-print available at arXiv:1312.0712v2. This paper is accepted by IEEE Transactions on Signal Processing;
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
%**************************************************************************
%--------------------------------------------------------------------------
[N,r,R,R] = size(C); % N by NK by R by R
K = r/N;             % N: dimentionality of target matrices; 
                     % K: number of target matrices in each set of JD problem
                     % R: number of multiple datasets
ERR = 1e-6;          % Stopping threshold
RBALANCE = 1;        % Normalization of the results every RBALANCE times
ITER = 100;          % Upper bound of iteration number
for i=1:R            % joint diagonalizer
    A(:,:,i)=eye(N);
end
err = ERR*N + 1;     % error
k = 0;               % Number of sweeps
err_old = 0;
err_new = zeros(R,1);
switch opt
    case 0
       %% GNJD in its original sequential version
        while err>ERR && k<ITER
              k = k + 1;
              L = eye(N); L = repmat(L,[R,1]);
              U = eye(N); U = repmat(U,[R,1]);   
              %%
              % compute off-norm
              M = C;
              rr = 1;
              sum_off(k) = 0;
              sum_diag(k) = 0;
              for r1 = 1 : R
                  for r2 = r1 : R
                      Mtx = M(:,:,r1,r2);
                      [c1,c2,c] = offdiag(Mtx);
                      sum_off(k) = sum_off(k) + c1;
                      sum_diag(k)= sum_diag(k)+ c2; 
                      sub_off(k,rr) = c1;
                      sub_diag(k,rr) = c2;
                      rr = rr + 1;
                  end
              end
              %% U stage
              for i = 1:N-1 % l in the document
                 for j = i+1:N % k in the document 
                    cindex = 1:r;
                    cindex(i:N:r) = [];    
                    % NK-K by R by R
                    I=i:N:r; J=j:N:r;  
                    W = reshape(C,[N,N,K,R,R]);
                    W = permute(W,[2,1,3,4,5]); 
                    W = reshape(W,[N,N*K,R,R]);                                    
                    ten_row_j = C(j,cindex,:,:); 
                    ten_row_i = C(i,cindex,:,:);
                    ten_col_j = W(j,cindex,:,:); 
                    ten_col_i = W(i,cindex,:,:);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[size(cindex,2),R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[size(cindex,2),R-ii+1]); 
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[size(cindex,2),ii]); 
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[size(cindex,2),ii]);                          
                        v = [v, temp_col_j, temp_row_i];
                        w = [w, temp_col_i, temp_row_j]; 
                        u = [u, temp_col_j, temp_row_j];
                    end
                    f1 = sum(reshape(sum(v.*conj(w),1),[R+1,R]));
                    f2 = sum(reshape(sum(u.*conj(u),1),[R+1,R]));
                    a = -f1./f2;
                    for ii = 1:R    
                        if abs(a(ii))>1, a(ii) = sign(a(ii))*1;end;
                        %% N by NK by R by R
                        C(i,:,ii,ii:R) = a(ii)*C(j,:,ii,ii:R) + C(i,:,ii,ii:R);  % update by left multiply 
                        C(:,I,1:ii,ii) = a(ii)'*C(:,J,1:ii,ii) + C(:,I,1:ii,ii); % update by right multiply 
                    end
                    U(i:N:R*N,:) = repmat(a.',[1,N]).*U(j:N:R*N,:) + U(i:N:R*N,:);            
                 end
              end     

              %% L
              for i = 1:N-1
                 for j = i+1:N
                    cindex = 1:r;
                    cindex(j:N:r) = [];            
                    I=i:N:r; J=j:N:r;            
                    W = reshape(C,[N,N,K,R,R]);
                    W = permute(W,[2,1,3,4,5]); 
                    W = reshape(W,[N,N*K,R,R]);                                    
                    ten_row_j = C(j,cindex,:,:); 
                    ten_row_i = C(i,cindex,:,:);
                    ten_col_j = W(j,cindex,:,:); 
                    ten_col_i = W(i,cindex,:,:);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[size(cindex,2),R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[size(cindex,2),R-ii+1]); 
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[size(cindex,2),ii]); 
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[size(cindex,2),ii]);                          
                        v = [v, temp_col_i, temp_row_j]; 
                        w = [w, temp_col_j, temp_row_i];
                        u = [u, temp_col_i, temp_row_i];
                    end
                    f1 = sum(reshape(sum(v.*conj(w),1),[R+1,R]));
                    f2 = sum(reshape(sum(u.*conj(u),1),[R+1,R]));
                    a = -f1./f2;
                    for ii = 1:R   
                        if abs(a(ii))>1, a(ii) = sign(a(ii))*1;end;
                        %% N by NK by R by R
                        C(j,:,ii,ii:R) = a(ii)*C(i,:,ii,ii:R) + C(j,:,ii,ii:R);  % update by left multiply 
                        C(:,J,1:ii,ii) = a(ii)'*C(:,I,1:ii,ii) + C(:,J,1:ii,ii); % update by right multiply             
                    end            
                    L(j:N:R*N,:) = repmat(a.',[1,N]).*L(i:N:R*N,:) + L(j:N:R*N,:);                        
                 end
              end

              for ii = 1:R           
                  A(:,:,ii) = L((ii-1)*N+1:ii*N,:)*U((ii-1)*N+1:ii*N,:)*A(:,:,ii); %
                  err_new(ii) = norm(L((ii-1)*N+1:ii*N,:)*U((ii-1)*N+1:ii*N,:)-eye(N,N),'fro')^2;
              end
              err = abs(err_old-max(err_new));
              err_old = max(err_new);       
              % Normalization     
              if rem(k,RBALANCE)==0
                  d = zeros(R,1);          
                  for ii = 1:R              
                      d(ii) = norm(A(:,:,ii),'fro');
                      A(:,:,ii) = A(:,:,ii)/d(ii);              
                      for t = 1:K
                          C(:,(t-1)*N+1:t*N,ii,ii:R) = d(ii)*C(:,(t-1)*N+1:t*N,ii,ii:R);
                          C(:,(t-1)*N+1:t*N,1:ii,ii) = d(ii)*C(:,(t-1)*N+1:t*N,1:ii,ii);                  
                      end
                  end               
              end
        end
    case 1
       %% Parallelization in columnwise fashion
        while err>ERR && k<ITER
              k = k + 1;
              L = eye(N); L = repmat(L,[1,R]);
              U = eye(N); U = repmat(U,[1,R]);   
              %%
              % compute off-norm
              M = C;
              rr = 1;
              sum_off(k) = 0;
              sum_diag(k) = 0;
              for r1 = 1 : R
                  for r2 = r1 : R
                      Mtx = M(:,:,r1,r2);
                      [c1,c2,c] = offdiag(Mtx);
                      sum_off(k) = sum_off(k) + c1;
                      sum_diag(k)= sum_diag(k)+ c2; 
                      sub_off(k,rr) = c1;
                      sub_diag(k,rr) = c2;
                      rr = rr + 1;
                  end
              end
              %% U stage
                 for j = 2:N % ;  
                    CC = repmat((ones(N)-eye(N)),[1,K,R,R]);
                    CC = CC.*C;
                    W = reshape(permute(reshape(C,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    WC = reshape(permute(reshape(CC,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    ten_row_j = repmat(C(j,:,:,:),[j-1,1]); 
                    ten_row_i = CC(1:j-1,:,:,:);
                    ten_col_j = repmat(W(j,:,:,:),[j-1,1]); 
                    ten_col_i = WC(1:j-1,:,:,:);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[j-1,N*K,R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[j-1,N*K,R-ii+1]);
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[j-1,N*K,ii]);
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[j-1,N*K,ii]);
                        v = cat(3,v,cat(3,temp_col_j,temp_row_i));
                        w = cat(3,w,cat(3,temp_col_i,temp_row_j)); 
                        u = cat(3,u,cat(3,temp_col_j,temp_row_j));
                    end
                    f1 = reshape(squeeze(sum(v.*conj(w),2)),[size(v,1),R+1,R]);
                    f1 = reshape(sum(f1,2),[size(v,1),R]);
                    [~,q2,q3] = size(u);
                    qq = repmat((ones(N)-eye(N)),[1,q2/N,q3]);
                    f2 = reshape(squeeze(sum(u.*qq(1:j-1,:,:).*conj(u),2)),[size(u,1),R+1,R]);
                    f2 = reshape(sum(f2,2),[size(u,1),R]);
                    a = -f1./f2;
                    for ii = 1:R    
                        if abs(a(:,ii))>1, a(:,ii) = sign(a(:,ii))*1;end;
                        [~,a2,a3,a4] = size(C(j,:,ii,ii:R));
                        aa = repmat(a(:,ii),[1,a2,a3,a4]);
                        %% N by NK by R by R
                        C(1:j-1,:,ii,ii:R) = aa.*repmat(C(j,:,ii,ii:R),[j-1,1]) + C(1:j-1,:,ii,ii:R);   % update by left multiply 
                        I = [];
                        J = [];
                        for r1 = 1:j-1
                            I = [I,r1:N:r];
                            J = [J,j:N:r];
                        end
                        I = sort(I);
                        J = sort(J);
                        [~,~,a3,a4] = size(C(:,J,1:ii,ii));
                        C(:,I,1:ii,ii) = (repmat(a(:,ii)',[N,K,a3,a4])).*C(:,J,1:ii,ii) + C(:,I,1:ii,ii);% update by right multiply 
                    end
                    aa = reshape(permute(repmat(a,[1,1,N]),[1,3,2]),[size(a,1),N*R]);
                    U(1:j-1,:) = U(1:j-1,:) + aa.*repmat(U(j,:),[j-1,1]);         
                 end  

              %% L
              for j = 1:N-1         
                    CC = repmat((ones(N)-eye(N)),[1,K,R,R]);
                    CC = CC.*C;
                    W = reshape(permute(reshape(C,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    WC = reshape(permute(reshape(CC,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    ten_row_j = CC(j+1:N,:,:,:); 
                    ten_row_i = repmat(C(j,:,:,:),[N-j,1]);
                    ten_col_j = WC(j+1:N,:,:,:); 
                    ten_col_i = repmat(W(j,:,:,:),[N-j,1]);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[N-j,N*K,R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[N-j,N*K,R-ii+1]);
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[N-j,N*K,ii]);
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[N-j,N*K,ii]);
                        v = cat(3,v,cat(3,temp_col_i,temp_row_j));
                        w = cat(3,w,cat(3,temp_col_j,temp_row_i)); 
                        u = cat(3,u,cat(3,temp_col_i,temp_row_i));
                    end
                    f1 = reshape(squeeze(sum(v.*conj(w),2)),[size(v,1),R+1,R]);
                    f1 = reshape(sum(f1,2),[size(v,1),R]);
                    [~,q2,q3] = size(u);
                    qq = repmat((ones(N)-eye(N)),[1,q2/N,q3]);
                    f2 = reshape(squeeze(sum(u.*qq(j+1:N,:,:).*conj(u),2)),[size(u,1),R+1,R]);
                    f2 = reshape(sum(f2,2),[size(u,1),R]);
                    a = -f1./f2;
                    for ii = 1:R   
                        if abs(a(:,ii))>1, a(:,ii) = sign(a(:,ii))*1;end;
                        %% N by NK by R by R
                        [~,a2,a3,a4] = size(C(j,:,ii,ii:R));
                        aa = repmat(a(:,ii),[1,a2,a3,a4]);
                        C(j+1:N,:,ii,ii:R) = aa.*repmat(C(j,:,ii,ii:R),[N-j,1]) + C(j+1:N,:,ii,ii:R);  % update by left multiply      
                        I = [];
                        J = [];
                        for r1 = j+1:N
                            I = [I,r1:N:r];
                            J = [J,j:N:r];
                        end
                        I = sort(I);
                        J = sort(J);
                        [~,~,a3,a4] = size(C(:,J,1:ii,ii));
                        C(:,I,1:ii,ii) = (repmat(a(:,ii)',[N,K,a3,a4])).*C(:,J,1:ii,ii) + C(:,I,1:ii,ii);% update by right multiply 
                    end             
                    aa = reshape(permute(repmat(a,[1,1,N]),[1,3,2]),[size(a,1),N*R]);
                    L(j+1:N,:) = L(j+1:N,:) + aa.*repmat(L(j,:),[N-j,1]);
              end

              for ii = 1:R           
                  A(:,:,ii) = L(:,(ii-1)*N+1:ii*N)*U(:,(ii-1)*N+1:ii*N)*A(:,:,ii); %
                  err_new(ii) = norm(L(:,(ii-1)*N+1:ii*N)*U(:,(ii-1)*N+1:ii*N)-eye(N,N),'fro')^2;
              end
              err = abs(err_old-max(err_new));
              err_old = max(err_new);       
              % Normalization     
              if rem(k,RBALANCE)==0
                  d = zeros(R,1);          
                  for ii = 1:R              
                      d(ii) = norm(A(:,:,ii),'fro');
                      A(:,:,ii) = A(:,:,ii)/d(ii);              
                      for t = 1:K
                          C(:,(t-1)*N+1:t*N,ii,ii:R) = d(ii)*C(:,(t-1)*N+1:t*N,ii,ii:R);
                          C(:,(t-1)*N+1:t*N,1:ii,ii) = d(ii)*C(:,(t-1)*N+1:t*N,1:ii,ii);                  
                      end
                  end               
              end
        end
    case 2
       %% Parallelization in diagonalwise fashion
        while err>ERR && k<ITER
              k = k + 1;
              L = eye(N); L = repmat(L,[1,R]);
              U = eye(N); U = repmat(U,[1,R]);   
              %%
              % compute off-norm
              M = C;
              rr = 1;
              sum_off(k) = 0;
              sum_diag(k) = 0;
              for r1 = 1 : R
                  for r2 = r1 : R
                      Mtx = M(:,:,r1,r2);
                      [c1,c2,c] = offdiag(Mtx);
                      sum_off(k) = sum_off(k) + c1;
                      sum_diag(k)= sum_diag(k)+ c2; 
                      sub_off(k,rr) = c1;
                      sub_diag(k,rr) = c2;
                      rr = rr + 1;
                  end
              end
              %% U stage
                 for j = 2:N % ;  
                    CC = repmat((ones(N)-eye(N)),[1,K,R,R]);
                    CC = CC.*C;
                    W = reshape(permute(reshape(C,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    WC = reshape(permute(reshape(CC,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    ten_row_j = C(j:N,:,:,:); 
                    ten_row_i = CC(1:N+1-j,:,:,:);
                    ten_col_j = W(j:N,:,:,:); 
                    ten_col_i = WC(1:N+1-j,:,:,:);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[N+1-j,N*K,R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[N+1-j,N*K,R-ii+1]);
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[N+1-j,N*K,ii]);
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[N+1-j,N*K,ii]);
                        v = cat(3,v,cat(3,temp_col_j,temp_row_i));
                        w = cat(3,w,cat(3,temp_col_i,temp_row_j)); 
                        u = cat(3,u,cat(3,temp_col_j,temp_row_j));
                    end
                    f1 = reshape(squeeze(sum(v.*conj(w),2)),[size(v,1),R+1,R]);
                    f1 = reshape(sum(f1,2),[size(v,1),R]);
                    [~,q2,q3] = size(u);
                    qq = repmat((ones(N)-eye(N)),[1,q2/N,q3]);
                    f2 = reshape(squeeze(sum(u.*qq(1:N+1-j,:,:).*conj(u),2)),[size(u,1),R+1,R]);
                    f2 = reshape(sum(f2,2),[size(u,1),R]);
                    a = -f1./f2;
                    for ii = 1:R    
                        if abs(a(:,ii))>1, a(:,ii) = sign(a(:,ii))*1;end;
                        [~,a2,a3,a4] = size(C(j:N,:,ii,ii:R));
                        aa = repmat(a(:,ii),[1,a2,a3,a4]);
                        %% N by NK by R by R
                        C(1:N+1-j,:,ii,ii:R) = aa.*C(j:N,:,ii,ii:R) + C(1:N+1-j,:,ii,ii:R);  % update by left multiply 
                        I = [];
                        J = [];
                        for r1 = 1:N+1-j
                            I = [I,r1:N:r];
                        end
                        for r2 = j:N
                            J = [J,r2:N:r];
                        end
                        I = sort(I);
                        J = sort(J);
                        [~,~,a3,a4] = size(C(:,J,1:ii,ii));
                        C(:,I,1:ii,ii) = (repmat(a(:,ii)',[N,K,a3,a4])).*C(:,J,1:ii,ii) + C(:,I,1:ii,ii);% update by right multiply 
                    end
                    aa = reshape(permute(repmat(a,[1,1,N]),[1,3,2]),[size(a,1),N*R]);
                    U(1:N+1-j,:) = U(1:N+1-j,:) + aa.*U(j:N,:);           
                 end  

              %% L
              for i = 2:N         
                    CC = repmat((ones(N)-eye(N)),[1,K,R,R]);
                    CC = CC.*C;
                    W = reshape(permute(reshape(C,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    WC = reshape(permute(reshape(CC,[N,N,K,R,R]),[2,1,3,4,5]),[N,N*K,R,R]); 
                    ten_row_j = CC(i:N,:,:,:); 
                    ten_row_i = C(1:N+1-i,:,:,:);
                    ten_col_j = WC(i:N,:,:,:); 
                    ten_col_i = W(1:N+1-i,:,:,:);
                    v = [];
                    w = [];
                    u = [];
                    for ii = 1:R   
                        temp_row_j = reshape(ten_row_j(:,:,ii,ii:R),[N+1-i,N*K,R-ii+1]);
                        temp_row_i = reshape(ten_row_i(:,:,ii,ii:R),[N+1-i,N*K,R-ii+1]);
                        temp_col_j = reshape(ten_col_j(:,:,1:ii,ii),[N+1-i,N*K,ii]);
                        temp_col_i = reshape(ten_col_i(:,:,1:ii,ii),[N+1-i,N*K,ii]);
                        v = cat(3,v,cat(3,temp_col_i,temp_row_j));
                        w = cat(3,w,cat(3,temp_col_j,temp_row_i)); 
                        u = cat(3,u,cat(3,temp_col_i,temp_row_i));
                    end
                    f1 = reshape(squeeze(sum(v.*conj(w),2)),[size(v,1),R+1,R]);
                    f1 = reshape(sum(f1,2),[size(v,1),R]);
                    [~,q2,q3] = size(u);
                    qq = repmat((ones(N)-eye(N)),[1,q2/N,q3]);
                    f2 = reshape(squeeze(sum(u.*qq(i:N,:,:).*conj(u),2)),[size(u,1),R+1,R]);
                    f2 = reshape(sum(f2,2),[size(u,1),R]);
                    a = -f1./f2;
                    for ii = 1:R   
                        if abs(a(:,ii))>1, a(:,ii) = sign(a(:,ii))*1;end;
                        %% N by NK by R by R
                        [~,a2,a3,a4] = size(C(1:N+1-i,:,ii,ii:R));
                        aa = repmat(a(:,ii),[1,a2,a3,a4]);
                        C(i:N,:,ii,ii:R) = aa.*C(1:N+1-i,:,ii,ii:R) + C(i:N,:,ii,ii:R);  % update by left multiply    
                        I = [];
                        J = [];
                        for r1 = i:N
                            I = [I,r1:N:r];
                        end
                        for r2 = 1:N-i+1
                            J = [J,r2:N:r];
                        end
                        I = sort(I);
                        J = sort(J);
                        [~,~,a3,a4] = size(C(:,J,1:ii,ii));
                        C(:,I,1:ii,ii) = (repmat(a(:,ii)',[N,K,a3,a4])).*C(:,J,1:ii,ii) + C(:,I,1:ii,ii);% update by right multiply 
                    end             
                    aa = reshape(permute(repmat(a,[1,1,N]),[1,3,2]),[size(a,1),N*R]);
                    L(i:N,:) = L(i:N,:) + aa.*L(1:N+1-i,:);
              end

              for ii = 1:R           
                  A(:,:,ii) = L(:,(ii-1)*N+1:ii*N)*U(:,(ii-1)*N+1:ii*N)*A(:,:,ii); %
                  err_new(ii) = norm(L(:,(ii-1)*N+1:ii*N)*U(:,(ii-1)*N+1:ii*N)-eye(N,N),'fro')^2;
              end
              err = abs(err_old-max(err_new));
              err_old = max(err_new);       
              % Normalization     
              if rem(k,RBALANCE)==0
                  d = zeros(R,1);          
                  for ii = 1:R              
                      d(ii) = norm(A(:,:,ii),'fro');
                      A(:,:,ii) = A(:,:,ii)/d(ii);              
                      for t = 1:K
                          C(:,(t-1)*N+1:t*N,ii,ii:R) = d(ii)*C(:,(t-1)*N+1:t*N,ii,ii:R);
                          C(:,(t-1)*N+1:t*N,1:ii,ii) = d(ii)*C(:,(t-1)*N+1:t*N,1:ii,ii);                  
                      end
                  end               
              end
       end
end
Y = C;
A = reshape(permute(A,[1,3,2]),[N*R,N]);
P1 = [sub_off; sub_diag];
P2 = [sum_off; sum_diag];
end
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
    c = c + sum(mtx1(:));             % norm square sum of off-diag
    c1 = c1 + norm(diag(diag(mtx)))^2;% norm square sum of diag     
end
off_d = c; diag_d = c1; c = c/c1;
end
function [ISI,J_ISI] = f_ISI(A,W)
% A is the true loading matrices arranged as A = [A1; A2; A3;...AK]
% W is the estimated unloading matrices arranged as W = [W1; W2; ... WK]
[L,N] = size(A);
K = L/N;
ISI = zeros(K,1);
A = permute(reshape(A,[N,K,N]),[1,3,2]); 
W = permute(reshape(W,[N,K,N]),[1,3,2]);
abs_G_2 = zeros(N,N);
for kk = 1:K
    temp_A = A(:,:,kk);     
    temp_W = W(:,:,kk);
    G = temp_W * temp_A;
    abs_G = abs(G);
    abs_G_2 = abs_G_2 + abs_G;
    term_1 = abs_G./repmat(max(abs_G,[],2),[1,N]);
    term_1 = sum(sum(term_1,2) - ones(N,1));
    term_2 = abs_G./repmat(max(abs_G,[],1),[N,1]);
    term_2 = sum(sum(term_2,1) - ones(1,N));
    ISI(kk) = (term_1 + term_2)/(2*N*(N-1));    
end
term_1 = abs_G_2./repmat(max(abs_G_2,[],2),[1,N]);
term_1 = sum(sum(term_1,2) - ones(N,1));
term_2 = abs_G_2./repmat(max(abs_G_2,[],1),[N,1]);
term_2 = sum(sum(term_2,1) - ones(1,N));
J_ISI = (term_1 + term_2)/(2*N*(N-1));
end