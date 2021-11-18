function [Y,B,P] = f_LUCJD(X,opt)
%% This function realizes LU decomposition based Complex-valued 
%  non-orthogonal Joint Diagonalization (LUCJD) proposed in:
%  Ke Wang, Xiao-Feng Gong, Qiu-Hua Lin, "Complex Non-Orthogonal Joint
%  Diagonalization Based on LU and LQ Decomposition", Proc. LVA/ICA2012
%% Inputs
%  X: Target matrices arranged as X = [X1,X2,...,XN]
%  opt: options for parallelization
%       0 - without parallelization
%       1 - columnwise parallelization
%       2 - diagonalwise parallelization
%% Outputs
%  Y: The unloaded target matrices
%  B: The unmixing matrix
%  P: Performance indices for all sweeps
%% Information 
%  Copyright: Xiao-Feng Gong, Dalian University of Technology
%  Author 	: Xiao-Feng Gong: xfgong@dlut.edu.cn; Xiu-Lin Wang
%  Date     : 2014-12-30
%  Citation : Please cite the above references if you use this program. 
%  Addition : Comments, bug reports, etc are welcome.
%% 
[n,m] = size(X);
N = m/n;       % N is the number of target matrices, n is the dimensionality of matrices 
tol = 1*10^-8; % Stopping threshold
Niter = 200;   % The maximal number of iterations
B = eye(n,n);  % Joint diagonalizer
err = tol*n+1; % Error 
k = 0;         % Number of sweeps
%% 
err_old = 0;
switch opt
    case 0
       %% LUCJD in its original sequential version
        while err>tol && k<Niter % 
            k = k+1;
            [c1(k),c2(k),c3(k)] = offdiag(X);
            L = eye(n);
            U = eye(n);
            % U stage 
            for i = 1:n
                for j = i+1:n
                    cindex = 1:m;
                    cindex(i:n:m) = [];
                    W = [];
                    for g = 1:N
                        W = [W X(:,(g-1)*n+1:g*n).'];
                    end
                    a = -((X(i,cindex)*X(j,cindex)')+W(j,cindex)*W(i,cindex)')/(X(j,cindex)*X(j,cindex)'+W(j,cindex)*W(j,cindex)');% Optimal parameter for current index pair
                    if abs(a)>1, a=sign(a)*1;end;
                    I = i:n:m;
                    J = j:n:m;
                    X(i,:) = a*X(j,:)+X(i,:); 
                    X(:,I) = a'*X(:,J)+X(:,I); % Update X
                    U(i,:) = U(i,:)+a*U(j,:);  % Update U
                end
            end            
            % L stage
            for i = 1:n
                for j = i+1:n
                    cindex = 1:m;
                    cindex(j:n:m) = [];
                    W = [];
                    for g = 1:N
                        W = [W X(:,(g-1)*n+1:g*n).'];
                    end
                    a = -((X(j,cindex)*X(i,cindex)')+W(i,cindex)*W(j,cindex)')/(X(i,cindex)*X(i,cindex)'+W(i,cindex)*W(i,cindex)');
                    if abs(a)>1, a = sign(a)*1; end;
                    I = i:n:m;
                    J = j:n:m;
                    X(j,:) = a*X(i,:)+X(j,:);                    
                    X(:,J) = a'*X(:,I)+X(:,J);
                    L(j,:) = L(j,:)+a*L(i,:);
                end
            end                  
            B = L*U*B; % Update of B after U and L stages are completed
            err_new = norm(L*U-eye(n,n),'fro')^2;
            err = abs(err_old-err_new);
            err_old = err_new;
            % Row balancing
            d = sum(abs(X'));
            D = diag(1./d*N); 
            for t = 1:N
                X(:,(t-1)*n+1:t*n) = D*X(:,(t-1)*n+1:t*n)*D;
            end;
            B = D*B;            
        end            
    case 1
       %% Parallelization in columnwise fashion
        while err>tol && k< Niter 
            k = k+1;
            L = eye(n);
            U = eye(n);
            [c1(k),c2(k),c3(k)] = offdiag(X);
            % U stage 
            for i = 2:n
                C = X;
                for ii = 1:n
                    cindex = ii:n:m;
                    C(ii,cindex) = 0;
                end
                WC = reshape(permute(reshape(C,[n,n,N]),[2,1,3]),[n,n*N]);
                W = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                reX = repmat(X(i,:),[i-1,1]);
                reW = repmat(W(i,:),[i-1,1]);
                a1 = sum(C(1:i-1,:).*conj(reX),2);                                 
                a2 = sum(reW.*conj(WC(1:i-1,:)),2);
                cc = ones(n)-eye(n);
                cc = repmat(cc,[1,N]);
                a3 = sum(reX.*cc(1:i-1,:).*conj(reX),2);
                a4 = sum(reW.*cc(1:i-1,:).*conj(reW),2);
                a = -(a1+a2)./(a3+a4);
                if abs(a)>1, a = sign(a)*1;end;
                X(1:i-1,:) = repmat(a,[1,m]).*reX+X(1:i-1,:);
                CQ = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);                
                CQ(1:i-1,:) = repmat(conj(a),[1,m]).*repmat(CQ(i,:),[i-1,1])+CQ(1:i-1,:);
                X = reshape(permute(reshape(CQ,[n,n,N]),[2,1,3]),[n,n*N]);
                U(1:i-1,:) = U(1:i-1,:)+repmat(a,[1,n]).*repmat(U(i,:),[i-1,1]);
            end            
            % L stage 
            for i = 1:n-1
                C = X;
                for ii = 1:n
                    cindex = ii:n:m;
                    C(ii,cindex) = 0;
                end
                WC = reshape(permute(reshape(C,[n,n,N]),[2,1,3]),[n,n*N]);
                W = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                reX = repmat(X(i,:),[n-i,1]);
                reW = repmat(W(i,:),[n-i,1]);
                a1 = sum(C(i+1:n,:).*conj(reX),2);                                 
                a2 = sum(reW.*conj(WC(i+1:n,:)),2);
                cc = ones(n)-eye(n);
                cc = repmat(cc,[1,N]);
                a3 = sum(reX.*cc(i+1:n,:).*conj(reX),2);
                a4 = sum(reW.*cc(i+1:n,:).*conj(reW),2);
                a = -(a1+a2)./(a3+a4);
                if abs(a)>1, a = sign(a)*1; end;
                X(i+1:n,:) = repmat(a,[1,m]).*reX+X(i+1:n,:);
                CQ = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                CQ(i+1:n,:) = repmat(conj(a),[1,m]).*repmat(CQ(i,:),[n-i,1])+CQ(i+1:n,:);
                X = reshape(permute(reshape(CQ,[n,n,N]),[2,1,3]),[n,n*N]);
                L(i+1:n,:) = L(i+1:n,:)+repmat(a,[1,n]).*repmat(L(i,:),[n-i,1]);
            end            
            B = L*U*B;
            err_new = norm(L*U-eye(n,n),'fro')^2;
            err = abs(err_old-err_new);
            err_old = err_new;
            % Row balancing 
            d = sum(abs(X')); 
            D = diag(1./d*N);
            for t = 1:N
                X(:,(t-1)*n+1:t*n) = D*X(:,(t-1)*n+1:t*n)*D;
            end;
            B = D*B;            
        end
    case 2
       %% Parallelization in diagonal wise fashion
        while err>tol && k<Niter % 
            k=k+1;
            L=eye(n);     
            U=eye(n);
            [c1(k),c2(k),c3(k)] = offdiag(X);
            % U stage
            for i=2:n
                C = X;
                for ii = 1:n
                    cindex = ii:n:m;
                    C(ii,cindex) = 0;
                end                
                WC = reshape(permute(reshape(C,[n,n,N]),[2,1,3]),[n,n*N]);
                W = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                a1 = sum(C(1:n+1-i,:).*conj(X(i:n,:)),2);                                 
                a2 = sum(W(i:n,:).*conj(WC(1:n+1-i,:)),2);
                cc = ones(n)-eye(n);
                cc = repmat(cc,[1,N]);
                a3 = sum(X(i:n,:).*cc(1:n+1-i,:).*conj(X(i:n,:)),2);
                a4 = sum(W(i:n,:).*cc(1:n+1-i,:).*conj(W(i:n,:)),2);
                a = -(a1+a2)./(a3+a4);
                if abs(a)>1, a = sign(a)*1;end;
                X(1:n+1-i,:) = repmat(a,[1,m]).*X(i:n,:)+X(1:n+1-i,:);               
                CQ = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                CQ(1:n+1-i,:) = repmat(conj(a),[1,m]).*CQ(i:n,:)+CQ(1:n+1-i,:);
                X = reshape(permute(reshape(CQ,[n,n,N]),[2,1,3]),[n,n*N]);
                U(1:n+1-i,:) = U(1:n+1-i,:)+repmat(a,[1,n]).*U(i:n,:);
            end            
            % L stage
            for i = 2:n
                C = X;
                for ii = 1:n
                    cindex = ii:n:m;
                    C(ii,cindex) = 0;
                end
                WC = reshape(permute(reshape(C,[n,n,N]),[2,1,3]),[n,n*N]);
                W = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                a1 = sum(C(i:n,:).*conj(X(1:n+1-i,:)),2);                                 
                a2 = sum(W(1:n+1-i,:).*conj(WC(i:n,:)),2);
                cc = ones(n)-eye(n);
                cc = repmat(cc,[1,N]);
                a3 = sum(X(1:n+1-i,:).*cc(i:n,:).*conj(X(1:n+1-i,:)),2);
                a4 = sum(W(1:n+1-i,:).*cc(i:n,:).*conj(W(1:n+1-i,:)),2);
                a = -(a1+a2)./(a3+a4);
                if abs(a)>1, a = sign(a)*1; end;
                X(i:n,:) = repmat(a,[1,m]).*X(1:n+1-i,:)+X(i:n,:);
                CQ = reshape(permute(reshape(X,[n,n,N]),[2,1,3]),[n,n*N]);
                CQ(i:n,:) = repmat(conj(a),[1,m]).*CQ(1:n+1-i,:)+CQ(i:n,:);
                X = reshape(permute(reshape(CQ,[n,n,N]),[2,1,3]),[n,n*N]);
                L(i:n,:) = L(i:n,:)+repmat(a,[1,n]).*L(1:n+1-i,:);
            end            
            B = L*U*B;
            err_new = norm(L*U-eye(n,n),'fro')^2;
            err = abs(err_old-err_new);
            err_old = err_new;            
            d = sum(abs(X')); 
            D = diag(1./d*N); 
            for t = 1:N
                X(:,(t-1)*n+1:t*n) = D*X(:,(t-1)*n+1:t*n)*D;
            end;
            B = D*B;                  
        end
end
Y = X;
P = [c1;c2;c3];
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
    c = c + sum(mtx1(:));%sum of squared norms of off-diagonal elements
    c1 = c1 + norm(diag(diag(mtx)))^2;%sum of squared norms of diagonal elements
end
off_d = c; diag_d = c1; c = c/c1;
end