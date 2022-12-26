function W = SNMTF_Esraa(X, r,Si,Hi,a, max_iter, tol)
%
% compute SNMTF [2]
% input:
%   X          nonnegative data input (m times n)
%   r          number of SNMTF components
%   max_iter   maximum number of iterations (defaut 5000)
%   tol        convergence tolerance (default 1e-5)
% output:
%   W          the factorizing matrix (m times r)
%
% [1]Chris Ding, Tao Li, Wei Peng, and Haesun Park. 2006. 
%     Orthogonal nonnegative matrix t-factorizations for clustering. 
%     In KDD '06. ACM, New York, NY, USA, 126-135.
% [2] Non-Negative Matrix Factorizations for Multiplex
%Network Analysis
%Vladimir Gligorijevi´c, Yannis Panagakis, and Stefanos Zafeiriou, Member, IEEE

a=0.8;
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter =500;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
check_step = 100;



%initialize W
Aavg_i=zeros(size(X{1,1},1),size(X{1,1},1));
for i=1:size(X,2)
Aavg_i=Aavg_i+X{1,i};
end
Aavg_i=Aavg_i/size(X,3);
[~,~,W,~,~,~]=orthnmfrule_Esraa(Aavg_i,r);

% m = size(X,1);
% W = rand(m,r);

% XX = X * X';
for iter=1:max_iter
    W_old = W;
    if mod(iter,check_step)==0
%         fprintf('iter=% 5d ', iter);
    end
    
%     W = W .* (XX*W) ./ (W*(W'*XX*W) + XX*W*(W'*W));
num=zeros(size(W));
% den=zero(size(W));
for j=1:size(X,2)
num=num+(X{1,j}*W*Si{1,j}+(a/2)*(Hi{1,j}*Hi{1,j}')*W);
end

    W = W .* (num) ./ ((W*W')*num);
%     W = W ./ norm(W);
    
    diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    if diffW<tol
        fprintf('converged after %d steps.\n', iter);
        break;
    end
    
%     if mod(iter,check_step)==0
%         fprintf('diff=%.10f, ', diffW);
%         fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
%         fprintf('\n');
%     end
end
