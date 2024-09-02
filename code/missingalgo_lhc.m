function [preY,V,A,W,P,iter,obj] = missingalgo_lhc(X,Y,lambda,beta,d,numanchor,ind)
% m      : the number of anchor. the size of Z is m*n.
% lambda : the hyper-parameter of regularization term.

% X      : n*di

%% initialize
maxIter = 50 ; % the number of iterations
m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);
delta = 1 / numview* ones(numview,1);
W = cell(numview,1);  % di * d
A = cell(numview,1);  % d  * m
Z = cell(numview,1);  % m  * n
missingindex = constructA(ind);


for i = 1:numview
    %%X{i} = mapstd(X{i}',0,1); % turn into d*n
    di = size(X{i},1);
    W{i} = zeros(di,d);
    A{i} = zeros(d,m);   % d  * m
    Z{i}= zeros(m,numsample); % m  * n
    Z{i}(:,1:m) = eye(m);
    delta = 1 / numview * ones(numview,1);
end
P = delta(1) * Z{1}; %Init P
for iv = 2:numview
    P = P + delta(iv) * Z{iv};
end


opt.disp = 0;

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;
    %disp(['number of iter:',num2str(iter)])
    %% optimize W_i
    for iv=1:numview
        C = X{iv}*Z{iv}'*A{i}';
        [U,~,V] = svd(C,'econ');
        W{iv} = U*V';
    end
   
    %% optimize A2
    
    for ia = 1:numview
        part1 = W{ia}' * X{ia} * Z{i}';
        [Unew,~,Vnew] = svd(part1,'econ');
        A{ia} = Unew*Vnew';
    end

    
   
    %% optimize P
    [preY,P]=iteration_P(X,delta,Z,numclass);
   
    
    
    
    %% optimize1 Z
    sum_qh=delta(1)*Z{1};
    for i=2:numview
        sum_qh=sum_qh+delta(i)*Z{i};
    end
    for i=1:numview
        C2=A{i}'*W{i}'*X{i};
        C2=C2-beta*delta(i)*(sum_qh-Z{i}-P);
        C1=ind(:,i)'+lambda*ones(1,numsample)+beta*delta(i)*delta(i)*ones(1,numsample);
        for ii=1:numsample
            idx=1:numanchor;
            ut=C2(idx,ii)./C1(ii);
            Z{i}(idx,ii)=EProjSimplex_new(ut');
            
        end
        
    end
    
    %
%     [UU,~,V]=svd(Z','econ');
%     ts{iter} = UU(:,1:numclass);
    
    %         C1=0;
    %         C2=0;
    %         for a=1:numview
    %             C1=C1+alpha(a)^2*ind(:,a)';%%ind是用来看是不是选这个
    %             C2=C2+alpha(a)^2*A'*W{a}'*X{a};
    %         end
    %         C1=C1+lambda*ones(1,numsample);
    %
    %         for ii=1:numsample
    %             for idx=1:numanchor
    %                 ut=C2(idx,ii)./C1(ii);
    %                 Z(idx,ii)=EProjSimplex_new(ut');
    %             end
    %         end
    

    %% optimize delta
    delta=iteration_delta(X,Z,P,m);
 
    %%
    term1 = 0;
    for iv = 1:numview
        term1 = term1 + norm(X{iv} - W{iv} * A{iv} * (Z{iv}.*repmat(missingindex{iv},m,1)),'fro')^2+beta*norm((sum_qh-P),'fro')^2;
    end
    obj(iter) = term1;
    
    
% % % %     if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
% % % %         [UU,~,V]=svd(Z','econ');
% % % %         UU = UU(:,1:numclass);
% % % %         flag = 0;
% % % %     end
    if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
         flag = 0;
    end
    
end





end


