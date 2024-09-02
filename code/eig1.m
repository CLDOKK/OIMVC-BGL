
function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)
%eigvec对应前c个最大特征值的特征向量 
%eigval前c个最大特征值
%eigval_full全部特征值，降序排列
if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end;

if nargin < 3
    isMax = 1;
    isSym = 1;
end;

if nargin < 4
    isSym = 1;
end;
%默认值 c = size(A,1); isMax = 1; isSym = 1.
if isSym == 1
    A = max(A,A');
end;
[v d] = eig(A); %v特征向量 d特征值
d = diag(d); %d本来是矩阵，取成向量
%d = real(d);
if isMax == 0
    [d1, idx] = sort(d);
else
    [d1, idx] = sort(d,'descend');
end;

idx1 = idx(1:c);
eigval = d(idx1); %前c个最大特征值
eigvec = v(:,idx1); %对应前c个最大特征值的特征向量

eigval_full = d(idx);

end