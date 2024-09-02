function [preY,P]= iteration_P(X,delta,Z,c)
  %Update P
    
    [m,n]=size(Z{1});
    View=length(X);
    sum_Z = sparse(n,m);
    for i1 = 1:View
        Z{i1}=Z{i1}'; 
        sum_Z = sum_Z + delta(i1) * Z{i1};
    end
    IterMaxP=50;
    [preY,~,P,~,~,~] = coclustering_bipartite_fast1(sum_Z, c, IterMaxP);
    P=P';
end

