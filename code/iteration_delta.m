function delta= iteration_delta(X,Z,P,m)
View=length(X);
n=size(X{1},2);   
   %Update delta
    for iv = 1:View
        Z{iv} = sparse(Z{iv});
    end
    ZZ = sparse(n * m,View);
    for i2 = 1:View
        ZZ(:,i2) = reshape(Z{i2},[n*m 1]);
    end
    newZ = ZZ'*ZZ;
    p = reshape(P,[n*m 1]);
    s = 2*ZZ'*p;
    %QP_options
    QP_options = 'quadprog';
    %QP_options = 'SimplexQP_acc';
    switch lower(QP_options)
        case {lower('SimplexQP_acc')}
            delta = SimplexQP_acc(newZ, s);
            
        case {lower('quadprog')}
            options = optimset( 'Algorithm','interior-point-convex','Display','off');
            delta = quadprog(2*newZ,-s,[],[],ones(1,View),1,zeros(View,1),ones(View,1),[],options);            
    end 
end

