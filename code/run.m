clear;
clc;
addpath('measure');
path = '../';
addpath(path);
resultdir1 = 'Results/';
if (~exist('Results', 'file'))
    mkdir('Results');
    addpath(genpath('Results/'));
end

resultdir2 = 'totalResults/';
if (~exist('totalResults', 'file'))
    mkdir('totalResults');
    addpath(genpath('totalResults/'));
end
addpath('');% datasets 
datadir='';% datasets

%%dataname={'MSRCV1_3v', 'yale_mtv_2', 'ORL','COIL20-3v','handwritten_3v'};
% dataname = {'BBCSport','Caltech101-7','Caltech101-20','CCV','cifar10','cifar100','handwritten','NUSWIDEOBJ','ORL_mtv','proteinFold','Reuters','SUNRGBD_fea','YoutubeFace_sel_fea'};
dataname = {'Caltech101-7','NUSWIDEOBJ','prokaryotic','proteinFold','Wiki_fea'};
numdata = length(dataname); % number of the test datasets
numname = {'_Per0.1', '_Per0.2', '_Per0.3', '_Per0.4','_Per0.5', '_Per0.6', '_Per0.7', '_Per0.8', '_Per0.9','_Per1'};
%
% dataname = {'3sourceIncomplete', 'bbcsportIncomplete','Mnist_5K_Per0.5'};
% numdata = length(dataname); % number of the test datasets
% numname = {''};
lee=length(dataname);
for idata =3
    ResBest = zeros(10, 8);
    ResStd = zeros(10, 8);
    % result = [Fscore Precision Recall nmi AR Entropy ACC Purity];
    for dataIndex =  1:1:10
        %     for dataIndex = 1 : 1
        datafile = [datadir, cell2mat(dataname(idata)), cell2mat(numname(dataIndex)), '.mat'];
        load(datafile);
        %data preparation...
        gt = truelabel{1};
        cls_num = length(unique(gt));
        k= cls_num;
        tic;
        %%[X1, ind] = DataPreparing(data, index);
        [X1, ind] = findindex(data, index);
        
        
        time1 = toc;
        maxAcc = 0;
        %         TempLambda1 = 0.01;
        
        %%%%%%%%%%%%%%%%%%%%%%参数%%%%%%%%%%%%%%%%
        pp=-3:0.5:3;
        TempLambda1 = [10.^pp];
        %TempLambda1=0.1;
        %         TempLambda1 = [0.001,0.1,1,10];
        TempLambda2 = [10.^pp];
        %TempLambda2=0.001;
        
        annum=[k,2*k,4*k,6*k,8*k];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        idx = 1;
        lamget=zeros(1,3);
        runtimeget=0.1;
        
        for LambdaIndex1 = 1 : length(TempLambda1)
            lambda1 = TempLambda1(LambdaIndex1);
            
            for LambdaIndex2 = 1 : length(TempLambda2)
                lambda2 = TempLambda2(LambdaIndex2);
                for ann=1:length(annum)
                    numanchor=annum(ann);
disp([char(dataname(idata)), char(numname(dataIndex)), '-lam1=', num2str(lambda1),'-lam2=', num2str(lambda2),'-lam3=', num2str(ann)]);
                    tic;
                    [preY,V,A,W,Z,iter,obj] = missingalgo_lhc(X1,gt,lambda1,lambda2,k,numanchor,ind); % X,Y,lambda,d,numanchor
                    time2 = toc;

                    res = Clustering8Measure(gt, preY);
                    runtime(idx) = time2;
                    disp(['runtime:', num2str(runtime(idx))])
                    idx = idx + 1;
                    
                    tempResBest(dataIndex, : ) = res;
                    tempResStd(dataIndex, : ) = res;
                    
                    if tempResBest(dataIndex, 1) > ResBest(dataIndex, 1)
                        lamget(1)=LambdaIndex1;
                        lamget(2)=LambdaIndex2;
                        lamget(3)=ann;                     
                    end
                    
                    for tempIndex = 1 : 8
                        if tempResBest(dataIndex, tempIndex) > ResBest(dataIndex, tempIndex)
                            ResBest(dataIndex, tempIndex) = tempResBest(dataIndex, tempIndex);
                        end
                    end
                    
                    save([resultdir1, char(dataname(idata)), char(numname(dataIndex)), '-l1=', num2str(lambda1), '-l2=', num2str(lambda2), ...
                        '-l3=', num2str(ann), '_result.mat'], 'tempResBest','time2');
                    
                end
            end
        end
        runtimeget=mean(runtime);
        PResBest = ResBest(dataIndex, :);
        save([resultdir2, char(dataname(idata)), char(numname(dataIndex)), '_result.mat'], 'PResBest','lamget','runtimeget');
        
    end
    
    
end
