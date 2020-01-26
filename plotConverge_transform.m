% Convergence data resave
clc;
close all;
clear ;
addpath(pwd);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)*randi(100,1)));
dbstop if error

datadir = 'F:/DATAt';
% addpath(genpath('RNN/BBOB09'));  % should point to fgeneric.m etc.
% rmpath(genpath('RNN/BBOB09/MCS'));
% opt.algName = 'PUT ALGORITHM NAME';
% opt.comments = 'PUT MORE DETAILED INFORMATION, PARAMETER SETTINGS ETC';
% maxfunevals = '1e4 * dim';  % SHORT EXPERIMENT, takes overall three minutes

algorithmSet = {
    'PSO_Bounds', ...
        'LPSO', ...
        'DMSPSO', ...
        'FIPS', ...
        'RLPSO', ...
        'new',...
    };

t0 = clock;
% global TYPE、
for TYPE =  1
    if TYPE == 1
        lower = -5;upper = 5;
        FUNNUMS = [9,10,11,12,13,14];
        FUN = 'fgeneric';
    else % our function
        lower = -20;upper = 20;
        FUNNUMS = 1:9;
        FUN = 'fgenericCBG';
    end
    
    
    for dim = [20]  % small dimensions first, for CPU reasons ,2,3,5,10,20,
        MAX_FES = 1e3 * dim;
        ERT = [];
        mERT = [];
        bestfStat = [];
        
        for ifun = FUNNUMS %[1:10]%1:24 %benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
            for mid = 1:numel(algorithmSet) % 调换了循环层级
                algorithm = algorithmSet{mid};
                
                maxRepeats = 3;
                bestf = [];
                num_succ = 0;
                num_unsucc = 0;
                succ_fes = [];
                unsucc_fes = [];
                convergeStat = zeros(maxRepeats*10, MAX_FES);
                
                for iinstance = 1:3 % [3,5,7]% [1:5, 1:5, 1:5]  % first 5 fct instances, three times                    
                    for repeat = 1:maxRepeats
                        T = 5; N = 10; s=1;
                        if TYPE == 1
                            datapath = sprintf('%s/%dD_BBOB_EXPERIMENTS/f%d/%s',datadir,dim,ifun,algorithm);
                            datafile = sprintf('%s/exp_%02u_%02u_f%d_DIM%d.tdat', datapath, iinstance, repeat, ifun, dim);
                        else
                            datapath = sprintf('%s/%dD_CBG_EXPERIMENTS/f%d/%s',datadir,dim,ifun,algorithm);
                            datafile = sprintf('%s/exp_%02u_%02u_f%d_DIM%d.tdat', datapath, iinstance, repeat, ifun, dim); % example: F:/DATA/2D_BBOB_EXPERIMENTS/f1/NELDER/exp_09_01_f4_DIM2.tdat
                        end
                        
                        if exist(datafile,'file')
                            load(datafile);
                            str1 = sprintf('convergeData = exp_%02u_%02u_f%d_DIM%d;', iinstance, repeat, ifun, dim);
                            str2 = sprintf('clear exp_%02u_%02u_f%d_DIM%d', iinstance, repeat, ifun, dim);
                            eval(str1);
                            eval(str2);
                            
                            converges = convergeData(:,[1,3]);
                            [m,D] = size(converges);
                            
                            A = [];
                            A(1,1) = converges(1,2);
                            k = 1;
                            if m > 1
                                for i = 2:m
                                    if converges(i,1)-converges(i-1,1) > 1 
                                        A(k,converges(i-1,1):converges(i,1)) = linspace(converges(i-1,2), converges(i,2), converges(i,1)-converges(i-1,1)+1);
                                    elseif converges(i,1)-converges(i-1,1) == 1 
                                        A(k,converges(i,1)) = converges(i,2);
                                    end
                                end
                            end
                            
                            converge = A;
                            
                            % 为画收敛曲线准备
                            if numel(converge) >= MAX_FES % 说明超过了最大FES， 要截断
                                Cconverge = converge(1:MAX_FES);
                            else % 说明提前达到最优就停止了
                                bestf = converge(end);
                                L = numel(converge);
                                tempC = zeros(MAX_FES,1);
                                tempC(1:L) = converge;
                                tempC(L:end) = bestf;
                                Cconverge = tempC;
                            end                             
                            convergeStat((iinstance-1) * maxRepeats + repeat, :) = Cconverge;
                            
                            % 为计算整个函数的ERT准备
                            [rows, cols] = find(converge>0);
                            if numel(rows) >= MAX_FES % 说明超过了最大FES， 要截断
                                num_unsucc = num_unsucc + 1;
                                unsucc_fes(num_unsucc) = MAX_FES;
                            else % 说明提前达到最优就停止了
                                num_succ = num_succ + 1;
                                succ_fes(num_succ) = numel(rows);
                            end
                            
                            % 最好适应值
                            if numel(converge) >= MAX_FES % 说明超过了最大FES， 要截断
                                bestfs((iinstance-1) * maxRepeats + repeat) = converge(MAX_FES);
                            else % 说明提前达到最优就停止了
                                bestfs((iinstance-1) * maxRepeats + repeat) = converge(end);
                            end
                            
                            outtime = formatTime(etime(clock, t0));
                            fprintf('TYPE%d, dim = %d, ifun = %d， method = %s, itrial = %d \t total time: %s\n',TYPE, dim, ifun, algorithm, iinstance, outtime);                            
                        else
                            error('file not found!');
                        end
                    end
                end % end of one method
                
                meanConv = mean(convergeStat);
                medianConv = median(convergeStat);
%                 L = numel(meanConv);
%                 for i = 1:L
%                     if meanConv(i) < 1e-6
%                         meanConv(i:L) = zeros(1,L-i+1);
%                         break;
%                     end
%                 end
%                 L = numel(medianConv);
%                 for i = 1:L
%                     if medianConv(i) < 1e-6
%                         medianConv(i:L) = zeros(1,L-i+1);
%                         break;
%                     end
%                 end
                dim_dir = sprintf('%s/%dD-TYPE%d-conv-mean-data',datadir,dim,TYPE);
                fun_dir = sprintf('%s/f%d',dim_dir,ifun);
                if ~exist(fun_dir, 'dir'); mkdir(fun_dir); end
                matname = sprintf('%s/%s.mat',fun_dir,algorithm);
                save(matname,'meanConv','medianConv');
                
                bestfStat(mid,ifun) = median(bestfs);
                
                if num_succ == 0
                    ERT(mid,ifun) = inf;
                else
                    if num_unsucc == 0
                        ERT(mid,ifun) = mean(succ_fes);
                    else
                        ERT(mid,ifun) = mean(succ_fes) + num_unsucc/num_succ * mean(unsucc_fes);
                    end
                end
                
            end % end of all methods
            
            mERT = ERT(:,ifun);
            
            if any(isinf(mERT))
                [rows,cols] = find(isinf(mERT)); % 重点关注
                for k = 1:numel(rows)
                    mmid = rows(k);
                    mERT(mmid,1) = 1e10 + bestfStat(mmid,ifun);
                end
            end
            
            [m,I] = sort(mERT);
            for i = 1:numel(algorithmSet)
                rankERTs(I(i),ifun) = i;
            end
            
            fprintf('---- T%d dimension %d-D, fun%d - all instances done, %s ---- \n', TYPE, dim, ifun, datetime);
        end
        
        matname = sprintf('%s/%dDERT-TYPE%d.mat',datadir,dim,TYPE);
        save(matname,'ERT','rankERTs','bestfStat');
        
        fprintf('-----***------T%d dimension %d-D, all functions done , %s------***------\n', TYPE, dim, datetime);
    end
    fprintf('Type%d done !!!! \n', TYPE);    
end
fprintf('over!!!\n');