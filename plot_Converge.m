% Convergence plots
clc;
close all;
clear ;
addpath(pwd);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)*randi(100,1)));
dbstop if error

datadir = 'F:/DATAt';
statdir =  'F:/STAT';
% addpath(genpath('RNN/BBOB09'));  % should point to fgeneric.m etc.

algorithmSet = {
    'PSO_Bounds', ...
        'LPSO', ...
        'DMSPSO', ...
        'FIPS', ...
        'RLPSO', ...
         'new',...
    };

t0 = clock;
% global TYPE
for TYPE =  1 %1:2
    if TYPE == 1
        lower = -5;upper = 5;
        FUNNUMS = [9,10,11,12,13,14];
        FUN = 'fgeneric';
    else % our function
        lower = -20;upper = 20;
        FUNNUMS = 1:9;
        FUN = 'fgenericCBG';
    end
    
    ERT = [];
    mERT = [];
    bestfStat = [];
    
    for dim = [20]  % small dimensions first, for CPU reasons ,2,3,5,10,20,
        MAX_FES = 1e3 * dim;
        for fid = 1:6%1:9 %[1:10]%1:24 %benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
            ifun = FUNNUMS(fid);
            for mid = 1:numel(algorithmSet) % 调换了循环层级
                algorithm = algorithmSet{mid};
                
                dim_dir = sprintf('%s/%dD-TYPE%d-conv-mean-data',datadir,dim,TYPE);
                fun_dir = sprintf('%s/f%d',dim_dir,ifun);
                if ~exist(fun_dir, 'dir'); mkdir(fun_dir); end
                matname = sprintf('%s/%s.mat',fun_dir,algorithm);
                %                 save(matname,'meanConv','medianConv');
                load(matname);               

                L = numel(meanConv);
                for i = 1:L
                    if meanConv(i) < 0
                        meanConv(i:L) = zeros(1,L-i+1);
                        break;
                    end
                end
                L = numel(medianConv);
                for i = 1:L
                    if medianConv(i) < 0
                        medianConv(i:L) = zeros(1,L-i+1);
                        break;
                    end
                end
                
                if strcmp(algorithm, 'PSO_Bounds')
                    legendmethod{mid} = 'PSO\_Bounds';
                else
                    legendmethod{mid} = algorithm;
                end
                if strcmp(algorithm, 'RLPSO')
                    s = 's'; c = 'red'; style ='-'; marker = 'none';  %
               
                elseif strcmp(algorithm, 'LPSO')
                    s = '^'; c = 'blue'; style =':'; marker = 'none';    
                elseif strcmp(algorithm, 'PSO_Bounds')
                    s = '^'; c = 'blue'; style ='--'; marker = 'none';
                elseif strcmp(algorithm, 'DMSPSO')
                   s = '.'; c = 'blue'; style ='-.';marker = 'none';
                elseif strcmp(algorithm, 'FIPS')
                    s = 's'; c = 'black'; style ='-'; marker = 'none';  %   
                elseif strcmp(algorithm, 'new')
                    s = 's'; c = 'red'; style ='-.'; marker = 'none';  %   
                end
     
%             end
                
                plot_interval = dim;
                %             semilogy(meanConv, 'linewidth',1.5);
                loglog(meanConv(1:plot_interval:end), 'LineWidth',3, 'Color', c, 'LineStyle', style );
                grid on
%                 loglog(meanConv(1:plot_interval:end), 'LineWidth',3, 'Color', c, 'LineStyle', style );
                %             plot(meanConv(1:plot_interval:end), 'LineWidth',1.5, 'Color', c, 'LineStyle', style );
                %             plot(meanConv);
                hold on;
%                 ylim([1e-6 1e1])
                %     set(gca,'XTick',[1,4,7,10]);
                %     set(gca,'XTickLabel',[300,1200,2100,3000]);
                set(gca,'FontSize',18,'Xcolor','k','Ycolor','k');
                if TYPE == 1
                    titlename = sprintf('f_{%d}',fid);
                else
                    titlename = sprintf('AF_{%d}',fid);
                end
                %         title(titlename);
                 xlabel('#FEs / D','FontSize',20,'Color','k');
                 ylabel('Mean best fitness','FontSize',20,'Color','k');
                %     xlabel('Generations','FontSize',20,'Color','k');
                %     ylabel('MSE','FontSize',20,'Color','k');
                %     legend('CGP','GP','Our')%,'GE','GP');
                if fid == 1
                    legend(legendmethod, 'FontSize',15,'TextColor','k','northeast');
                end
                
            end
            
%             capben = sprintf('f_{%d}', ifun);
            title(titlename, 'FontSize',22,'Color','k');
            hold off;
            convFigs = sprintf('%s/convFigs/TYPE%d-%dD',statdir,TYPE,dim);
            if ~exist(convFigs, 'dir'); mkdir(convFigs);end
            pname=sprintf('%s\\f%d.pdf',convFigs, ifun);
            saveas(gcf,pname);
            
%             pname=sprintf('%s\\f%d.pdf',convFigs, ifun);
%             saveas(gcf,pname);
            
            fprintf('---- T%d dimension %d-D,  - all functions done, %s ---- \n', TYPE, dim, datetime);
        end
        
        fprintf('-----***------T%d dimension %d-D, all functions done , %s------***------\n', TYPE, dim, datetime);
    end
    fprintf('Type%d done !!!! \n', TYPE);
end
fprintf('over!!!\n');

