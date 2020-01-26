function M_GPSO11(FUN, dim, lower, upper, ftarget, maxfunevals)
% clc;
% clear ;
% close all;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)*randi(100,1)));
global nname 
global FUN_NUM FITCOUNT Fvalues MAX_FES %ftarget
global xx yy z IS_PLOT

% addpath(genpath('RNN/BBOB09'));  % should point to fgeneric.m etc.
% rmpath(genpath('RNN\BBOB09\MCS'));
dbstop if error

DIMENSION = dim;
opt.repeat = 1;
% funcID = 9;
% ftarget = fgeneric('initialize',... 
%                     funcID,... % funID 
%                     2 ,... % dim
%                     1,... % inst
%                     1, 'F:\\STAT', opt);
ft = 1e-8;
nname=sprintf('D:\\桌面\\ANN优化\\递归新\\DATA\\nn\\111122.mat');
FITCOUNT = 0;
Fvalues = [];
FUN_NUM = 11;

POP_SIZE = 40;

INIT_MIN = lower;INIT_MAX = upper;
MAX_FES = maxfunevals;%10000*DIMENSION;
MAX_GENERATION=floor(MAX_FES/POP_SIZE);
PLOT_INTERVAL=POP_SIZE*5;
iwt=0.9-(1:MAX_GENERATION)*(0.5/MAX_GENERATION);%inertia weight

% c_sigma = 0.5-(1:MAX_GENERATION)*(0.01/MAX_GENERATION);
c_sigma = 0.1-(1:MAX_GENERATION)*(0.1/MAX_GENERATION);
v_sigma = 1.5-(1:MAX_GENERATION)*(1.5/MAX_GENERATION);
s_sigma = 100-(1:MAX_GENERATION)*(100/MAX_GENERATION);

xinit = rand(1,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;
% xinit = ones(1,DIMENSION)*4 + normrnd(0,0.3,[1,DIMENSION]);
% x = repmat(xinit,POP_SIZE,1) + 1 * normrnd(0,1, [POP_SIZE,DIMENSION]);

x=rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;

v=(rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;

% initialize(DIMENSION,50,INIT_MAX,INIT_MIN);

IS_PLOT = 0;

if DIMENSION>2
    IS_PLOT=0;
end

if IS_PLOT
    [pos,xx,yy]=PositionMatrix(INIT_MAX,INIT_MIN,99);
    d=length(yy);
    f=fitnessf(pos);
    z=vec2mat(f,d)';
    FITCOUNT = 0;
    
    figure(100);
    % colormap jet;
    % shading interp;
    contour(xx,yy,z,50);
    hold on;
    plot(x(:,1),x(:,2),'k.','markersize',10);
    pause(0.5)
end

f=fitnessf(x);
fitcount=POP_SIZE;

pv=zeros(POP_SIZE,DIMENSION);
ppv=zeros(1,DIMENSION);
pf = ones(1,POP_SIZE)*1.0E100;
ppf = 1.0E100;


p=zeros(POP_SIZE,DIMENSION);
pp=zeros(1,DIMENSION);
generation=zeros(floor(MAX_FES/PLOT_INTERVAL),1);
bestf = inf;
p_delta_f = ones(1,POP_SIZE)*1.0E100;
pp_delta_f = 1.0E100;

x=rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;
v=(rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;

X=[];
F(1,:) = f;
meanTheta = [];
stuck = 0;
unstuck = 0;
for g = 1:MAX_GENERATION
%     if(fitcount>MAX_FES || bestf < ft); break;end; 
      if(fitcount>MAX_FES) || feval(FUN, 'fbest') < ftarget; break;end
     f = feval(FUN, x')';
%     f=fitnessf(x);
    fitcount=fitcount+POP_SIZE;
    F(g,:) = f;
    if g > 10
        delta_f = F(g,:) - min(F(g-10:g,:));
    else
        delta_f = F(g,:) - min(F(1:g,:));
    end
    
%     f_old = f;
    if mod(g,5) == 0
        p_delta_f = ones(1,POP_SIZE)*1.0E100;
        pp_delta_f = 1.0E100;        
    end
    
    if mod(g,5) == 0
        pf = ones(1,POP_SIZE)*1.0E100;
        ppf = 1.0E100;        
    end
    
    X=[X;x];
%     F=[F;f];
    bestflag = 1;
    for j=1:POP_SIZE
        
        if delta_f(j) < p_delta_f(j)
            p(j,:) = v(j,:);
            p_delta_f(j) = delta_f(j);
                        
            if p_delta_f(j) < pp_delta_f
                pp(:)=p(j,:);
                pp_delta_f = p_delta_f(j);
                bestvi = j;             
            end
        end
    end
    
    for j=1:POP_SIZE
        if(f(j)<pf(j))
            pv(j,:) = v(j,:);
            pf(j)=f(j);
            if(pf(j)<ppf)
                ppv(:)=pv(j,:);
                ppf=pf(j); 
                
                xbest = x(j,:);
                bestflag = 0;
                besti = j;
        
            end
        end
    end
    
%     [tempppf,besti] = min(f);
%         
%     if tempppf < ppf
%         ppf = tempppf;
%         xbest = x(besti,:);
%         bestflag = 0;
%     end
    
    if bestflag
        stuck = stuck + 1;
        unstuck = 0;
    else        
        unstuck = unstuck + 1;
        stuck = 0;
    end
    
    if g == 1
        xbestold = xbest;
    end
    
    c_cov = 1/(3 *sqrt(DIMENSION) + 5);
    
    T = randi(POP_SIZE);
    sigma = 0.3;
    for j=1:POP_SIZE
        
        % 2 + 4 慢，探索更多
        % 2 + 2 快，易早熟
        % sigma 很重要
%         c = 2*(pp/norm(pp)-p(j,:)); 
%         c = symmetricVec3(s_sigma(g), p(j,:), pp);
        c = symmetricVec3(s_sigma(g), v(j,:), p(j,:));
        P1 = xbest+v(j,:);
        P2 = xbest+p(j,:);
        P3 = xbest+c;
        VC = (P1 + P2 + P3)/3; % 几合中心
        radius = (norm(P1-VC) + norm(P2-VC) + norm(P3-VC)) / 3;
        VC_shift = VC + normrnd(0,c_sigma(g),[1,DIMENSION])*radius;
        vL1 = VC_shift-xbest;
        
        c = symmetricVec3(s_sigma(g), v(j,:), pp);
        P1 = xbest+v(j,:);
        P2 = xbest+pp;
        P3 = xbest+c;
        VC = (P1 + P2 + P3)/3; % 几合中心
        radius = (norm(P1-VC) + norm(P2-VC) + norm(P3-VC)) / 3;
        VC_shift = VC + normrnd(0,c_sigma(g),[1,DIMENSION])*radius;
        vL2 = VC_shift-xbest;
        
        c = symmetricVec3(s_sigma(g), v(j,:), pv(j,:));
        P1 = xbest+v(j,:);
        P2 = xbest+pv(j,:);
        P3 = xbest+c;
        VC = (P1 + P2 + P3)/3; % 几合中心
        radius = (norm(P1-VC) + norm(P2-VC) + norm(P3-VC)) / 3;
        VC_shift = VC + normrnd(0,c_sigma(g),[1,DIMENSION])*radius;
        vL3 = VC_shift-xbest;
        
        c = symmetricVec3(s_sigma(g), v(j,:), ppv);
        P1 = xbest+v(j,:);
        P2 = xbest+ppv;
        P3 = xbest+c;
        VC = (P1 + P2 + P3)/3; % 几合中心
        radius = (norm(P1-VC) + norm(P2-VC) + norm(P3-VC)) / 3;
        VC_shift = VC + normrnd(0,c_sigma(g),[1,DIMENSION])*radius;
        vL4 = VC_shift-xbest;
        
%         c = (2*(p(j,:)*pp'/norm(pp)^2)*pp-p(j,:))*((norm(pp)+norm(p(j,:)))/(2*norm(p(j,:)))); 
%         v(j,:) = rand*v(j,:) + rand * pp + rand * p(j,:); %1
%         v(j,:) = sigma * normrnd(0,1,[1,1]) * v(j,:) + rand * c + rand * p(j,:); %2
%         v(j,:) = v_sigma(g) * normrnd(0,1,[1,1]) * v(j,:) + rand(1,DIMENSION) .* vL1 + rand(1,DIMENSION) .* vL2 + rand(1,DIMENSION) .* vL3 + rand(1,DIMENSION) .* vL4; %2
        v(j,:) = v_sigma(g) * normrnd(0,1,[1,1]) * v(j,:) + (rand(1,DIMENSION) .* vL1 + rand(1,DIMENSION) .* vL2 + rand(1,DIMENSION) .* vL3 + rand(1,DIMENSION) .* vL4)/4; %2
%         v(j,:) = v_sigma(g) * normrnd(0,1,[1,1]) * v(j,:) + rand(1,DIMENSION) .* vL3 + rand(1,DIMENSION) .* vL4; %2

%         v(j,:) = rand * v(j,:) + rand(1,DIMENSION) .* vL1 + rand(1,DIMENSION) .* vL2; %2.1
%         v(j,:) = normrnd(0,v_sigma(g),[1,DIMENSION]) .* v(j,:) + vL; %2.5
%         v(j,:) = rand * v(j,:) + vL; %2
%         v(j,:) = normrnd(0,0.3,[1,1]) * v(j,:) + rand * c + rand * p(j,:) + rand * (xbestold - x(j,:)); %3
%         v(j,:) = normrnd(0,0.3,[1,DIMENSION]) .* v(j,:) + rand * c + rand * p(j,:); %4
        
        if stuck > 1
            s = mean(abs(v(j,:)),2)*10;
            v(j,:) = v(j,:) + normrnd(0,s,[1,DIMENSION]); %4
            stuck = 0;
        end
        if unstuck > 1 % && (j == besti || j == bestvi)
            v(j,:) = v(j,:)*2; %5
            unstuck = 0;
        end
        
        for k = 1:DIMENSION
            if abs(v(j,k)) > INIT_MAX*0.2
                vsign = sign(v(j,k));
                vecV = abs(reshape(v,1,POP_SIZE*DIMENSION));
                [row, col] = find(vecV < INIT_MAX*0.2);
                vecV = vecV(col);
                m = mean(vecV) + std(vecV) * randn * 10;
%                     v(j,k) = normrnd(0,1,[1,1]);
                v(j,k) = vsign * m; %(rand(1,1)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;
            end
        end
   
        
%         x(j,:) = x(j,:) + normrnd(0,0.3,[1,1])*v(j,:); %0
        x(j,:) = xbest + v(j,:); %1
%         x(j,:) = xbest + v_sigma(g) * normrnd(0,1,[1,1]) * v(j,:); %2    
%         x(j,:) = xbest +  normrnd(0,0.3,[1,DIMENSION]) .* v(j,:) ;%3
%         x(j,:) = xbest +  normrnd(0,0.3,[1,DIMENSION]) .* v(j,:) + rand * (x(j,:)-xbestold );%4
%         x(j,:) = xbest + normrnd(0,0.3,[1,1]) * v(j,:) + normrnd(0,0.3,[1,DIMENSION]);%5
%         x(j,:) = xbest + sqrt(c_cov)*normrnd(0,1,[1,1]) * v(j,:) + sqrt(1-c_cov)*normrnd(0,0.3,[1,DIMENSION]);%6
%         if rand < 0.1 %j == T
%             x(j,:) = xbest + 0.2*rand * v(j,:) + normrnd(0,0.3,[1,DIMENSION]);
%         end
        for k = 1:DIMENSION
            if abs(x(j,k)) > INIT_MAX
                xsign = sign(x(j,k));
                vecX = abs(reshape(x,1,POP_SIZE*DIMENSION));
                [row, col] = find(vecX < INIT_MAX);
                vecX = vecX(col);
                m = mean(vecX) + std(vecX) * randn * 10;
                x(j,k) = xsign * m; 
%                     v(j,k) = normrnd(0,1,[1,1]);
%                 x(j,k) = (rand(1,1)*(INIT_MAX-INIT_MIN)+INIT_MIN);
            end
        end
        
        cost1 = sum(p(j,:).*pp)/(norm(p(j,:))*norm(pp)); 
        pgtheta(j) = abs(acos(cost1));
        
        for i = 1:POP_SIZE
            for k = 1:POP_SIZE
                cost2 = sum(v(j,:).*v(k,:))/(norm(v(j,:))*norm(v(k,:))); 
                vvtheta(j,k) = abs(acos(cost2));
            end
        end

    end
    
    xbestold  = xbest;
    
    meanPGtheta(g) = mean(pgtheta);
    meanVVtheta(g) = mean(mean(vvtheta));
    
    xbestCoverg(g,:) = xbest;
    
    xCoverg(g,:) = mean(x,2);
    
    vCoverg(g,:) = mean(abs(v),2);
    
%     if IS_PLOT
%         figure(100);   
%         contour(xx,yy,z,50);
%         hold on;
%         plot(x(:,1),x(:,2),'k.','markersize',10);
%         hold on;
%         for i = 1:POP_SIZE
%             plot([xbest(1),x(i,1)], [xbest(2),x(i,2)],'g-','markersize',10);            
%             hold on
%         end
% %         plot(X(:,1),X(:,2),'k.','markersize',10);
%         hold on;
%         plot(xbest(1),xbest(2),'r.','markersize',20);
%         hold on;
%         plot([xbest(1),x(besti,1)], [xbest(2),x(besti,2)],'r-','markersize',10);
%         hold on;
%         plot([xbest(1),xbest(1)+pp(1)], [xbest(2),xbest(2)+pp(2)],'b-','markersize',10);
%         titlename = sprintf('fes=%d, bestf = %.2e',fitcount,bestf);
%         title(titlename);
%         hold off
%         pause(0.01)
%     end

%     if IS_PLOT && mod(g,50)==0
%         pause
%     end

    bestf = min(min(F));
%     fprintf('generation=%d , fes= %d, bestfit=%e  \n',g, fitcount, bestf);
    
end

fitness=bestf;

% figure(1)
% if ~IS_PLOT
%     BestFvalues = [];
%     for i = 1:numel(Fvalues)
%         BestFvalues(i) = min(Fvalues(1:i));
%     end
%     Bestfs = BestFvalues(1:10:end);
% else
%     BestFvalues = [];
%     for i = 1:numel(Fvalues)-10000
%         BestFvalues(i) = min(Fvalues(10000+1:10000+i));
%     end
%     Bestfs = BestFvalues(1:10:end);
% end
% loglog(Bestfs);
% xlabel('FES / D');
% ylabel('Fitness');
% title(sprintf('f%d',funcID));
% ylim([1e-2 1e4]);

% 
% figure(100)
% plot(meanPGtheta);
% xlabel('Generation');
% ylabel('theta');
% title(sprintf('f%d , vpbest vgbest',funcID));
% 
% figure(101)
% plot(meanVVtheta);
% xlabel('Generation');
% ylabel('theta');
% title(sprintf('f%d , v v',funcID));
% 
% figure(102)
% plot(xbestCoverg);
% xlabel('Generation');
% ylabel('value');
% title(sprintf('f%d , xbest dim',funcID));
% 
% figure(103)
% plot(xCoverg);
% xlabel('Generation');
% ylabel('value');
% title(sprintf('f%d , x dim',funcID));
% 
% figure(104)
% plot(vCoverg);
% xlabel('Generation');
% ylabel('value');
% title(sprintf('f%d , v dim',funcID));


