clc;
clear ;
close all;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)*randi(100,1)));
global nname
global FUN_NUM FITCOUNT Fvalues MAX_FES ftarget
global xx yy z IS_PLOT

addpath(genpath('RNN/BBOB09'));  % should point to fgeneric.m etc.
rmpath(genpath('RNN\BBOB09\MCS'));
dbstop if error

DIMENSION = 10;
opt.repeat = 1;
funcID = 13;
ftarget = fgeneric('initialize',...
    funcID,... % funID
    2 ,... % dim
    1,... % inst
    1, 'F:\\STAT', opt);
ft = 1e-8;
nname=sprintf('D:\\桌面\\ANN优化\\递归新\\DATA\\nn\\111122.mat');
FITCOUNT = 0;
Fvalues = [];
FUN_NUM = 11;

POP_SIZE = 40;

INIT_MIN = -5;INIT_MAX = 5;
MAX_FES=5000*DIMENSION;
MAX_GENERATION=floor(MAX_FES/POP_SIZE);
PLOT_INTERVAL=POP_SIZE*5;
iwt=0.9-(1:MAX_GENERATION)*(0.5/MAX_GENERATION);%inertia weight

% c_sigma = 0.5-(1:MAX_GENERATION)*(0.01/MAX_GENERATION);
c_sigma = 0.1-(1:MAX_GENERATION)*(0.1/MAX_GENERATION);
v_sigma = 1.5-(1:MAX_GENERATION)*(1.5/MAX_GENERATION);
s_sigma = 15-(1:MAX_GENERATION)*(10/MAX_GENERATION);

xinit = rand(1,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;
% xinit = ones(1,DIMENSION)*4 + normrnd(0,0.3,[1,DIMENSION]);
% x = repmat(xinit,POP_SIZE,1) + 1 * normrnd(0,1, [POP_SIZE,DIMENSION]);

x=rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;

v=(rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;

% initialize(DIMENSION,50,INIT_MAX,INIT_MIN);

IS_PLOT = 1;

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

pf=ones(1,POP_SIZE)*1.0E100;
ppf=1.0E100;

p=zeros(POP_SIZE,DIMENSION);
pp=zeros(1,DIMENSION);
fitcount=0;
generation=zeros(floor(MAX_FES/PLOT_INTERVAL),1);
bestf = inf;

X=[];
F=[];
meanTheta = [];

sigma = 1;

for g=1:MAX_GENERATION
    if(fitcount>MAX_FES || bestf < ft); break;end;
    f=fitnessf(x);
    fitcount=fitcount+POP_SIZE;
    
    X=[X;x];
    F=[F;f];
    bestflag = 1;
    for j=1:POP_SIZE
        if(f(j)<pf(j))
            p(j,:) = v(j,:);
            pf(j)=f(j);
            if(pf(j)<ppf)
                pp(:)=p(j,:);
                ppf=pf(j);
                
                xbest = x(j, :);
                bestflag = 0;
                besti = j;
            end
        end
    end
    
    if g == 1
        xbestold = xbest;
    end
    
    if bestflag == 0
        sigma = sigma * 1.3;
    else
        sigma = sigma * 0.95;
    end
    
    I = eye(DIMENSION);
    Mu0 = zeros(1,DIMENSION);
    N0 = mvnrnd(Mu0, I, POP_SIZE);
    
    cos_theta = [];
    theta = [];
    for k = 1:POP_SIZE
        cos_theta(k) = sum(N0(k,:).*pp)/(norm(N0(k,:))*norm(pp));
        theta(k) = acos(cos_theta(k));
    end
    
    alpha = 1;
    sig = 1;
    for k = 1:POP_SIZE
        if theta(k) <= pi/2
            scale(k,1) = alpha * exp(-theta(k)^2/(2*sig));
        elseif theta(k) > pi/2
            scale(k,1) = alpha * exp(-(pi-theta(k))^2/(2*sig));
        else
            error('ddd');
        end
    end
    
    Nc = [];
    for k = 1:POP_SIZE
        if rand < 0.7
            flg = 1;
        else            
            flg = 0; 
        end
        Nc(k,:) = N0(k,:) * norm(pp)/norm(N0(k,:)) * scale(k) + flg * randn * pp;
    end
    
%     Nc = N0 .* repmat(scale,1,DIMENSION);
    
    x = repmat(xbest,POP_SIZE,1) + sigma * Nc;
    
        
    if IS_PLOT
        figure(100);
        contour(xx,yy,z,50);
        hold on;
        plot(x(:,1),x(:,2),'k.','markersize',10);
        hold on;
        for i = 1:POP_SIZE
            plot([xbest(1),x(i,1)], [xbest(2),x(i,2)],'g-','markersize',10);
            hold on
        end
        %         plot(X(:,1),X(:,2),'k.','markersize',10);
        hold on;
        plot(xbest(1),xbest(2),'r.','markersize',20);
        hold on;
        plot([xbest(1),x(besti,1)], [xbest(2),x(besti,2)],'r-','markersize',10);
        hold on;
        plot([xbest(1),xbest(1)+pp(1)], [xbest(2),xbest(2)+pp(2)],'b-','markersize',10);
        titlename = sprintf('fes=%d, bestf = %.2e',fitcount,bestf);
        title(titlename);
        hold off
        pause(0.1)
    end
    
    %     if mod(g,20)==0
    %         pause
    %     end
    
    bestf = ppf;
    fprintf('generation=%d , fes= %d, bestfit=%e  \n',g, FITCOUNT, bestf);
    
end

fitness=ppf;

figure(1)
if ~IS_PLOT
    BestFvalues = [];
    for i = 1:numel(Fvalues)
        BestFvalues(i) = min(Fvalues(1:i));
    end
    Bestfs = BestFvalues(1:DIMENSION:end);
else
    BestFvalues = [];
    for i = 1:numel(Fvalues)-10000
        BestFvalues(i) = min(Fvalues(10000+1:10000+i));
    end
    Bestfs = BestFvalues(1:DIMENSION:end);
end
loglog(Bestfs);
xlabel('FES / D');
ylabel('Fitness');
title(sprintf('f%d',funcID));
% ylim([1e-2 1e4]);

