function [fitness,generation]=fips_uring(FUN,POP_SIZE,DIMENSION,INIT_MIN,INIT_MAX,ftarget,MAX_FES)
PLOT_INTERVAL=64*5; %3200;
MAX_GENERATION = floor(MAX_FES/POP_SIZE);
x=rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;
v=(rand(POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;
pf=ones(1,POP_SIZE)*1.0E100;
ppf=ones(1,POP_SIZE)*1.0E100;
pppf=1.0E100;

p=zeros(POP_SIZE,DIMENSION);
pp=zeros(POP_SIZE,DIMENSION);
fitcount=0;
generation=zeros(floor(MAX_FES/PLOT_INTERVAL),1);


for g=1:MAX_GENERATION
    if(fitcount>MAX_FES) || feval(FUN, 'fbest') < ftarget; break;end
    f = feval(FUN, x')';
    fitcount=fitcount+POP_SIZE;
    
%     for j=1:POP_SIZE
%         f(j)=fitnessf(x(j,:),FUNNUM);
%         fitcount=fitcount+1;
%         %»æÍ¼
%         
%         if(fitcount==1||mod(fitcount,PLOT_INTERVAL)==0)
%             if(fitcount~=1)
%                 generation(fitcount/PLOT_INTERVAL+1)=pppf;
%                 best=pppf;
%             else
%                 generation(1)=f(1,1);
%                 best=f(1,1);
%             end;
%             if(IS_PLOT==1)
%                 hold on;
%                 plot(fitcount,best,'.','MarkerSize',5);
%                 title(sprintf('Fitness=%e',best));
%                 xlim([0 MAX_FES]);
%                 pause(0.01);
%             end
%         end;
%     end
    
    
    for j=1:POP_SIZE
        if(f(j)<pf(j))
            p(j,:)=x(j,:);
            pf(j)=f(j);
        end
    end
    
    for j=1:POP_SIZE
        for k=-1:1
            if(pf(mod((j-1)+k,POP_SIZE)+1)<ppf(j))
                pp(j,:)=p(mod((j-1)+k,POP_SIZE)+1,:);
                ppf(j)=pf(mod((j-1)+k,POP_SIZE)+1);
            end
        end
    end
    
    for j=1:POP_SIZE
        if(ppf(j)<pppf)
            pppf=ppf(j);
        end
    end
    
    for j=1:POP_SIZE
        for m=1:DIMENSION
            fi1=rand*4.1/2;
            fi2=rand*4.1/2;
            k1=mod((j-1)-1,POP_SIZE)+1;
            k2=mod((j-1)+1,POP_SIZE)+1;
            Pm=(pf(k1)*fi1*p(k1,m)+pf(k2)*fi2*p(k2,m))/(pf(k1)*fi1+pf(k2)*fi2);
            v(j,m)= 0.7298*(v(j,m)+(fi1+fi2)*(Pm-x(j,m)));
            x(j,m)= v(j,m)+x(j,m);
        end
    end
end
fitness=pppf;



