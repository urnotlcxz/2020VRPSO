function [fitness,generation]=LPSO(FUN,POP_SIZE,DIMENSION,INIT_MIN,INIT_MAX,ftarget,MAX_FES)
PLOT_INTERVAL=64*5; %3200;%
MAX_GENERATION = floor(MAX_FES/POP_SIZE);
iwt=0.9-(1:MAX_GENERATION)*(0.5/MAX_GENERATION);%inertia weight
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
        for k=1:DIMENSION
            v(j,k)= iwt(g)*v(j,k)+2*rand*(p(j,k)-x(j,k))+2*rand*(pp(j,k)-x(j,k));
            if(v(j,k)>(INIT_MAX-INIT_MIN)*0.2)
                v(j,k)=(INIT_MAX-INIT_MIN)*0.2;
            else
                if(v(j,k)<(INIT_MAX-INIT_MIN)*-0.2)
                    v(j,k)=(INIT_MAX-INIT_MIN)*-0.2;
                end;
            end;
            x(j,k)= v(j,k)+x(j,k);
        end
    end
    
end
fitness=pppf;



