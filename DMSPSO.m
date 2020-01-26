function [fitness,generation]=DMSPSO(FUN,POP_SIZE,PPOP_SIZE,DIMENSION,INIT_MIN,INIT_MAX,ftarget,MAX_FES)
PLOT_INTERVAL=64*5;%3200;%
MAX_GENERATION = floor(MAX_FES/(POP_SIZE*PPOP_SIZE/2));
iwt=0.9-(1:MAX_GENERATION)*(0.7/MAX_GENERATION);%inertia weight

x=rand(PPOP_SIZE*POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN;
v=(rand(PPOP_SIZE*POP_SIZE,DIMENSION)*(INIT_MAX-INIT_MIN)+INIT_MIN)*0.2;
pf=ones(1,PPOP_SIZE*POP_SIZE)*1.0E100;
ppf=ones(1,PPOP_SIZE)*1.0E100;

pppf=1.0e100;
p=zeros(PPOP_SIZE*POP_SIZE,DIMENSION);
pp=zeros(PPOP_SIZE,DIMENSION);
ppp=zeros(1,DIMENSION);
fitcount=0;
generation=zeros(floor(MAX_FES/PLOT_INTERVAL),1);


for g=1:floor(0.9*MAX_GENERATION)
    if(fitcount>MAX_FES) || feval(FUN, 'fbest') < ftarget; break;end
    if(mod(g,5)==0)
        for j=1:POP_SIZE*PPOP_SIZE
            order(j)=j;
        end;
        for j=0:POP_SIZE*PPOP_SIZE-1
            r=j+mod(ceil(rand*65535),(POP_SIZE*PPOP_SIZE-j));
            temp=order(j+1);
            order(j+1)=order(r+1);
            order(r+1)=temp;
        end;
        
        for i=1:PPOP_SIZE*POP_SIZE
            x1(i,:)=x(order(i),:);
            v1(i,:)=v(order(i),:);
            p1(i,:)=p(order(i),:);
            pf1(i)=pf(order(i));
        end;
        
        for i=1:PPOP_SIZE*POP_SIZE
            x(i,:)=x1(i,:);
            v(i,:)=v1(i,:);
            p(i,:)=p1(i,:);
            pf(i)=pf1(i);
        end
        
        for i=1:PPOP_SIZE
            ppf(i)=1.0e100;
        end
        pp=zeros(PPOP_SIZE,DIMENSION);
    end;
    
    f = feval(FUN, x')';
    fitcount=fitcount+PPOP_SIZE*POP_SIZE;
%     for i=1:PPOP_SIZE*POP_SIZE
%         temp(1,1:DIMENSION)=x(i,1:DIMENSION);
%         f(i) = feval(FUN, temp')';
%         fitcount=fitcount+1;
% %         f(i)=fitnessf(temp,FUNNUM);
% %         fitcount=fitcount+1;
% %         %»æÍ¼
% %         if(fitcount==1||mod(fitcount,PLOT_INTERVAL)==0)
% %             if(fitcount~=1)
% %                 generation(fitcount/PLOT_INTERVAL+1)=pppf;
% %                 best=pppf;
% %             else
% %                 generation(1)=f(1,1);
% %                 best=f(1,1);
% %             end;
% %             if(IS_PLOT==1)
% %                 hold on;
% %                 plot(fitcount,best,'.','MarkerSize',5);
% %                 title(sprintf('Fitness=%e',best));
% %                 xlim([0 MAX_FES]);
% %                 pause(0.01);
% %             end
% %         end;
%     end
    
    for i=1:PPOP_SIZE*POP_SIZE
        if(f(i)<pf(i))
            p(i,:)=x(i,:);
            pf(i)=f(i);
        end
    end
    
    for i=1:PPOP_SIZE*POP_SIZE
        j=floor((i-1)/POP_SIZE)+1;
        if(pf(i)<ppf(j))
            pp(j,:)=p(i,:);
            ppf(j)=pf(i);
        end
    end
    
    for i=1:PPOP_SIZE
        if(ppf(i)<pppf)
            ppp(:)=pp(i,:);
            pppf=ppf(i);
        end
    end
    
    for i=1:PPOP_SIZE*POP_SIZE
        for k=1:DIMENSION
            j=floor((i-1)/POP_SIZE)+1;
            v(i,k)= iwt(g)*v(i,k)+2*rand*(p(i,k)-x(i,k))+2*rand*(pp(j,k)-x(i,k));
            if(v(i,k)>(INIT_MAX-INIT_MIN)*0.2)
                v(i,k)=(INIT_MAX-INIT_MIN)*0.2;
            else
                if(v(i,k)<(INIT_MAX-INIT_MIN)*-0.2)
                    v(i,k)=(INIT_MAX-INIT_MIN)*-0.2;
                end;
            end;
            x(i,k)= v(i,k)+x(i,k);
        end
    end
end
clear ppf pp
ppf=pppf;
pp=ppp;
clear pppf ppp
for g=floor(0.9*MAX_GENERATION):MAX_GENERATION
    if(fitcount>MAX_FES) || feval(FUN, 'fbest') < ftarget; break;end
    f = feval(FUN, x')';
    fitcount=fitcount+POP_SIZE*PPOP_SIZE;
%     for j=1:POP_SIZE*PPOP_SIZE
%         
%         f(j)=fitnessf(x(j,:),FUNNUM);
%         fitcount=fitcount+1;
%         %»æÍ¼
%         if(fitcount==1||mod(fitcount,PLOT_INTERVAL)==0)
%             if(fitcount~=1)
%                 generation(fitcount/PLOT_INTERVAL+1)=ppf;
%                 best=ppf;
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
    
    
    for j=1:POP_SIZE*PPOP_SIZE
        if(f(j)<pf(j))
            p(j,:)=x(j,:);
            pf(j)=f(j);
        end
    end
    
    
    for j=1:POP_SIZE*PPOP_SIZE
        if(pf(j)<ppf)
            pp(:)=p(j,:);
            ppf=pf(j);
        end
    end
    
    
    for j=1:POP_SIZE*PPOP_SIZE
        for k=1:DIMENSION
            v(j,k)= iwt(g)*v(j,k)+2*rand*(p(j,k)-x(j,k))+2*rand*(pp(k)-x(j,k));
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
end;
fitness=ppf;

