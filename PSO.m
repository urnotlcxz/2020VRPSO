function PSO(FUN, DIM, lower, upper, ftarget, maxfunevals)
% Set algorithm parameters
popsize = 40;
c1 = 1.4944;%2;
c2 = 1.4944;%2;
w = 0.792;
xbound = upper;
vbound = upper;
% Allocate memory and initialize
xmin = -xbound * ones(1,DIM);
xmax = xbound * ones(1,DIM);
vmin = -vbound * ones(1,DIM);
vmax = vbound * ones(1,DIM);
x = 2 * xbound * rand(popsize,DIM) - xbound;
v = 2 * vbound * rand(popsize,DIM) - vbound;
pbest = x;
% update pbest and gbest
cost_p = feval(FUN, pbest');
[cost,index] = min(cost_p);
gbest = pbest(index,:);
maxfunevals = min(1e5 * DIM, maxfunevals);
%maxfunevals = min(10000, maxfunevals);
maxiterations = ceil(maxfunevals/popsize);
for iter = 2 : maxiterations
    % Update inertia weight
    % w = 0.9 - 0.8*(iter-2)/(maxiterations-2);
    % Update velocity
    v = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).*(repmat(gbest,popsize,1)-x);
    % Clamp veloctiy
    s = v < repmat(vmin,popsize,1);
    v = (1-s).*v + s.*repmat(vmin,popsize,1);
    b = v > repmat(vmax,popsize,1);
    v = (1-b).*v + b.*repmat(vmax,popsize,1);
    % Update position
    x = x + v;
    % Clamp position - Absorbing boundaries
    % Set x to the boundary
    s = x < repmat(xmin,popsize,1);
    x = (1-s).*x + s.*repmat(xmin,popsize,1);
    b = x > repmat(xmax,popsize,1);
    x = (1-b).*x + b.*repmat(xmax,popsize,1);
    % Clamp position - Absorbing boundaries
    % Set v to zero
    b = s | b;
    v = (1-b).*v + b.*zeros(popsize,DIM);
    % Update pbest and gbest if necessary
    cost_x = feval(FUN, x'); %%%%%%%%   ”¶÷µ∆¿π¿
    s = cost_x<cost_p;
    cost_p = (1-s).*cost_p + s.*cost_x;
    s = repmat(s',1,DIM);
    pbest = (1-s).*pbest + s.*x;
    [cost,index] = min(cost_p);
    gbest = pbest(index,:);
    % Exit if target is reached
    if feval(FUN, 'fbest') < ftarget
        break;
    end
end