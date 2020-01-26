function [x, ilaunch, f] = NELDER(FUN, DIM, lower, upper, FTarget, maxfunevals)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed

options = optimset('MaxFunEvals', min(1e9*DIM, maxfunevals), ...
    'MaxIter', 2e9*DIM, ... % overwritten later
    'Tolfun', 1e-11, ...
    'TolX', 1e-11, ...
    'OutputFcn', @callback, ...
    'Display', 'off');
ilocal = 0;
% multistart such that ftarget is reached with reasonable prob.
for ilaunch = 1:1e5; % relaunch optimizer up to 1e5 times，最少重启10次
    % set initial conditions
    ilocal = ilocal + 1;
    if ilocal == 1 % (re-)start from scratch
        xstart = (upper - lower) * rand(DIM, 1) + lower; % random start solution
        usual_delta = 2 * ones(DIM, 1); % %%%%%%%% 修改过
        options = optimset(options, 'MaxIter', floor(200*sqrt(DIM)));
    else % refining restart run
        xstart = x; % try to improve found solution
        usual_delta = 10 * myrange(v,2); % a bit of regularization in given coordinate sys
        if rand(1,1) < 0.2 % a bit of desperatation
            usual_delta = usual_delta + (1/ilocal) * (0.1/ilocal).^rand(DIM,1);
        end
        if rand(1,1) < 0.1 * (ilocal-10)/sqrt(DIM) % final run
            options = optimset(options, 'MaxIter', 500*DIM); % long run
            ilocal = 0; % real restart after this run
        end
    end
    % try fminsearch from Matlab, modified to take usual_delta as arg
    [x,f,e,o,v] = fminsearch_mod(FUN, xstart, usual_delta, options); % MaxIter 或者 MaxFunEvals 只要有一个限制，就会停止
%     fminsearch
    % FTarget 这里就等于 Fopt + DeltaFtarget
    % DeltaFtarget 就是 precision 精确度
    % Fopt 就是CEC以前的fbias 
    % 所以函数的最优值就是以某个精度接近偏置
    % 函数值减去函数最优值 就得到了真正的最优值，接近0的
    % 这里的f等于FTrue + Fopt
    % feval(FUN, 'fbest') 是f里最好的那个
%     fprintf( 'ilocal = %d, \t fes = %d, \t ftrue-preci = %.4e, \t min-delta = %.4e, \t mxd/mnd = %.4e \n', ...
%             ilocal, feval(FUN, 'evaluations'), f-FTarget, min(usual_delta), max(usual_delta)/min(usual_delta) );
%     feval(FUN, 'fbest')
%     feval(FUN, 'evaluations')
    if feval(FUN, 'fbest') < FTarget || feval(FUN, 'evaluations') >= maxfunevals % 达到 maxfunevals 或精度条件 FTrue < DeltaFtarget 就必须停止，不再重启
%         fprintf('bestX: %.2f, %.2f, %.2f, %.2f, %.2f\n', x(1), x(2), x(3), x(4), x(5));
        break;
    end
    % if useful, modify more options here for next launch
end
    function stop = callback(x, optimValues, state)
        stop = false;
        if optimValues.fval < FTarget || feval(FUN, 'evaluations') >= maxfunevals
            stop = true;
        end
    end
end

function y = myrange(x,dim)
%RANGE  Sample range.
%   Y = RANGE(X) returns the range of the values in X.  For a vector input,
%   Y is the difference between the maximum and minimum values.  For a
%   matrix input, Y is a vector containing the range for each column.  For
%   N-D arrays, RANGE operates along the first non-singleton dimension.
%
%   RANGE treats NaNs as missing values, and ignores them.
%
%   Y = RANGE(X,DIM) operates along the dimension DIM.
%
%   See also BOUNDS, MIN, MAX, IQR, MAD, STD.

%   Copyright 1993-2016 The MathWorks, Inc.


if nargin < 2
    y = max(x) - min(x);
else
    y = max(x,[],dim) - min(x,[],dim);
end

end


