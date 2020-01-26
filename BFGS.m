function [x, ilaunch] = BFGS(FUN, DIM, lower, upper, ftarget, maxfunevals)
% minimizes FUN in D dimensions by independent restarts of fminunc (BFGS).
% ftarget and maxfunevals are additional external termination conditions.
% Search space is [-5, 5]^D
% set options, make sure we always terminate
options = optimset('fminunc');
options = optimset(options, 'LargeScale','off', 'DerivativeCheck','on'); % BFGS algorithm
options = optimset(options, 'MaxIter', inf, 'Tolfun', 1e-11, 'TolX', 0, ...
    'OutputFcn', @callback, 'Display', 'off');
% maxfunevals = min(1e4*DIM, maxfunevals);
% multistart such that ftarget is reached with reasonable prob.
for ilaunch = 1:10000 % relaunch optimizer up to 100 times
    options = optimset(options, 'MaxFunEvals', ...
        maxfunevals - feval(FUN, 'evaluations'));
    try
        x = fminunc(FUN, (upper-lower)*rand(DIM,1)+lower, options);
    catch
%         fprintf('error,  ilaunch %d \n',ilaunch);
    end
    if (feval(FUN, 'fbest') < ftarget || feval(FUN, 'evaluations') >= maxfunevals)
        break;
    end
end
% if ilaunch >= 100
%     disp('11111\n')
% end
    function stop = callback(x, optimValues, state)
        stop = false;
        if optimValues.fval < ftarget
%             stop = true;
        end
        if feval(FUN, 'evaluations') >= maxfunevals
            stop = true;
        end
    end
end % function