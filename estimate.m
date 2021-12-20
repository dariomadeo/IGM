function thetaOptim = estimate(W, Tw, Idata_est, Rdata_est, Ddata_est, lb, ub, theta0, par_obj, M)


    NT = numel(Idata_est);
    tFinal = NT-1;
    tSpan = 0:tFinal;
    tBreaks = (Tw:Tw:Tw*W)-1; 
    
    % Estimation ...
    fun_handler = @(theta, tSpan)foroptim(theta, tSpan, Idata_est, Rdata_est, Ddata_est, M, tBreaks, par_obj);
    
    optopt = optimoptions('lsqcurvefit', ...
                          'Display', 'iter', ... 
                          'Algorithm', 'trust-region-reflective', ...
                          'MaxFunctionEvaluations',20000, ...
                          'FunctionTolerance', 1e-14, ...
                          'OptimalityTolerance', 1e-14, ...
                          'StepTolerance', 1e-14);

    [thetaOptim, resnorm, ~, exitflag, output] = lsqcurvefit(fun_handler, theta0, tSpan, [Idata_est; Rdata_est; Ddata_est]./[Idata_est; Rdata_est; Ddata_est], lb, ub, optopt);



end

function out = foroptim(theta, tSpan, Idata_win, Rdata_win, Ddata_win, M, tBreaks, par_obj)

    [S, U, I, R, D, x] = simode(theta, Idata_win, Rdata_win, Ddata_win, M, tBreaks, par_obj); %simode(theta, Idata_win, Rdata_win, Ddata_win, M, beta0, gamma, tBreaks);
    
    out = [I; R; D]./[Idata_win; Rdata_win; Ddata_win];

end


