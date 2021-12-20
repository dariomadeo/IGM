function [S, U, I, R, D, x] = simode(theta, Idata, Rdata, Ddata, M, tBreaks, par_obj)
    
    ntBreaks = numel(tBreaks);
    
    
    S0 = getParameter(par_obj, ntBreaks, theta, 'S0', 1);
    x0 = getParameter(par_obj, ntBreaks, theta, 'x0', 1);    
    
    S = [];
    U = [];
    I = [];
    R = [];
    D = [];
    x = [];
   
    % Initial condition
    I0 = Idata(1);
    R0 = Rdata(1);
    D0 = Ddata(1);
    U0 = M - (S0+I0+R0+D0);    
    next_ic = [S0, U0, I0, R0, D0, x0]';
    
    odeopt = [];
%     odeopt = odeset('MaxStep', 1e-1);
%     'RelTol', 1e-5, 'AbsTol', 1e-8, 

    
    t1 = 0;
    for i=1:ntBreaks
        t2 = tBreaks(i);
        
        beta0       = getParameter(par_obj, ntBreaks, theta, 'beta0', i);
        gamma       = getParameter(par_obj, ntBreaks, theta, 'gamma', i);
        sigmaC0     = getParameter(par_obj, ntBreaks, theta, 'sigmaC0', i);
        sigmaN0     = getParameter(par_obj, ntBreaks, theta, 'sigmaN0', i);
        alpha       = getParameter(par_obj, ntBreaks, theta, 'alpha', i);
        tau         = getParameter(par_obj, ntBreaks, theta, 'tau', i); 
        omega       = getParameter(par_obj, ntBreaks, theta, 'omega', i);
        zeta        = getParameter(par_obj, ntBreaks, theta, 'zeta', i);
        e           = getParameter(par_obj, ntBreaks, theta, 'e', i);         
        a           = getParameter(par_obj, ntBreaks, theta, 'a', i); 
%         par_list = {'alpha', 'tau', 'omega', 'zeta', 'e', 'a'}; 
        
        % Simulate
        [~, sol] = ode45(@model, t1:t2, next_ic, odeopt, M, beta0, e, gamma, tau, omega, zeta, alpha, sigmaC0, sigmaN0, a);

        % Extract solution
        S = [S; sol(1:end-1,1)];
        U = [U; sol(1:end-1,2)];
        I = [I; sol(1:end-1,3)];
        R = [R; sol(1:end-1,4)];
        D = [D; sol(1:end-1,5)];
        x = [x; sol(1:end-1,6)];
        
        next_ic = sol(end,:)';
        
        t1 = t2;
    end
    
    S = [S; next_ic(1)];
    U = [U; next_ic(2)];
    I = [I; next_ic(3)];
    R = [R; next_ic(4)];
    D = [D; next_ic(5)];
    x = [x; next_ic(6)];    
end