function dy = model(t, y, M, beta0, e, gamma, tau, omega, zeta, alpha, sigmaC0, sigmaN0, a)

    %% RHS of the ODE system
    
    S = y(1);
    U = y(2);
    I = y(3);
    R = y(4);
    D = y(5);
    x = y(6);

    beta = beta0*(1-e*x);
    
    dS =   -beta*S*U/M                                 +alpha*R;
    dU =   +beta*S*U/M -gamma*U -tau*U;
    dI =                        +tau*U -omega*I -zeta*I;                                                                                   
    dR =               +gamma*U        +omega*I        -alpha*R;
    dD =                                        +zeta*I;                                                                   
    
    dx = x*(1-x)*((sigmaC0+sigmaN0)*x - sigmaN0)*(I+D-a)/M;

    dy = [dS; dU; dI; dR; dD; dx];

end
