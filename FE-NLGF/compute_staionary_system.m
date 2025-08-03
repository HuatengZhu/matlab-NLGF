function [U, num_updates, iter, res] = compute_staionary_system(cell_v, ncell, nvert, vertex, dt, epsilon, itermax, tol, relaxmin, num_updates, M, b, U_pre);

rhs = M * U_pre + dt * b;
iter = 0;
res = 1;

%% Newton
X_pre = U_pre; %test_cases(idt*dt,vertex, ucase)'; % initial guess of the solution
relax = 0.1;
res_pre = 1.;
while ( iter < itermax && res > tol && relax > relaxmin )
    [J, F] = assemble_jacobian_system(cell_v,ncell,nvert,vertex, X_pre, epsilon);
   
    Aglob = M + dt .* J;
    RHS = rhs - ( M + dt .* F)* X_pre ;

    delta_X = Aglob\RHS;

    X = X_pre +  relax * delta_X;
    iter = iter + 1;

    %% Residues
    res = norm(RHS, Inf);% residual of nonlinear eq at previous step: ||F(X_pre) - rhs||
    % fprintf("res %e, relax %e, idt/Ndt %d/%d\n", res, relax, idt, Ndt(imesh));

    %% If the residual is too large, we don't update and we reduce relax
    if (res > 1.2*res_pre)
        iter = iter-1;
        relax = relax/5;
        num_updates=num_updates+1;
    else
        X_pre = X;
        res_pre = res;
        relax = min(1,relax*1.2);
    end;

end % end nonlinear iterations

if (iter==itermax)
    res
    iter
    error('no convergence')
end;
% Newton loop end: X outputed here is the solution at the current time step idt*dt

U = X; % current time solution is retireing