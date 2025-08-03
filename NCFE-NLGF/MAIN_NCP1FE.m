%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonconforming P1 code for  u_t = div(\nabla u/ \sqrt{1 + |nabla u|^2) + f with  Neuman BC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long;
tic
%% Final times
T= 1; %1e-3;
itermax=1000;
tol=1e-12;

%% Epsilon
epsilon = 1;% minimal surface flow

%% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_4.mat';'mesh1_5.mat';'mesh1_6.mat'};
nbmeshes=size(meshes,1);

%% Initiations
% Errors
MAXL2error = zeros(nbmeshes, 1);
L1W11error = zeros(nbmeshes, 1);
L2H1error = zeros(nbmeshes, 1);

MAXL2norm = zeros(nbmeshes, 1);
L1W11norm = zeros(nbmeshes, 1);
L2H1norm = zeros(nbmeshes, 1);

% Order of convergence
ocMAXL2error = zeros(nbmeshes - 1, 1);
ocL1W11error = zeros(nbmeshes - 1, 1);
ocL2H1error = zeros(nbmeshes - 1, 1); 

ave_newton = zeros(nbmeshes, 1);

% To see the results printed in file
fid = fopen('results.txt','w');

%% Test case number
ucase = 1;

% Min relax
relaxmin = 1e-5;

str = sprintf('Final time is %d. \n', T);
forkprint(fid,str);
str = sprintf('Tolerance is %4.2e. \n',tol);
forkprint(fid,str);
str = sprintf('Epsilon is %4.2e. \n', epsilon);
forkprint(fid,str);

%% Loop over each mesh in the sequence
for imesh=1:nbmeshes
 
    %% Load mesh here!
    loadmesh=strcat('load ../matlab_meshes/', meshes{imesh});
    disp(loadmesh);
    eval(loadmesh);
    % disp('mesh loaded');

    %% Compute real centers of mass, mesh size, and midpoint of all edges
    cg=gravity_centers(ncell, cell_v, vertex, area); %centers of mass
    h(imesh)=max(abs(diam)); %mesh size
    mpe=midpoints_edges(ncell,nedge,cell_e,cell_v,vertex); %midpoint of all edges

    %% Time steps
    Ndt(imesh) = ceil(T/h(imesh)); %k = O(h)
    % Ndt(imesh) = ceil(T/h(imesh)^2); %k = O(h^2)
    dt = T/Ndt(imesh);
    ave_newton(imesh) = 0;

    %% Initial condition
    U_pre = test_cases(0,mpe, ucase)';

    %% Assemble mass matrices
    M = assemble_mass_system(area, ncell, cell_e, nedge);

    % Print the scheme solution and exact solution at initial condition to view in Paraview
    if imesh == nbmeshes
    write_solution_vtk_ncP1(U_pre, strcat('VTKout/ncp1_solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    write_solution_vtk_ncP1(U_pre, strcat('VTKout/exact_solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    end

    %% Error norms initiation
    L2error = zeros(Ndt(imesh), 1);
    W11error = zeros(Ndt(imesh), 1);
    H1error = zeros(Ndt(imesh), 1);

    L2norm = zeros(Ndt(imesh),1);
    W11norm = zeros(Ndt(imesh),1);
    H1norm = zeros(Ndt(imesh),1);

    %% Time stepping starts here!
    ITER = 0;
    num_updates=0;
    Res = 0;
    for idt = 1 : Ndt(imesh)
        b = assemble_source(cell_e, ncell, nedge, area, cg, idt * dt, epsilon, ucase);

        [U, num_updates, iter, res] = compute_staionary_system(cell_v, cell_e, ncell, nedge, vertex, area, mpe, dt, idt, imesh, epsilon, itermax, tol, relaxmin, num_updates, M, b, U_pre);
       
        %% Newton data
        ave_newton(imesh) = ave_newton(imesh) + iter;
        ITER = ITER + iter;
        Res = Res + abs(res);

        %% Errors and norm of the exact solution
        [L2error(idt), W11error(idt), H1error(idt)] = compute_norms(U - test_cases(idt * dt, mpe, ucase)', area, mpe, cell_e);
        [L2norm(idt), W11norm(idt), H1norm(idt)] = compute_norms(test_cases(idt * dt, mpe, ucase)', area, mpe, cell_e);

        str = sprintf('Solution computed, iter=%d, res=%4.2e \n', iter, res);
        % forkprint(fid,str);

         U_pre = U;
        % Create files for visualising the approximate and exact solution
        if imesh == nbmeshes
        write_solution_vtk_ncP1(U,strcat('VTKout/ncp1_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the scheme solution at current time in Paraview
        write_solution_vtk_ncP1(test_cases(idt*dt, mpe, ucase)',strcat('VTKout/exact_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the exact solution at current time in Paraview
        write_solution_vtk_ncP1(zeros(size(U)),'VTKout/grid',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
        end
    end

    ITER=ceil(ITER/Ndt(imesh));
    Res=(Res/Ndt(imesh));
    num_updates=ceil(num_updates/Ndt(imesh));
    ave_newton(imesh) = ave_newton(imesh)/Ndt(imesh);

    str = sprintf('Mesh %i, num_time_steps=%d, Avg_iter=%d with Avg_relax=%4.2e, Avg_Newt=%d, Avg_residue=%4.2e\n',imesh,Ndt(imesh),ITER,num_updates, ave_newton(imesh),Res);
    forkprint(fid,str);
    clear ITER Res num_updates ave_newton

    %% Time norms: L infinity, L1 norm, and L2 norm
    MAXL2error(imesh) = max(L2error);
    L1W11error(imesh) = dt * sum(W11error);
    L2H1error(imesh) = sqrt(dt * sum(H1error.^2));

    MAXL2norm(imesh) = max(L2norm);
    L1W11norm(imesh) = dt * sum(W11norm);
    L2H1norm(imesh) = sqrt(dt * sum(H1norm.^2));

    %% Relative errors
    MAXL2error(imesh) = MAXL2error(imesh)/MAXL2norm(imesh);
    L1W11error(imesh) = L1W11error(imesh)/L1W11norm(imesh);
    L2H1error(imesh) = L2H1error(imesh)/L2H1norm(imesh);

    str = sprintf('Mesh %i. MAXL2error=%4.2e, L1W11error:%4.2e, L2H1error:%4.2e \n',imesh, MAXL2error(imesh), L1W11error(imesh), L2H1error(imesh));
    forkprint(fid,str);

    Time(imesh)=toc;
end
str = sprintf(['\nElapsed time is ' num2str(toc) 'seconds']);
forkprint(fid,str);

%% Orders of convergence
ocL2error=zeros(nbmeshes-1,1);
ocH10error=zeros(nbmeshes-1,1);
for imesh = 1:nbmeshes-1
    ocMAXL2error(imesh)=log(MAXL2error(imesh)/MAXL2error(imesh+1)) / log(h(imesh)/h(imesh+1));
    ocL1W11error(imesh)=log(L1W11error(imesh)/L1W11error(imesh+1)) / log(h(imesh)/h(imesh+1));
    ocL2H1error(imesh)=log(L2H1error(imesh)/L2H1error(imesh+1)) / log(h(imesh)/h(imesh+1));
end

str = sprintf('\nMaxL2 error and convergence rate:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',MAXL2error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',MAXL2error(imesh),ocMAXL2error(imesh-1));
        forkprint(fid,str);
    end
end

str = sprintf('\nL1W11 error and convergence rate:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',L1W11error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',L1W11error(imesh), ocL1W11error(imesh-1));
        forkprint(fid,str);
    end
end

str = sprintf('\nL2H1 error and convergence rate:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',L2H1error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',L2H1error(imesh), ocL2H1error(imesh-1));
        forkprint(fid,str);
    end
end

% Write data file
fid = fopen('data_rates.dat','w');
fprintf(fid,'meshsize timestep MaxL2error L2W11error Time\n');
for i=1:nbmeshes
    fprintf(fid,'%f %f %f %f %f %f\n', h(i), T/Ndt(i), MAXL2error(i), L1W11error(i), L2H1error(i) , Time(i));
end;

%% Plot of the MAXL2 error of u and L1W11 error Gradu with respect of h in a logarithmic scale
subplot(1,3,1);
loglog(h,MAXL2error,'*-');
axis([5e-3 5e-1 1e-3 7e-1]);
hold on
% loglog(h,h); % I compute the vector in order to have a slope equal to 1
ax = axis;
xmin = ax(1); xmax = ax(2);
ymin = ax(3); ymax = ax(4);
% Triangle scale (slope 1: same factor in log-scale)
r = 3;  % scaling ratio (how large the triangle is)
% Choose a point at bottom-right of the plot
x0 = 9e-2;     % starting x (closer to the right)
y0 = 6e-3;     % starting y (closer to bottom)
% Triangle vertices (flipped horizontally)
x1 = x0 * r;         % move right (flipped direction)
y1 = y0;
x2 = x0 * r;         
y2 = y0 * r;         % move up to keep slope 1
% Draw the triangle
plot([x0, x1], [y0, y1], 'k', 'LineWidth', 1,'HandleVisibility','off');         % horizontal
plot([x1, x2], [y1, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % vertical
plot([x0, x2], [y0, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % diagonal (hypotenuse)
%Label horizontal edge with "1"
text(sqrt(x0 * x1), y0 * 0.8, '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'center');
% Label vertical edge with "1"
text(x1 * 1.1, sqrt(y1 * y2), '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

title('Errors in $L^{\infty}L^2$ norm on the functions','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
% ylabel('Error','Interpreter','latex')
legend('NCP1FEM','Location','southeast');
% grid on
subplot(1,3,2);
loglog(h,L1W11error,'*-');
axis([5e-3 5e-1 1e-3 7e-1]);
hold on
% loglog(h,h);
% Draw the triangle
plot([x0, x1], [y0, y1], 'k', 'LineWidth', 1,'HandleVisibility','off');         % horizontal
plot([x1, x2], [y1, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % vertical
plot([x0, x2], [y0, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % diagonal (hypotenuse)
%Label horizontal edge with "1"
text(sqrt(x0 * x1), y0 * 0.8, '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'center');
% Label vertical edge with "1"
text(x1 * 1.1, sqrt(y1 * y2), '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

title('Errors in $L^1L^{1}$ norm on the gradients','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
% ylabel('Error','Interpreter','latex')
legend('NCP1FEM','Location','southeast')
% grid on
subplot(1,3,3);
loglog(h,L2H1error,'*-');
axis([5e-3 5e-1 1e-3 7e-1]);
hold on
% loglog(h,h);
% Draw the triangle
plot([x0, x1], [y0, y1], 'k', 'LineWidth', 1,'HandleVisibility','off');         % horizontal
plot([x1, x2], [y1, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % vertical
plot([x0, x2], [y0, y2], 'k', 'LineWidth', 1,'HandleVisibility','off');         % diagonal (hypotenuse)
%Label horizontal edge with "1"
text(sqrt(x0 * x1), y0 * 0.8, '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'center');
% Label vertical edge with "1"
text(x1 * 1.1, sqrt(y1 * y2), '1', ...
     'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

title('Errors in $L^2L^{2}$ norm on the gradients','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
% ylabel('Error','Interpreter','latex')
legend('NCP1FEM','Location','southeast')
% hold on
% grid on


