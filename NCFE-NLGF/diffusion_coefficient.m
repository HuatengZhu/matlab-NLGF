function [A, grad_u] = diffusion_coefficient(u,i,epsilon,mpedges,cell_e)


%% Barycentric coordinates of cell i
midpoints = [mpedges(cell_e{i}(1),:)' mpedges(cell_e{i}(2),:)' mpedges(cell_e{i}(3),:)'];    
L1 = [ones(1,3); midpoints]'\ [1;0;0];
	% L1 are the coordinates of the 1st barycentric coordinate lambda_1 (associated to
	%		the 1st vertex). lambda_1 is a function of x,y, given by
	%		lambda_1(x,y)= alpha0 + alpha1 x + alpha2 y
	%		if L1=[alpha0 alpha1 alpha2]
L2 = [ones(1,3); midpoints]'\[0;1;0];
L3 = [ones(1,3); midpoints]'\[0;0;1];

% gradient at the gravity centre
grad_u= u(cell_e{i}(1))*L1(2:3) + u(cell_e{i}(2))*L2(2:3) + u(cell_e{i}(3))*L3(2:3); 

% u_cg = (u(cell_v{i}(1)) + u(cell_v{i}(2)) + u(cell_v{i}(3)))./3;
% B = sqrt(epsilon^2 + u_cg^2);

%B = sqrt(epsilon^2 + norm(grad_u_cg,2)^2);

B = sqrt(epsilon^2 + grad_u(1)^2 + grad_u(2)^2);

A = 1/B;

end
