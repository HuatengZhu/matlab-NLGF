%% Local stiffness matrix of P1 FE

function M = stima(mpedges,cell_e,area,r)

%% Area of cell r 
% mk = 1/2 * det([ones(1,3);vertices]);
mk = area(r);

%% Barycentric coordinates of cell r
midpoints = [mpedges(cell_e{r}(1),:)' mpedges(cell_e{r}(2),:)' mpedges(cell_e{r}(3),:)'];
L1 = [ones(1,3); midpoints]'\ [1;0;0];
	% L1 are the coordinates of the 1st barycentric coordinate lambda_1 (associated to
	%		the 1st vertex). lambda_1 is a function of x,y, given by
	%		lambda_1(x,y)= alpha0 + alpha1 x + alpha2 y
	%		if L1=[alpha0 alpha1 alpha2]
L2 = [ones(1,3); midpoints]'\[0;1;0];
L3 = [ones(1,3); midpoints]'\[0;0;1];

%% Element stiffness matrix
M = mk * [ [L1(2);L1(3)]'*[L1(2);L1(3)] [L1(2);L1(3)]'*[L2(2);L2(3)] [L1(2);L1(3)]'*[L3(2);L3(3)]
         [L2(2);L2(3)]'*[L1(2);L1(3)] [L2(2);L2(3)]'*[L2(2);L2(3)] [L2(2);L2(3)]'*[L3(2);L3(3)]
         [L3(2);L3(3)]'*[L1(2);L1(3)] [L3(2);L3(3)]'*[L2(2);L2(3)] [L3(2);L3(3)]'*[L3(2);L3(3)]];
end

