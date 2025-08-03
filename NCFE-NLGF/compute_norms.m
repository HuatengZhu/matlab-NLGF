%% Compute L2 norm in space of vector v in the non-conforming P1 space

function [L2norm, W11norm, H1norm] = compute_norms(v, area, mpedges, cell_e)

L2norm = 0;
W11norm = 0;
H1norm = 0;
ncell = size(area, 1);

for r = 1:ncell
    % The L2 norm is approximated by taking the value at the centers of
    % mass, which is the average values at the vertices
    v_center = (1/3)*(v(cell_e{r}(1)) + v(cell_e{r}(2)) + v(cell_e{r}(3)));
    L2norm = L2norm + area(r) * v_center^2;

    % Compute coefficients of barycentric coordinate functions (see stima.m)
    midpoints = [mpedges(cell_e{r}(1),:)' mpedges(cell_e{r}(2),:)' mpedges(cell_e{r}(3),:)'];
    L1=[ones(1,3);midpoints]'\[1;0;0];
    L2=[ones(1,3);midpoints]'\[0;1;0];
    L3=[ones(1,3);midpoints]'\[0;0;1];

    grad_v=v(cell_e{r}(1))*L1(2:3) + v(cell_e{r}(2))*L2(2:3) + v(cell_e{r}(3))*L3(2:3);

    W11norm = W11norm + area(r) * norm(grad_v); 
    H1norm = H1norm + area(r) * norm(grad_v)^2;
end

%% Square root to get the actual L2 norm
L2norm = sqrt(L2norm);
H1norm = sqrt(H1norm);