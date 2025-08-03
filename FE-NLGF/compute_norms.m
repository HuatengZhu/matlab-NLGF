%% Compute L2 norm in space of vector v in the P1 space and L1 norm in space of \nabla v

function [L2norm, W11norm, H1norm] = compute_norms(v, area, vertex, cell_v)

L2norm = 0;
W11norm = 0;
H1norm = 0;
ncell = size(area, 1);

for r = 1:ncell
    % The L2 norm is approximated by taking the value at the centers of
    % mass, which is the average values at the vertices
    v_center = (1/3)*(v(cell_v{r}(1)) + v(cell_v{r}(2)) + v(cell_v{r}(3)));
    L2norm = L2norm + area(r) * v_center^2;

    % Compute coefficients of barycentric coordinate functions (see stima.m)
    vertices=[vertex(cell_v{r}(1),1) vertex(cell_v{r}(2),1) vertex(cell_v{r}(3),1);
              vertex(cell_v{r}(1),2) vertex(cell_v{r}(2),2) vertex(cell_v{r}(3),2)];
    L1=[ones(1,3);vertices]'\[1;0;0];
    L2=[ones(1,3);vertices]'\[0;1;0];
    L3=[ones(1,3);vertices]'\[0;0;1];

    grad_v=v(cell_v{r}(1))*L1(2:3) + v(cell_v{r}(2))*L2(2:3) + v(cell_v{r}(3))*L3(2:3);

    W11norm = W11norm + area(r) * norm(grad_v); 
    H1norm = H1norm + area(r) * norm(grad_v)^2;
end

%% Square root to get the actual L2 norm
L2norm = sqrt(L2norm);
H1norm = sqrt(H1norm);