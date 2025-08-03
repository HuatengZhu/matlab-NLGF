% Assemble the global diffusion matrix
function [J, F] = assemble_jacobian_system(cell_v, cell_e, ncell, vertex, nedge, area, mpedges, u,epsilon)


%% Initialise vectors for sparse matrix
% Evaluate number of non-zeros entries (cf below how many times we do "pos=pos+1")
nz=9*ncell;
IS=zeros(nz,1);
JS=zeros(nz,1);
VS1=zeros(nz,1);
VS2=zeros(nz,1);


% "pos"=position inside the vectors IA, JA, VA that store the entries of A
pos=0;

for i=1:ncell
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

    B = sqrt(epsilon^2 + grad_u(1)^2 + grad_u(2)^2);

    A = 1/B;

    D = eye(2)/B - (grad_u * grad_u')/B^3;

    Xloc = A * stima(mpedges,cell_e,area,i);

    % element stiffness matrix
    Sloc=area(i)*[(D*[L1(2);L1(3)])'*[L1(2);L1(3)] (D*[L1(2);L1(3)])'*[L2(2);L2(3)] (D*[L1(2);L1(3)])'*[L3(2);L3(3)]
        (D*[L2(2);L2(3)])'*[L1(2);L1(3)] (D*[L2(2);L2(3)])'*[L2(2);L2(3)] (D*[L2(2);L2(3)])'*[L3(2);L3(3)]
        (D*[L3(2);L3(3)])'*[L1(2);L1(3)] (D*[L3(2);L3(3)])'*[L2(2);L2(3)] (D*[L3(2);L3(3)])'*[L3(2);L3(3)]];

    % Loop over edges
    for jj=1:3
        jedge = cell_e{i}(jj);
        for kk=1:3
            kedge = cell_e{i}(kk);
            pos=pos+1;
            IS(pos) = jedge;
            JS(pos) = kedge;
            VS1(pos) = Sloc(jj,kk);
            VS2(pos) = Xloc(jj,kk);
        end
    end
end
J = sparse(IS(1:pos),JS(1:pos),VS1(1:pos),nedge,nedge);
F = sparse(IS(1:pos),JS(1:pos),VS2(1:pos),nedge,nedge);

