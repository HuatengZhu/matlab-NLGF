
% Assemble the global diffusion matrix
function [J, F] = assemble_jacobian_system(cell_v,ncell,nvert,vertex, u,area, epsilon)


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
    % We are collecting vertices for each cell in the following to compute local stiffness matrix
    vertices=[vertex(cell_v{i}(1),1) vertex(cell_v{i}(2),1) vertex(cell_v{i}(3),1);
        vertex(cell_v{i}(1),2) vertex(cell_v{i}(2),2) vertex(cell_v{i}(3),2)]; % this is a 2 X 3 matrix whose j-th column is the coordinate of vertex cell_v{i}(j).

    %% Barycentric coordinates of cell i
    L1 = [ones(1,3); vertices]'\ [1;0;0];
    % L1 are the coordinates of the 1st barycentric coordinate lambda_1 (associated to
    %		the 1st vertex). lambda_1 is a function of x,y, given by
    %		lambda_1(x,y)= alpha0 + alpha1 x + alpha2 y
    %		if L1=[alpha0 alpha1 alpha2]
    L2 = [ones(1,3); vertices]'\[0;1;0];
    L3 = [ones(1,3); vertices]'\[0;0;1];

    % gradient at the gravity centre
    grad_u = u(cell_v{i}(1))*L1(2:3) + u(cell_v{i}(2))*L2(2:3) + u(cell_v{i}(3))*L3(2:3); 

    % A = diffusion_coefficient(u, i, epsilon, vertices, cell_v);

    % Weighted local stiffness matrix
    
    B = sqrt(epsilon^2 + grad_u(1)^2 + grad_u(2)^2);

    A = 1/B;

    D1 = eye(2)/B;
    D2 = - (grad_u * grad_u')/B^3;
    D = D1 + D2;

    % element stiffness matrix
    % Sloc=area(i)*[(D*[L1(2);L1(3)])'*[L1(2);L1(3)] (D*[L1(2);L1(3)])'*[L2(2);L2(3)] (D*[L1(2);L1(3)])'*[L3(2);L3(3)]
    %     (D*[L2(2);L2(3)])'*[L1(2);L1(3)] (D*[L2(2);L2(3)])'*[L2(2);L2(3)] (D*[L2(2);L2(3)])'*[L3(2);L3(3)]
    %     (D*[L3(2);L3(3)])'*[L1(2);L1(3)] (D*[L3(2);L3(3)])'*[L2(2);L2(3)] (D*[L3(2);L3(3)])'*[L3(2);L3(3)]];
    Sloc1=area(i) * A .*[
        [L1(2);L1(3)]'*[L1(2);L1(3)] [L1(2);L1(3)]'*[L2(2);L2(3)] [L1(2);L1(3)]'*[L3(2);L3(3)]
        [L2(2);L2(3)]'*[L1(2);L1(3)] [L2(2);L2(3)]'*[L2(2);L2(3)] [L2(2);L2(3)]'*[L3(2);L3(3)]
        [L3(2);L3(3)]'*[L1(2);L1(3)] [L3(2);L3(3)]'*[L2(2);L2(3)] [L3(2);L3(3)]'*[L3(2);L3(3)]
        ];
    Sloc2=area(i) * A^3 .*[
        (grad_u'*[L1(2);L1(3)]).*(grad_u'*[L1(2);L1(3)]) (grad_u'*[L1(2);L1(3)]).*(grad_u'*[L2(2);L2(3)]) (grad_u'*[L1(2);L1(3)]).*(grad_u'*[L3(2);L3(3)])
        (grad_u'*[L2(2);L2(3)]).*(grad_u'*[L1(2);L1(3)]) (grad_u'*[L2(2);L2(3)]).*(grad_u'*[L2(2);L2(3)]) (grad_u'*[L2(2);L2(3)]).*(grad_u'*[L3(2);L3(3)])
        (grad_u'*[L3(2);L3(3)]).*(grad_u'*[L1(2);L1(3)]) (grad_u'*[L3(2);L3(3)]).*(grad_u'*[L2(2);L2(3)]) (grad_u'*[L3(2);L3(3)]).*(grad_u'*[L3(2);L3(3)])
        ];
    Sloc = Sloc1 - Sloc2;
    % Sloc = epsilon^2 * A^3 * stima(vertices,i);
    
     Xloc = A * stima(vertices,i);

    % Loop over vertices
    for jj=1:3
        jvert = cell_v{i}(jj);
        for kk=1:3
            kvert = cell_v{i}(kk);
            pos=pos+1;
            IS(pos) = jvert;
            JS(pos) = kvert;
            VS1(pos) = Sloc(jj,kk);%It does not overwrite by definition
            VS2(pos) = Xloc(jj,kk);%It does not overwrite by definition
        end
    end
end

%% Creation of the sparse matrix
J=sparse(IS(1:pos),JS(1:pos),VS1(1:pos),nvert,nvert);
F=sparse(IS(1:pos),JS(1:pos),VS2(1:pos),nvert,nvert);
end

