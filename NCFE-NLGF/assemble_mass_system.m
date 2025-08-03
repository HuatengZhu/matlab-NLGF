% Assemble the global mass matrix
function MassMat = assemble_mass_system(area, ncell, cell_e, nedge);


%% Initialise vectors for sparse matrix
% Evaluate number of non-zeros entries (cf below how many times we do "pos=pos+1")
nz=9*ncell;
IM=zeros(nz,1);
JM=zeros(nz,1);
VM=zeros(nz,1);

% "pos"=position inside the vectors IA, JA, VA that store the entries of A
pos=0;

%% Loop over cells on index i
for i=1:ncell
  %% Mass matrix on each cell
  Mloc = area(i)/12 .* [2, 1, 1; 
                        1, 2, 1; 
                        1, 1, 2];
  
  %% Loop over vertices on each cell
  for jj=1:3
    jedge = cell_e{i}(jj); 
      for kk=1:3
        kedge = cell_e{i}(kk);
        pos=pos+1;
        IM(pos) = jedge;
        JM(pos) = kedge;
        VM(pos) = Mloc(jj,kk);%It does not overwrite by definition
      end
  end
end

%% Creation of the sparse matrix
MassMat=sparse(IM(1:pos),JM(1:pos),VM(1:pos),nedge,nedge); % VM is a vector that assigns values to  IM \times JM

