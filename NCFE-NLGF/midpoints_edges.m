% Computes the midpoint of all edges
function mpedges = midpoints_edges(ncell,nedge,cell_e,cell_v,vertex);

mpedges=zeros(nedge,2);

for i=1:ncell
  for j=1:size(cell_e{i},2)
    mpedges(cell_e{i}(j),:) = (1/2) * (vertex(cell_v{i}(j),:) + vertex(cell_v{i}(j+1),:));
  end;
end;


