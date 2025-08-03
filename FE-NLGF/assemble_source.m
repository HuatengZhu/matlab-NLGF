% Assemble the global diffusion matrix
function F = assemble_source(cell_v, ncell, nvert, area, cg, t, epsilon, ucase)

F = zeros(nvert,1);  

    for i=1:ncell 
        F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+(area(i)/3)* fs(cg(i,:),t,epsilon,ucase);
    end
end