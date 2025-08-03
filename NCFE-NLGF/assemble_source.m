% Assemble the global diffusion matrix
function F = assemble_source(cell_e, ncell, nedge, area, cg, t, epsilon, ucase)

F = zeros(nedge,1);  

    for i=1:ncell 
        F(cell_e{i}(1:3))=F(cell_e{i}(1:3))+(area(i)/3)* fs(cg(i,:),t,epsilon,ucase);
    end
end