%% Compute the interpolate in the P1 space of the exact solution
function interp=interpolate_ex_sol(t,vertex, Ndt(imesh))

nvert = size(vertex,1);
interp_entry(t) = zeros(nvert,1);
interp = zeros(Ndt(imesh), nvert);
for j = 1 : Ndt(imesh) + 1
    for i = 1 : nvert
        interp_entry=test_cases((j-1) * dt,vertex(i,:));
    end
    interp(j,:) = interp_entry;
end
