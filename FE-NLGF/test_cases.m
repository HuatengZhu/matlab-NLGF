function u =test_cases(t,z,ucase)

x=z(:,1);
y=z(:,2);
u=zeros(size(z,1),1); 

if (ucase==1)
    u=(cos(t)*(cos(pi*x).*cos(pi*y)))'; 

elseif (ucase==2)
        u=(cos(t)* abs(x).^3 .* (1-x).^2 .* abs(y).^3 .* (1-y).^2)'; 

end
end
