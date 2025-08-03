function f = fs(z,t,epsilon,ucase)

x = z(:,1);
y = z(:,2);

if (ucase==1)
    f = (2*pi^2*cos(pi*x)*cos(pi*y)*cos(t))/(pi^2*cos(x*pi)^2*sin(y*pi)^2*cos(t)^2 + pi^2*cos(y*pi)^2*sin(x*pi)^2*cos(t)^2 + epsilon^2)^(1/2) - cos(pi*x)*cos(pi*y)*sin(t) - (pi*cos(pi*y)*sin(pi*x)*cos(t)*(2*pi^3*cos(pi*x)*cos(pi*y)^2*sin(pi*x)*cos(t)^2 - 2*pi^3*cos(pi*x)*sin(pi*x)*sin(pi*y)^2*cos(t)^2))/(2*(pi^2*cos(x*pi)^2*sin(y*pi)^2*cos(t)^2 + pi^2*cos(y*pi)^2*sin(x*pi)^2*cos(t)^2 + epsilon^2)^(3/2)) - (pi*cos(pi*x)*sin(pi*y)*cos(t)*(2*pi^3*cos(pi*x)^2*cos(pi*y)*sin(pi*y)*cos(t)^2 - 2*pi^3*cos(pi*y)*sin(pi*x)^2*sin(pi*y)*cos(t)^2))/(2*(pi^2*cos(x*pi)^2*sin(y*pi)^2*cos(t)^2 + pi^2*cos(y*pi)^2*sin(x*pi)^2*cos(t)^2 + epsilon^2)^(3/2));
elseif (ucase==2)
    f = (((2*x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 4*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)*(x^2*y^2*cos(t)*(2*x - 2)*(2*y - 2) + 4*x*y*cos(t)*(x - 1)^2*(y - 1)^2 + 2*x*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(2*x - 2)*(y - 1)^2) + (2*x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 4*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)*(2*x^2*y^2*cos(t)*(x - 1)^2 + 2*x^2*cos(t)*(x - 1)^2*(y - 1)^2 + 4*x^2*y*cos(t)*(2*y - 2)*(x - 1)^2))*(x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2))/(2*((x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 2*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)^2 + (x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)^2 + epsilon^2)^(3/2)) - (2*x^2*y^2*cos(t)*(y - 1)^2 + 2*y^2*cos(t)*(x - 1)^2*(y - 1)^2 + 4*x*y^2*cos(t)*(2*x - 2)*(y - 1)^2)/((x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 2*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)^2 + (x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)^2 + epsilon^2)^(1/2) - (2*x^2*y^2*cos(t)*(x - 1)^2 + 2*x^2*cos(t)*(x - 1)^2*(y - 1)^2 + 4*x^2*y*cos(t)*(2*y - 2)*(x - 1)^2)/((x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 2*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)^2 + (x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)^2 + epsilon^2)^(1/2) + (((2*x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 4*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)*(x^2*y^2*cos(t)*(2*x - 2)*(2*y - 2) + 4*x*y*cos(t)*(x - 1)^2*(y - 1)^2 + 2*x*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(2*x - 2)*(y - 1)^2) + (2*x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 4*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)*(2*x^2*y^2*cos(t)*(y - 1)^2 + 2*y^2*cos(t)*(x - 1)^2*(y - 1)^2 + 4*x*y^2*cos(t)*(2*x - 2)*(y - 1)^2))*(x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 2*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2))/(2*((x^2*y^2*cos(t)*(2*x - 2)*(y - 1)^2 + 2*x*y^2*cos(t)*(x - 1)^2*(y - 1)^2)^2 + (x^2*y^2*cos(t)*(2*y - 2)*(x - 1)^2 + 2*x^2*y*cos(t)*(x - 1)^2*(y - 1)^2)^2 + epsilon^2)^(3/2)) - x^2*y^2*sin(t)*(x - 1)^2*(y - 1)^2;
end

end
 

% u_t = div(\nabla u/ \sqrt{\epsilon^2 + |\nabla u|^2} + f
% syms t x y
% 
% u = cos(t)*(cos(pi*x).*cos(pi*y)); 
% 
% u_t = diff(u,t);
% u_x = diff(u,x);
% u_y = diff(u,y);
% frac_u = 1/sqrt(epsilon^2 + u_x^2 + u_y^2);
% 
% div = diff(u_x * frac_u, x) + diff(u_y * frac_u, y); 
% 
% f = u_t - div
% 
% syms t x y
% 
% u = cos(t) * x.^2 .* (1-x).^2 .* y.^2 .* (1-y).^2; 
% 
% u_t = diff(u,t);
% u_x = diff(u,x);
% u_y = diff(u,y);
% frac_u = 1/sqrt(epsilon^2 + u_x^2 + u_y^2);
% 
% div = diff(u_x * frac_u, x) + diff(u_y * frac_u, y); 
% 
% f = u_t - div

% % u_t = div(\nabla u)
% syms t x y
% 
% u = cos(t)*(cos(pi*x).*cos(pi*y)); 
% 
% u_t = diff(u,t);
% u_x = diff(u,x);
% u_y = diff(u,y);
% div = diff(u_x,x) + diff(u_y,y);
% f = u_t - div 
