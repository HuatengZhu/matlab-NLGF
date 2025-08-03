function u =test_cases(t,z,ucase)

% global ucase;

x=z(:,1);
y=z(:,2);
u=zeros(size(z,1),1);
%nargout is for the number of return values from [u,ut,gradu,D2u]
% if (nargout>1)
%     ut=zeros(size(z,1),1);
%     if (nargout>2)
%         gradu=zeros(size(z,1),2);
%         % if (nargout>3)
%         %     D2u=zeros(2,2,size(z,1));
%         % end
%     end
% end

if (ucase==1)
    u=(cos(t)*(cos(pi*x).*cos(pi*y)))';
    % if (nargout>1)
    %     ut(:,1)=-sin(t)*(cos(pi*x).*cos(pi*y));
    %     if (nargout>2)
    %         gradu(:,1)=cos(t)*(-pi*sin(pi*x).*cos(pi*y));
    %         gradu(:,2)=cos(t)*(-pi*cos(pi*x).*sin(pi*y));
    %         if (nargout>3)
    %             D2u(1,1,:)=cos(t)*(-pi^2*cos(pi*x).*cos(pi*y));
    %             D2u(1,2,:)=cos(t)*(pi^2*sin(pi*x).*sin(pi*y));
    %             D2u(2,1,:)=D2u(1,2,:);
    %             D2u(2,2,:)=D2u(1,1,:);
    %         end
    %     end
    % end

elseif (ucase==2)
     u=(exp(t)*(cos(pi*x).*cos(pi*y)))';
        % u = (cos(t) * x.^2 .* (1-x).^2 .* y.^2 .* (1-y).^2)'; 
        % if (nargout>1)
        %     ut(:,1)=-sin(t) * x.^2 .* (1-x).^2 .* y.^2 .* (1-y).^2;
        %     if (nargout>2)
        %         gradu(:,1)=cos(t) * (2*x.*(1-x).^2-2*x.^2.*(1-x)).*y.^2.*(1-y.^2);
        %         gradu(:,2)=cos(t) * (2*y.*(1-y).^2-2*y.^2.*(1-y)).*x.^2.*(1-x.^2);
        %     end
        % end

end
end
