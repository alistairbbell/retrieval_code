% x_hat=convolute_profile_with_avk(x,x_a,A,grid);
%
% x_hat=x_a+A*(x-x_a);
%
% x     = [p [Pa] / z [m], vmr [-]];
% x_a   = [p [Pa] / z [m], vmr [-]];
% x_hat = [p [Pa] / z [m], vmr [-]];
% grid  = 'p'-> pressure, 'z'-> altitude grid
%
% A has to be given on the grid of x_a and is
% interpolated to the grid of x, conserving
% the measurement response:  
%
% sum(A(i,:))=sum(A_i(i,:))

function x_hat=convolute_profile_with_avk(x,x_a,A,grid)

% interpolate x_a and A to the grid of x
if grid=='p'
    x_a_i=interp1(log(x_a(:,1)),x_a(:,2),log(x(:,1)));
    A_i=interp1(log(x_a(:,1)),A',log(x(:,1)))';
elseif grid=='z'
    x_a_i=interp1(x_a(:,1),x_a(:,2),x(:,1));
    A_i=interp1(x_a(:,1),A',x(:,1))';
else
    disp('you have to define the grid: grid=''z''/''p'' ');
    return
end

% remove NaN's
ind=find(isnan(x_a_i)==0 & isnan(x(:,2))==0);
x_a_i=x_a_i(ind);
A_i=A_i(:,ind);

% rescale each line of A_i in order to conserve the msm response
normfac=sum(A,2)./sum(A_i,2);

for i=1:length(x_a)
	A_i(i,:)=A_i(i,:)*normfac(i);
end

% convolute x with A
x_hat(:,1)=x_a(:,1);
x_hat(:,2)=x_a(:,2)+A_i*(x(ind,2)-x_a_i);
