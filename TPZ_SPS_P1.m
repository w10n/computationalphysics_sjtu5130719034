function [ lt,ls ] = TPZ_SPS_P1( N )

% initialize
x=linspace(0,pi/4,N+1);
y=sqrt(1+cos(x).^2);
dx=pi/4/N;

% trapezoid
lt=0;
for i=1:N
    % calculate trapezoid area
    lt=lt+(y(i)+y(i+1))*dx/2;
end

% Simpson
% construct vector a
a(1)=dx/3;
for i=1:N/2
    a(2*i)=4*dx/3;
    a(2*i+1)=2*dx/3;
end
a(N+1)=dx/3;
ls=sum(a.*y);


end

