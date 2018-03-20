function [ result ] = zerosofNo2( V,a )

%FIND_ZEROS_2 此处显示有关此函数的摘要
%   此处显示详细说明
%   V represents the depth of the well
%   a represents the width of the well

E=linspace(0,3*V,V*50);

evenfun=@(E)sqrt(V-E).*tan(a*sqrt(V-E))-sqrt(E);
oddfun=@(E)sqrt(V-E).*cot(a*sqrt(V-E))+sqrt(E);
evenfun1=feval(evenfun,E);
oddfun1=feval(oddfun,E);
EE_even=[0;0];
EE_odd=[0;0];
result=[0;0];


% find the approximate intervals of solution
j=1;k=1;
for i=1:length(E)-1
    
    if evenfun1(i)*evenfun1(i+1)<0
        EE_even(1,j)=E(i);
        EE_even(2,j)=E(i+1);
        j=j+1;
    end
    if oddfun1(i)*oddfun1(i+1)<0
        EE_odd(1,k)=E(i);
        EE_odd(2,k)=E(i+1);
        k=k+1;
    end
end

% use bisection method to solve the function in each interval
% even part
j=1;
for n=1:length(EE_even(1,:))
    for i=1:35
        mid=mean([EE_even(1,n),EE_even(2,n)]);
        if evenfun(mid)*evenfun(EE_even(2,n))<0
            EE_even(1,n)=mid;
        else
            EE_even(2,n)=mid;
        end
    end
    
    if abs(evenfun(mid))<10^-9
        result(1,j)=mid;
        j=j+1;
    end
end

% odd part
j=1;
for n=1:length(EE_odd(1,:))

    for i=1:35
        mid=mean([EE_odd(1,n),EE_odd(2,n)]);
        if oddfun(mid)*oddfun(EE_odd(2,n))<0
            EE_odd(1,n)=mid;
        else
            EE_odd(2,n)=mid;
        end
    end
   
    if abs(oddfun(mid))<10^-9
        result(2,j)=mid;
        j=j+1;
    end
end
    result(result==0)=NaN;

end

