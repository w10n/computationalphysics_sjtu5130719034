function [I] = Gauss_quad_4(f,a,b)
%this function use Gauss method with 4 points to calculate the Integration 
%from a to b
%f is the function you would like to do the interval
%a and b is the bound of your integration interval
%I returns the value of the intergration

%find the abscissas and its weight
A=roots([35,0,-30,0,3]);
A=sort(A);
A=A';
B=[];
for i =1:4
    B=[B;A(1)^(i-1),A(2)^(i-1),A(3)^(i-1),A(4)^(i-1)];
end
bb=[2,0,2/3,0]';
w=B\bb;
w=w';

AA=(b+a)/2+(b-a)/2.*A;
ww=(b-a)/2.*w;
I=sum(f(AA).*ww);