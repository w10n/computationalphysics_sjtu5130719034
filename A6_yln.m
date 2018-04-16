%A8
%% No.1
clc;clear;
% q.1
% 3 points, 2 intervals
[lt,ls]=TPZ_SPS_P1(2);
fprintf('the length of the curve is:\n %5.5f (trapezoid method)\nand %5.5f (simpson method)\n',lt,ls)

% q.2
% theoretical value
syms x y
y=sqrt(1+cos(x)^2);
l=eval(int(y,x,0,pi/4));
fprintf('\nthe theoretical value is %5.5f\n',l);
fprintf('Simpson is larger and Trepezoid is smaller\n\n');

% q.3
% construct two matrixes to present the result
outputhead=['     N','           ¦ÅT','           ¦ÅS'];
[lt,ls]=TPZ_SPS_P1(2);
output=[2,abs(lt-l)/l,abs(ls-l)/l];
for i=1:15
N=5*2^i;
[lt,ls]=TPZ_SPS_P1(N);
output(i+1,:)=[N,abs(lt-l)/l,abs(ls-l)/l];
end
format shortE
disp(outputhead)
disp(output)

% q.4
plot(log(output(:,1)),log(output(:,2)),log(output(:,1)),log(output(:,3)));
xlabel('log(N)')
ylabel('log(error)')
legend('Trapezoid method','Simpsom method')

% q.5
% from the figure we plotted in q.4,
% we can observe that it appears as a straight line on a log-log plot.
% with N increases, the logarithm of relative error decreases,
% and the error of Simpson method decreases faster.
% we can also see that,
% the error of simpson method stops to decrease because of the round off error.
% so, with a relative error at about 10^-15, or after 15 decimal places, 
% the machine precision has been met.

%% No.2
%(a)
f=@(x) sqrt(1+(cos(x)).^2);
format long
figure(1)
for n=1:1:5
    N(n)=2^(n-1);
    I(n)=Simpson2(f,0,pi/4,N(n));
    error(n)=I(n)-integral(f,0,pi/4);
    h4(n)=((pi/4)/(3*N(n)))^4;
end
Table=[N;I;error];
fprintf('Corresponding value matrix is :\n')
disp(Table)
p=polyfit(h4,error,1);x=linspace(0,h4(1),100);y=p(1)*x+p(2);
plot(x,y,'b',h4,error,'*r'),title('fitting curve between error and h^4')
%(b)
k=1;N=0;
while k>5*10^(-9)
    N=N+1;
    I0 = Simpson2(f,0,pi/4,N);
    I2 = integral(f,0,pi/4);
    k=abs((I0-I2)/I2);
end
h=(pi/4)/(3*N);
fprintf('N=%d and h=%d so that the composite Simpson 3/8 rule can be used to compute the intergral with an acccuracy of 5*10^(-9)\n\n',N,h)
%QUESTION?
I0 = Simpson2(f,0,pi/4,12);%3/8 simpson method 
I1 = Simpson(f,0,pi/4,36);
I2 = integral(f,0,pi/4); % matlab built-in command
n1=(I0-I2)/(I1-I2);

%% No.3
clc
clear
close
format long
%a
%Find the four abscissas
A=roots([35,0,-30,0,3]);%use the legendre polynomial of degree 4 to find the four abscissas
A=sort(A);%rearrange the sequence ofthe roots A
A=A';

%Find the weights of each abscissas
B=[];
for i =1:4
    B=[B;A(1)^(i-1),A(2)^(i-1),A(3)^(i-1),A(4)^(i-1)];
end
b=[2,0,2/3,0]';
w=B\b;
w=w';
fprintf('the four abscissas are%0.4f\n',A)
fprintf('the four weight are%0.4f\n',w)
%c
f=@(x) (1+cos(x).^2).^0.5;
I=Gauss_quad_4(f,0,pi/4);
fprintf('the intergation of problem1 is %0.10f\n',I)

%% No.4
close
clc
clear
format long
rv=0.062204682443299;
f=@(x) exp(-x).*sin(x)./(x.^3+1);
zz=1:10^-3:2;
zzy=f(zz);
plot(zz,zzy);
xlabel('x');ylabel('y');title('exp(-x)*sin(x)/(x^3+1)')
s=Simpson(f,1,2,4);
se=abs((s-rv)/rv);
g=Gauss_quad_4(f,1,2);
ge=abs((g-rv)/rv);
r=Romberg(f,1,2,2,3);
re=abs((r(1,3)-rv)/rv);
syms t;
i=double(int(exp(-t).*sin(t)./(t.^3+1),1,2));
fprintf('The intergation by Simpson method with n=4 is %0.15f\n',s)
fprintf('The intergation by Gauss_quad method with n=4 is %0.15f\n',g)
fprintf('The intergation by Romberg method with n=2 and nlevel=3 is %0.15f\n',r(1,3))
fprintf('The intergation by built-in function int is %0.15f\n',i)
fprintf('the relative error by Simpson method is %0.15f\n',se)
fprintf('the relative error by Gauss_quad method is %0.15f\n',ge)
fprintf('the relative error by Romberg method is %0.15f\n',re)


