 %Nose cone function and tangent approximation
clc
close all
clear all

L = 0.52; %nose length
R = 0.05 ; %nose base radius
C= 0; %haack series parameter, 0 for Von Karman
x = linspace(0,L,100);

%nosecone shape function
y = @(x) (R/sqrt(pi))*sqrt(acos(1-(2.*x)/L)-(sin(2.*acos(1-(2.*x)/L)))/2+C*(sin(acos(1-(2.*x)/L))).^3);
ys = y(x);

%slope of the tangent lines
dy_by_dx = @(x) -(R.*((2.*cos(2.*acos(1 - (2.*x)./L)))./(L.*(1 - ((2.*x)./L - 1).^2.).^(1/2)) - 2./(L.*(1 - ((2.*x)./L - 1).^2).^(1/2)) + (6.*C.*(1 - ((2.*x)./L - 1).^2).^(1/2).*((2.*x)./L - 1))./L))./(2.*sqrt(pi).*(acos(1 - (2.*x)./L) - sin(2.*acos(1 - (2.*x)./L))/2 + C.*(1 - ((2.*x)./L - 1).^2).^(3/2)).^(1/2));
y_f = dy_by_dx(x);

plot(x,ys,"k","LineWidth",1)
hold on
plot(x,-ys,"k","LineWidth",1)

S = 13; % Desidered Segments
N = S + 1; % Number of Nodes
delta = (max(y_f)-min(y_f))/(N);
x_t = zeros(N,1);
y_t = zeros(N,1);
x_t(1) = 0;
y_t(1) = 0;
k = 2;

for i = 1:length(y_f)
    
    if  y_f(i) <= max(y_f) - delta && k == 2
        x_t(k) = x(i);
        k = k + 1;
       
    elseif (y_f(i) <= (dy_by_dx(x_t(k - 1)) - delta)) && (k > 2)
        x_t(k) = x(i);
        k = k + 1;
    end
end

x_t(N) = L;
y_t(N) = R;
y_t = y(x_t);
y_t_negativo = -1.*y_t;

plot(x_t,y_t,"bo", "LineWidth",2)
plot(x_t,y_t_negativo, "bo", "LineWidth",2)
k = 2;

for i = 1:S
    x1 = x_t(i);
    x2 = x_t(i + 1);
    y_tangente = @(x) (x-x1)/(x2-x1)*(y(x2)-y(x1)) + y(x1);
    z_i = linspace(x1,x2);
    plot(z_i,y_tangente(z_i),"b","LineWidth",2)
    k = k + 1;
end







