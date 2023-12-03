%Nose cone function and tangent approximation
clc
close all
clear all

L = 0.52; %nose length
R = 0.05 ; %nose base radius
C= 0; %haack series parameter, 0 for Von Karman
n = 1/2; %power series parameter
x = linspace(0,L,1000);

%nosecone shape function
%von karman
%y = @(x) (R/sqrt(pi))*sqrt(acos(1-(2.*x)/L)-(sin(2.*acos(1-(2.*x)/L)))/2+C*(sin(acos(1-(2.*x)/L))).^3);
%power series
y = @(x) R*(x./L).^n;
ys = y(x);

%slope of the tangent lines
%von karman
%dy_by_dx = @(x) -(R.*((2.*cos(2.*acos(1 - (2.*x)./L)))./(L.*(1 - ((2.*x)./L - 1).^2.).^(1/2)) - 2./(L.*(1 - ((2.*x)./L - 1).^2).^(1/2)) + (6.*C.*(1 - ((2.*x)./L - 1).^2).^(1/2).*((2.*x)./L - 1))./L))./(2.*sqrt(pi).*(acos(1 - (2.*x)./L) - sin(2.*acos(1 - (2.*x)./L))/2 + C.*(1 - ((2.*x)./L - 1).^2).^(3/2)).^(1/2));
%power series
dy_by_dx = @(x) n*(R/L^n).*x.^(n-1);
y_f = dy_by_dx(x);

plot(x,ys,"k","LineWidth",1)
hold on
plot(x,-ys,"k","LineWidth",1)
axis([0 L -0.1 0.1])

S = 15; % Desidered Segments
N = S + 1; % Number of Nodes
delta = (max(y_f)-min(y_f))/(N); %gradient delta
x_t = zeros(N,1); %x position of the tangent nodes
y_t = zeros(N,1); %y position of the tangent nodes
x_t(1) = 0; %the tip is set to be a tangent node
y_t(1) = 0;
k = 2;%number of tangent nodes alreay found (tip+base)

for i = 1:length(y_f)
    
    if  y_f(i) <= max(y_f) - delta && k == 2
        x_t(k) = x(i);
        k = k + 1;
       
    elseif (y_f(i) <= (dy_by_dx(x_t(k - 1)) - delta)) && (k > 2)
        x_t(k) = x(i);
        k = k + 1;
    end
end

x_t(N) = L;%the base is set to be a tangent node
y_t(N) = R;
y_t = y(x_t);
%to solve numerical cancellation problem that might leave some zeros in tne
%x_t array:
k = 2;
while k<N
    if x_t(k) == 0 %if a 0 is found in the array (in a position different for 1, because k starts at 2)
        ns = N-k;%number of nodes to correct
        d = (x_t(N)-x_t(k-1))/(ns+1);%the 0 nodes are substituted by equispaced nodes of length d
        for j = 1:ns
            x_t(k+j-1)=x_t(k-1)+d*j;
            y_t(k+j-1) = y(x_t(k+j-1));
        end
        k = N+1;%to break the loop
    end
    k = k+1;
end
%Method 1: join the tangent nodes. Tangent nodes coincident to the
%final nodes. Red segments
plot(x_t,y_t,"bo", "LineWidth",2)
plot(x_t,-y_t, "bo", "LineWidth",2)
for i = 1:S
    x1 = x_t(i);
    x2 = x_t(i + 1);
    plot([x1;x2], [y(x1);y(x2)],'r','linewidth',2)
    plot([x1;x2], [-y(x1);-y(x2)],'r','linewidth',2)
end

%Method 2: find the final nodes by intersecting the tangents to the
%already found tangent nodes. Green segments.
x_n = zeros(N+1,1); %x position of the final nodes
y_n = zeros(N+1,1); %y position of the final nodes
for j = 1:(length(x_t)-1)
    x1 = x_t(j)
    x2 = x_t(j+1)
    y1 = y_t(j)
    y2 = y_t(j+1)
    m1 = n*(R/L^n)*x1^(n-1)%slopes of the tangents in the tangent nodes
    m2 = dy_by_dx(x2)
    x = (y2-y1+m1*x1-m2*x2)/(m1-m2);
    y = m1*(x-x1)+y1;
    x_n(j+1)=x;
    y_n(j+1)=y;
    
end
x_n(length(x_t)+1) = L;
y_n(length(y_t)+1) = R;

for k=1:length(x_n)-1
    plot([x_n(k);x_n(k+1)], [y_n(k);y_n(k+1)],'g','linewidth',2)
    plot([x_n(k);x_n(k+1)], [-y_n(k);-y_n(k+1)],'g','linewidth',2)
end
