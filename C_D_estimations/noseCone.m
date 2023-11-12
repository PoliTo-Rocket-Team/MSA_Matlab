%Nose cone function and tangent approximation
close all
clear all
L = 0.52; %nose length
R = 0.05 ; %nose base radius
C= 0; %haack series parameter, 0 for Von Karman
x = linspace(0,L,100);

%nosecon shape function
y = @(x) (R/sqrt(pi))*sqrt(acos(1-(2.*x)/L)-(sin(2.*acos(1-(2.*x)/L)))/2+C*(sin(acos(1-(2.*x)/L))).^3);
ys = y(x);

%slope of the tangent lines
dy_by_dx = @(x) -(R.*((2.*cos(2.*acos(1 - (2.*x)./L)))./(L.*(1 - ((2.*x)./L - 1).^2.).^(1/2)) - 2./(L.*(1 - ((2.*x)./L - 1).^2).^(1/2)) + (6.*C.*(1 - ((2.*x)./L - 1).^2).^(1/2).*((2.*x)./L - 1))./L))./(2.*sqrt(pi).*(acos(1 - (2.*x)./L) - sin(2.*acos(1 - (2.*x)./L))/2 + C.*(1 - ((2.*x)./L - 1).^2).^(3/2)).^(1/2));
y_f = dy_by_dx(x);

S = 4; %desired segments
d_t_n = S;
f_t_n = 2; %found tangent nodes
g_max = max(y_f)
g_min = min(y_f)
d_range = (g_max-g_min)/(S-1) 
i = 2; %counter
x_t_n = [min(x)]
y_t_n = [0]
y_f_t_n = [g_max]
x_n = [min(x)];
y_n = [0];
x_pr_t_n = min(x);
y_f_pr_t_n = max(y_f);
while d_t_n > f_t_n
    delta = abs(y_f_pr_t_n-y_f(i));
    if delta >d_range
        x_t_n(f_t_n) = x(i);
        y_f_t_n(f_t_n) = y_f(i);
        y_t_n(f_t_n) = ys(i);
        x_pr_t_n = x(i);
        y_f_pr_t_n = y_f(i);
        f_t_n = f_t_n + 1;
    end
    i = i+1;
end
x_t_n(f_t_n) = max(x);
y_f_t_n(f_t_n) = g_min;
y_t_n(f_t_n) = R;

figure 
plot(x,ys,"k",'linewidth',2)
hold on
plot(x,-ys,"k",'linewidth',2)
axis([0 L -0.1 0.1])

for j = 1:(length(x_t_n)-1)
    x1 = x_t_n(j);
    x2 = x_t_n(j+1);
    y1 = y_t_n(j);
    y2 = y_t_n(j+1);
    m1 = y_f_t_n(j);
    m2 = y_f_t_n(j+1);
    x = (y2-y1+m1*x1-m2*x2)/(m1-m2);
    y = m1*(x-x1)+y1;
    x_n(j+1)=x;
    y_n(j+1)=y;
    
end
x_n(length(x_t_n)+1) = L;
y_n(length(y_t_n)+1) = R;

for k=1:length(x_n)-1
    plot([x_n(k);x_n(k+1)], [y_n(k);y_n(k+1)],'r','linewidth',2)
    plot([x_n(k);x_n(k+1)], [-y_n(k);-y_n(k+1)],'r','linewidth',2)
end
hold off
figure
plot(x,y_f,'linewidth',2)
