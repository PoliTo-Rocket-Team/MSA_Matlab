%Nose cone function and tangent approximation
clc
close all
clear all
format long

L = 0.45; %nose length
R = 0.05 ; %nose base radius
A_ref = pi*R^2;
C= 0; %haack series parameter, 0 for Von Karman
n = 1/2; %power series parameter
c = 340;
M_inf = 1.2;
V_inf = c*M_inf;
Machs = [ 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.2 2.4 2.6];
%Machs = [1.8];
DatMachs = [1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6];
DatCds = [0.59 0.571 0.558 0.547 0.534 0.523 0.512 0.501 0.491 0.481 0.465 0.457 0.448 0.440 0.433];
Cdps = zeros(length(Machs),1);
g = 1.4; %ratio of specific heats
P_inf = 101325;%freestream pressure, atmosferic pressure
P_inf_bar = P_inf/(100000);
T_inf = 273; %freestream temperature
rho_inf = 1.225; %freestream density

x = linspace(0,L,1000);
%Ps =zeros(1000,1);

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

S = 50; % Desidered Segments
N = S + 1; % Number of Nodes
delta = (max(y_f)-min(y_f))/(N); %gradient delta
x_t = zeros(N,1); %x position of the tangent nodes
y_t = zeros(N,1); %y position of the tangent nodes
x_t(1) = 0; %the tip is set to be a tangent node
y_t(1) = 0;
k = 2;%number of tangent nodes alreay found (tip+base)

Mcs = ones(S,1);
nu2s = ones(S,1);
M2s = ones(S,1);
Pcs = ones(S,1);
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
    x1 = x_t(j);
    x2 = x_t(j+1);
    y1 = y_t(j);
    y2 = y_t(j+1);
    m1 = n*(R/L^n)*x1^(n-1);%slopes of the tangents in the tangent nodes
    m2 = dy_by_dx(x2);
    x = (y2-y1+m1*x1-m2*x2)/(m1-m2);
    y = m1*(x-x1)+y1;
    x_n(j+1)=x;
    y_n(j+1)=y;
    
end
x_n(length(x_t)+1) = L;
y_n(length(y_t)+1) = R;


%barrowman method + expansion waves 
NU = @(M) (atan((M^2 - 1)^(1/2)*((g - 1)/(g + 1))^(1/2))*(g + 1)^(1/2))/(g - 1)^(1/2) - atan((M^2 - 1)^(1/2));
dNU = @(M) (M*((g - 1)/(g + 1))^(1/2)*(g + 1)^(1/2))/((M^2 - 1)^(1/2)*(((M^2 - 1)*(g - 1))/(g + 1) + 1)*(g - 1)^(1/2)) - 1/(M*(M^2 - 1)^(1/2));
B = @(M) (g*M^2)/(2*(M^2-1));
Omega = @(M) (1/M)*((1+M^2*(g-1)/2)/((g+1)/2))^((g+1)/(2*(g-1)));
for l = 1:length(Machs)
M_inf = Machs(l);
%M_inf = 2;
V_inf = c*M_inf;
Cdp = 0;
Cdwls = zeros(length(x_t)-1);
for k=1:length(x_t)-1
    theta = atan((y_t(k+1)-y_t(k))./(x_t(k+1)-x_t(k)))
    theta_deg = theta *(180/pi)
    if k ==1
        [M0,P0_bar] = taylor_maccoll_solver(M_inf, theta_deg, g, P_inf_bar, T_inf, rho_inf)%P is constant
        P0 = 100000*P0_bar;
        Fx = @(x) P0*2*pi*sin(theta)*sin(theta)*x;   
        %Cdwl = (4/(g*M_inf^2))*int(P0*sin(theta)*2*pi,x,0,x_t(k+1)) given
        %by Barrowman
        Cdwl = (integral(Fx,0,x_t(k+1)))/(0.5*rho_inf*V_inf^2*A_ref) %local Cd, a part of the final integral
        if M0 <=1.1
            M0 = 1.1
        end
        M1 = M0; %mach number on the first frustrum is constant
        M2s(1) = M0;
        nu2s(1) = NU(M0);
        Mcs(1) = M0;
        Pcs(1)=P0;
        dpds1 = 0; 
        P2 = P0;
        xl = linspace(x_t(k),x_t(k+1),100);
        Ps(1:100)=P0;
        xps(1:100)=xl;
    else
        [MC,PC_bar] = taylor_maccoll_solver(M_inf, theta_deg,g, P_inf_bar, T_inf, rho_inf)%to find the equivalent pressure on a cone with the same half-angle
        PC = 100000*PC_bar;
        B1 = B(M1);
        nu1 = real(NU(M1));
        O1 = Omega(M1);
        f = @(M) (atan((M^2 - 1)^(1/2)*((g - 1)/(g + 1))^(1/2))*(g + 1)^(1/2))/(g - 1)^(1/2)- atan((M^2 - 1)^(1/2)) - (thetapr-theta) - nu1;
        df = @(M) (M*((g - 1)/(g + 1))^(1/2)*(g + 1)^(1/2))/((M^2 - 1)^(1/2)*(((M^2 - 1)*(g - 1))/(g + 1) + 1)*(g - 1)^(1/2)) - 1/(M*(M^2 - 1)^(1/2));
        [M2,num_it,M2_vec] = NewtonMethod(M_inf,f,df,-10,0.01,10);%getting the mack behind the shockwave using the newton method 
        M2 = real(M2)
        if M2<=1.1
            M2 = 1.1;
        end
        B2 = B(M2);
        nu2 = abs(thetapr-theta) + nu1;
        O2 = Omega(M2);
        M2s(k) = M2;
        nu2s(k) = nu2;
        Mcs(k) = MC;
        Pcs(k) = PC;
        dpds2 = (B2/y_t(k))*((O1/O2)*sin(thetapr)-sin(theta))+(B2/B1)*(O1/O2)*dpds1;
        %eta = dpds2 * (x-x2)/((PC-P2)*cos(theta))
        %P = @(x) (PC-(PC-P2)*exp(-dpds2 * (x-x_t(k))/((PC-P2)*cos(theta))));%pressure distribution on the frustrum
        Pe = P2*((1+(M1^2)*((g-1)/2))/(1+(M2^2)*((g-1)/2)))^(g/(g-1))
        %Fx = @(x) 2*pi.*sin(theta).*sin(theta).*x.*(PC-(PC-P2).*exp(-dpds2 .* (x-x_t(k))./((PC-P2).*cos(theta))));%drag force on the frustrum
        Fx = @(x) 2*pi.*sin(theta).*sin(theta).*x.*Pe;
        Cdwl = (integral(Fx,x_t(k),x_t(k+1)))/(0.5*rho_inf*V_inf^2*A_ref)%local cd component
        %updating the variables
        M1 = M2;
        dpds1 = dpds2;
        %P2 = P(x_t(k+1))
        P2 = Pe;
        xl = (linspace(x_t(k),x_t(k+1),100));
        %Pl = P(xl);
        Pl = transpose(ones(length(xl),1)*Pe);
        xps = cat(2,xps,xl);
        Ps = cat(2,Ps,Pl);
        
    end
    plot([x_n(k);x_n(k+1)], [y_n(k);y_n(k+1)],'g','linewidth',2)
    plot([x_n(k);x_n(k+1)], [-y_n(k);-y_n(k+1)],'g','linewidth',2)
    thetapr = theta;
    Cdp = Cdp + Cdwl
    Cdwls(k) = Cdwl;
end
Cdps(l)=Cdp;
Contribute = (Cdwls./Cdp).*100;

figure
plot(1:S,Mcs,'linewidth',2)
title(M_inf)
xlabel('frustrum n.')
ylabel('equivalent mach')

figure
plot(1:S,M2s,'linewidth',2)
title(M_inf)
xlabel('frustrum n.')
ylabel('mach at the base (M2)')

figure
plot(1:S,nu2s,'linewidth',2)
title(M_inf)
xlabel('frustrum n.')
ylabel('nu function')

figure
plot(xps,Ps,'linewidth',2)
title(M_inf)
xlabel('x')
ylabel('Pressure (Pa)')

figure
plot(1:S,Pcs,'linewidth',2)
title(M_inf)
xlabel('x')
ylabel('equivalent pressure (Pa)')

Ps = ones(1,100);
xps = ones(1,100);
end
figure
plot(Machs,Cdps,'linewidth',2)
hold on
plot(DatMachs, DatCds,'linewidth',2)
legend('Barrowman + Taylor-Maccoll (Nose cone only)','Missile Datcom (whole rocket)')
xlabel('Mach number')
ylabel('Pressure drag coefficient')
hold off
figure
plot(1:(length(x_t)-1),Contribute,'linewidth',2)
xlabel('Frustrum n.')
ylabel('Cd contribution (Percentage)')

% figure
% plot(xps,Ps,'linewidth',2)
% xlabel('x')
% ylabel('Pressure (Pa)')
