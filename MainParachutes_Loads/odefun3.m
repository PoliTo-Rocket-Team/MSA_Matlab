function dydt = odefun3(t,y)
rho = 0.9556;
m_p = 0.213;
m_V = 9.972;
p = 0.2;
CdS = 2.569;
CdS_d = 0.2467;
C_Dr = 0.5;
S_r = 8.495e-03;
CdS_r = C_Dr * S_r;
k_g = 1.068*(1-1.465*p-0.25975*p^2+1.2626*p^3);
k_a = ((4/3)*pi*0.6^3*k_g)/((CdS)^(3/2));
W_V = m_V*9.81;
W_p = m_p*9.81;
t_f = (1.18+1.62+0.76+0.75+0.92+1.03+0.43+0.7+0.72+0.46+0.44+0.6+0.5+1.25+0.72)/15;
t_f = 0.7;
dydt=zeros(4,1);
dydt(4)= 0; %%CdS finale / filling time
dydt(2)=(3/2).*rho.*k_a.*(y(4).^(1/2)).*(0);
dydt(1)=-(y(3)+(0.5.*rho.*y(1).^2.*(CdS_r+CdS_d))-W_V)./m_V;
dydt(3)= 0.5.*rho.*(2.*y(1).*dydt(1).*y(4)+((y(1)).^2).*dydt(4))
+dydt(1).*dydt(2)+0.75.*rho.*k_a.*((y(4)).^(-0.5)).*((dydt(4)).^2).*y(1)
+dydt(2).*dydt(1)-(1/m_V).*(y(2)+m_p).*(dydt(3)+rho.*y(1).*dydt(1).*(CdS_r+CdS_d));

end
