function dydt = odefun(t,y,C_Dr_in,p_in,t_f_in)
rho = 0.9556;
m_p = 0.213;
m_V = 9.972;
p = p_in;
CdS = 2.569;
CdS_d = 0.2467;
C_Dr = C_Dr_in;
S_r = 8.495e-03;
CdS_r = C_Dr * S_r;
k_g = 1.068*(1-1.465*p-0.25975*p^2+1.2626*p^3);
k_a = ((4/3)*pi*0.6^3*k_g)/((CdS)^(3/2));
W_V = m_V*9.81;
W_p = m_p*9.81;
t_f = t_f_in;
dydt=zeros(4,1);
dydt(4)=CdS/t_f; %%CdS finale / filling time
dydt(2)=(3/2).*rho.*k_a.*(y(4).^(1/2)).*(CdS/t_f);
dydt(1)=-(y(3)+(0.5.*rho.*y(1).^2.*(CdS_r+CdS_d))-W_V)./m_V;
dydt(3)= 0.5.*rho.*(2.*y(1).*dydt(1).*y(4)+((y(1)).^2).*dydt(4))
+dydt(1).*dydt(2)+0.75.*rho.*k_a.*((y(4)).^(-0.5)).*((dydt(4)).^2).*y(1)
+dydt(2).*dydt(1)-(1/m_V).*(y(2)+m_p).*(dydt(3)+rho.*y(1).*dydt(1).*(CdS_r+CdS_d));

end

% PROVA