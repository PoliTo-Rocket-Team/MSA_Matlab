function [x, num_it, x_vec] = NewtonMethod(x_0,f,df,tolf,tolr,k_max)
x_vec = zeros(k_max + 1,1);
k = 1;
x_vec(k) = x_0;
err_step = inf;

while k <= k_max && abs(f(x_vec(k))) >= tolf && err_step >= tolr ;
    x_vec(k+1) = x_vec(k) - f(x_vec(k))/df(x_vec(k));
    % aggiorno variabili
    err_step = abs(x_vec(k+1) - x_vec(k))/abs(x_vec(k+1));
    k = k + 1;
end
num_it = k;%k-1
x = x_vec(num_it);
x_vec(num_it + 1:end) = [];

end 