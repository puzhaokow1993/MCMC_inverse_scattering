function phi = link_func(t) 
M_0 = 1; 
t_1 = 8; 
t_2 = -8; 
a = 2; 

C_0 = (M_0 - t_1.^(-a+1)-(-t_2).^(-a+1))./(t_1-t_2); 
C_1 = (-t_2).^(-a+1) - t_2.*(M_0 - t_1.^(-a+1)-(-t_2).^(-a+1))./(t_1-t_2); 

% small perturbation here can avoid NaN issue 
phi = ((t > t_1).*(M_0 - (t+10^(-16)).^(-a+1)) + (t < t_2).*((-t-10^(-16)).^(-a+1)) + (t_2 <= t).*(t <= t_1).*(C_0.*t + C_1))./C_1; 
end 