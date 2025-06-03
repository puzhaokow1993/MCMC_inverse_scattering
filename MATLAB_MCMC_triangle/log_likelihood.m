function ell = log_likelihood(k,N,XY_coordinates,measurements,coefficient_F,angles,sigma)  
Gf = zeros(1,N); 
for h=1:N
    Gf(1,h) = param_forward_map(k,angles(2,h),angles(1,h),coefficient_F,XY_coordinates); 
end
ell = -0.5*((norm(measurements-Gf)/sigma)^2); 
end 