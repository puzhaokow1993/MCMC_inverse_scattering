function ell = log_likelihood_8_workers(k,N,XY_coordinates,measurements,coefficient_F,angles,sigma)  
Gf = zeros(1,N); 
for h=1:(N/8)
    Gf(1,8*h-7:8*h) = param_forward_map_8_workers(k,angles(2,8*h-7:8*h),angles(1,8*h-7:8*h),coefficient_F,XY_coordinates); 
end
ell = -0.5*((norm(measurements-Gf)/sigma)^2); 
end 