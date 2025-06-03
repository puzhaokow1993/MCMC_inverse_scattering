function z = basis(x,y,r,s,ell_max) 
z = (0<(2^(ell_max))*(x+1)-r).*((2^(ell_max))*(x+1)-r<1).*(0<(2^(ell_max))*(y+1)-s).*((2^(ell_max))*(y+1)-s<1); 
end