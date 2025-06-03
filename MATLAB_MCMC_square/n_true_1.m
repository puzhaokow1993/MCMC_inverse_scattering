function z = n_true_1(x,y)
z = 1 + (x>-1/2).*(x<1/2).*(y>-1/2).*(y<1/2)./10;  
end