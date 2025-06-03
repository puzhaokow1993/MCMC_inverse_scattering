function z = n_true_1(x,y)
z = 1 + (x<1/2).*(x+y>-1/4).*(y<1/2)./10;  
end