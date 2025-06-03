function [x,y] = circleRfunction(bs,s) 
global R 
if nargin == 0  
    x = 4; % 4 segments
    return 
end
if nargin == 1
    % Outer circle with radius R
    dl = [0      pi/2   pi       3*pi/2
          pi/2    pi     3*pi/2   2*pi
          1       1      1        1 % region label to left (anticlockwise)
          0       0      0        0]; % region label to right (anticlockwise)
    x = dl(:,bs);   
    return 
end 
x = zeros(size(s)); 
y = zeros(size(s)); 
if isscalar(bs) % Does bs need scalar expansion?
    bs = bs*ones(size(s)); % Expand bs
end
cbs = find(bs <= 4); % Outer circle with radius R
x(cbs) = R*cos(s(cbs));
y(cbs) = R*sin(s(cbs));
end