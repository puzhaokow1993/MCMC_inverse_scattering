function far = forward_map(k,thetaout,theta,coefficient_n) 

global R FEM_mesh_size

model = createpde(1); % Create a PDE model with a single equation
g = @circleRfunction;
geometryFromEdges(model,g); % Specify geometry from boundary


a = @(location,state) coefficient_n(location.x,location.y).*(-k^2);
f = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta))+location.y.*(-sin(theta)))));

specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model,'Hmax',FEM_mesh_size); % generate a mesh ;

% Solve the equation
result = solvepde(model);
uscattered = result.NodalSolution;

% Express the incident field exp(1i*k*(x*cos(theta) + y*sin(theta))) by using the mesh
uincident = zeros(length(model.Mesh.Nodes),1);
for n=1:length(model.Mesh.Nodes)
    uincident(n,1) = exp(1i*k*((model.Mesh.Nodes(1,n)*(-cos(theta)) + (model.Mesh.Nodes(2,n)*(-sin(theta))))));
end

% % ######################################################### 
% % ##### the following plots are for debugging purpose ##### 
% % ######################################################### 
% Plot the real part of the incident plane field
% figure
% subplot(2,2,1)
% pdeplot(model,'XYData',real(uincident),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(incident plane wave)'
% axis equal
% 
% % Plot the real part of the scattered field 
% subplot(2,2,2)
% pdeplot(model,'XYData',real(uscattered),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(scattered field)'
% axis equal
% 
% % Plot the imaginary part of the incident plane field
% subplot(2,2,3)
% pdeplot(model,'XYData',imag(uincident),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Im(incident plane wave)'
% axis equal
% 
% % Plot the imaginary part of the scattered field 
% subplot(2,2,4)
% pdeplot(model,'XYData',imag(uscattered),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Im(scattered field)'
% axis equal
% 
% % Plot the mesh 
% figure 
% pdemesh(generateMesh(model)) 
% % #########################################################

% compute far-field (estimate by the values at |x| = 3 = 0.75*R)
Rfar = 0.75*R;

XYData = model.Mesh.Nodes;
G = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered,'natural');
farcoef = exp(1i*pi/4 + 1i*k*Rfar)/(sqrt(8*pi*k*Rfar));
xfar = Rfar*cos(thetaout);
yfar = Rfar*sin(thetaout); 
far = G(xfar,yfar)./farcoef;

end