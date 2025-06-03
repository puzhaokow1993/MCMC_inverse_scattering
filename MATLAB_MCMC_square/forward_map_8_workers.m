function far = forward_map_8_workers(k,thetaout,theta,coefficient_n) 

% thetaout and theta should be vector of length 8 

global R FEM_mesh_size
g = @circleRfunction;
a = @(location,state) coefficient_n(location.x,location.y).*(-k^2);


model1 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model1,g); % Specify geometry from boundary
f1 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,1)))+location.y.*(-sin(theta(1,1))))));
specifyCoefficients(model1, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f1); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model1,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model1,'Hmax',FEM_mesh_size); % generate a mesh ;

model2 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model2,g); % Specify geometry from boundary
f2 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,2)))+location.y.*(-sin(theta(1,2))))));
specifyCoefficients(model2, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f2); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model2,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model2,'Hmax',FEM_mesh_size); % generate a mesh ;

model3 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model3,g); % Specify geometry from boundary
f3 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,3)))+location.y.*(-sin(theta(1,3))))));
specifyCoefficients(model3, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f3); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model3,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model3,'Hmax',FEM_mesh_size); % generate a mesh ;

model4 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model4,g); % Specify geometry from boundary
f4 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,4)))+location.y.*(-sin(theta(1,4))))));
specifyCoefficients(model4, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f4); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model4,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model4,'Hmax',FEM_mesh_size); % generate a mesh ;

model5 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model5,g); % Specify geometry from boundary
f5 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,5)))+location.y.*(-sin(theta(1,5))))));
specifyCoefficients(model5, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f5); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model5,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model5,'Hmax',FEM_mesh_size); % generate a mesh ;

model6 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model6,g); % Specify geometry from boundary
f6 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,6)))+location.y.*(-sin(theta(1,6))))));
specifyCoefficients(model6, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f6); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model6,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model6,'Hmax',FEM_mesh_size); % generate a mesh ;

model7 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model7,g); % Specify geometry from boundary
f7 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,7)))+location.y.*(-sin(theta(1,7))))));
specifyCoefficients(model7, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f7); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model7,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model7,'Hmax',FEM_mesh_size); % generate a mesh ;

model8 = createpde(1); % Create a PDE model with a single equation
geometryFromEdges(model8,g); % Specify geometry from boundary
f8 = @(location,state) complex((k^2).*(coefficient_n(location.x,location.y)-1).*exp(1i*k*(location.x.*(-cos(theta(1,8)))+location.y.*(-sin(theta(1,8))))));
specifyCoefficients(model8, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f8); % specify coefficiens (https://www.mathworks.com/help/pde/ug/pde.coefficientassignment-properties.html)
applyBoundaryCondition(model8,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model8,'Hmax',FEM_mesh_size); % generate a mesh ;

% Use parfeval to run both solvepde tasks in parallel
result1 = parfeval(@solvepde, 1, model1);  % 1st parallel task
result2 = parfeval(@solvepde, 1, model2);  % 2nd parallel task
result3 = parfeval(@solvepde, 1, model3);  % 3rd parallel task
result4 = parfeval(@solvepde, 1, model4);  % 4th parallel task
result5 = parfeval(@solvepde, 1, model5);  % 5th parallel task
result6 = parfeval(@solvepde, 1, model6);  % 6th parallel task
result7 = parfeval(@solvepde, 1, model7);  % 7th parallel task
result8 = parfeval(@solvepde, 1, model8);  % 8th parallel task

fetched_result1 = fetchOutputs(result1); 
fetched_result2 = fetchOutputs(result2); 
fetched_result3 = fetchOutputs(result3); 
fetched_result4 = fetchOutputs(result4); 
fetched_result5 = fetchOutputs(result5); 
fetched_result6 = fetchOutputs(result6); 
fetched_result7 = fetchOutputs(result7); 
fetched_result8 = fetchOutputs(result8); 

uscattered1 = fetched_result1.NodalSolution;
uscattered2 = fetched_result2.NodalSolution;
uscattered3 = fetched_result3.NodalSolution;
uscattered4 = fetched_result4.NodalSolution;
uscattered5 = fetched_result5.NodalSolution;
uscattered6 = fetched_result6.NodalSolution;
uscattered7 = fetched_result7.NodalSolution;
uscattered8 = fetched_result8.NodalSolution;

% compute far-field (estimate by the values at |x| = 3 = 0.75*R)
Rfar = 0.75*R;
far = zeros(1,8); 

XYData = model1.Mesh.Nodes;
farcoef = exp(1i*pi/4 + 1i*k*Rfar)/(sqrt(8*pi*k*Rfar));

xfar1 = Rfar*cos(thetaout(1,1));
yfar1 = Rfar*sin(thetaout(1,1)); 
G1 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered1,'natural');
far(1,1) = G1(xfar1,yfar1)./farcoef;

xfar2 = Rfar*cos(thetaout(1,2));
yfar2 = Rfar*sin(thetaout(1,2)); 
G2 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered2,'natural');
far(1,2) = G2(xfar2,yfar2)./farcoef;

xfar3 = Rfar*cos(thetaout(1,3));
yfar3 = Rfar*sin(thetaout(1,3)); 
G3 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered3,'natural');
far(1,3) = G3(xfar3,yfar3)./farcoef;

xfar4 = Rfar*cos(thetaout(1,4));
yfar4 = Rfar*sin(thetaout(1,4)); 
G4 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered4,'natural');
far(1,4) = G4(xfar4,yfar4)./farcoef;

xfar5 = Rfar*cos(thetaout(1,5));
yfar5 = Rfar*sin(thetaout(1,5)); 
G5 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered5,'natural');
far(1,5) = G5(xfar5,yfar5)./farcoef;

xfar6 = Rfar*cos(thetaout(1,6));
yfar6 = Rfar*sin(thetaout(1,6)); 
G6 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered6,'natural');
far(1,6) = G6(xfar6,yfar6)./farcoef;

xfar7 = Rfar*cos(thetaout(1,7));
yfar7 = Rfar*sin(thetaout(1,7)); 
G7 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered7,'natural');
far(1,7) = G7(xfar7,yfar7)./farcoef;

xfar8 = Rfar*cos(thetaout(1,8));
yfar8 = Rfar*sin(thetaout(1,8)); 
G8 = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered8,'natural');
far(1,8) = G8(xfar8,yfar8)./farcoef;
end