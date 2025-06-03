close all
clear all
clc 

tic

% meaningful parameters 
k = 5; % wave number
N = 560; % sample size (we choose multiply of 8)
sigma = 0.1; 

% technical parameters 
global R FEM_mesh_size
FEM_mesh_size = 0.3;
R = 4;
beta = 0.05; % learning rate 

% create mesh and extract the nodes 
model = createpde(1); % Create a PDE model with a single equation
g = @circleRfunction;
geometryFromEdges(model,g); % Specify geometry from boundary
mesh = generateMesh(model,'Hmax',FEM_mesh_size); % generate a mesh ;
XY_coordinates = mesh.Nodes; % Nodes is a 2-by-N array consists of the coordinates of the mesh nodes 
x = XY_coordinates(1,:); 
y = XY_coordinates(2,:); 
% plot(nodes(1,:),nodes(2,:),'.') % Plot of nodes 
clear model g mesh ans 

% simulating the sampling experiment 
angles = rand(2,N); % first (resp. second) row consists of random chosen N incident (resp. outgoing) angle 
measurements = zeros(1,N); 
for h=1:(N/8) 
    measurements(1,8*h-7:8*h) = forward_map_8_workers(k,angles(2,8*h-7:8*h),angles(1,8*h-7:8*h),@n_true_1); % simulating the measurements Y_i in the manuscript 
end 

% initial guess 
initial_F0 = zeros(1,length(XY_coordinates(1,:))); 
iter = 10000; % number of iterations 
ell_max = 1; % truncation level

array_log_likelihood = zeros(1,iter);
F_seq = cell(1,iter+1); % create an empty cell 
F_seq{1} = initial_F0; 
keys = cell(1,iter+1); % create an empty cell 
keys{1} = 'accept'; 

for tau=1:iter 
    F_initial = F_seq{tau}; 
    F_seq_next = sqrt(1-beta^2)*F_initial; 
    parfor t=1:(16^(ell_max))
        r = rem(t-1,4^(ell_max)); 
        s = 4^(-ell_max)*(t-1-r); 
        F_seq_next = F_seq_next + beta*randn*basis(x,y,r,s,ell_max); 
    end
    if (keys{tau} == 'accept') 
        current_log_likelihood = log_likelihood_8_workers(k,N,XY_coordinates,measurements,F_initial,angles,sigma); 
    end
    array_log_likelihood(1,tau) = current_log_likelihood; 
    alpha = log_likelihood_8_workers(k,N,XY_coordinates,measurements,F_seq_next,angles,sigma)-current_log_likelihood; 
    u = log(rand); 
    if (u < alpha) 
        F_seq{tau+1} = F_seq_next; % accept the proposal F_seq_next  
        keys{tau+1} = 'accept'; 
    else
        F_seq{tau+1} = F_seq{tau}; % reject the proposal F_seq_next  
        keys{tau+1} = 'reject'; 
    end
    disp(['Iteration: ', num2str(tau),'/', num2str(iter)]) % show the number of current iteration 
    disp(['current log likelihood: ', num2str(current_log_likelihood)]) % show the number of current iteration 
    if mod(tau,5000) == 0 
        save('backup') 
    end
    clear F_seq_next F_tilde F_initial ell_next 
end

key = 'accept'; % key to extract 
idx = strcmp(keys,key); % find indices of matching keys 
F_accepted = F_seq(idx); % extract corresponding values 

% Define the burn-in period (number of samples to discard) 
burnIn = floor(length(F_accepted)/2); 
F_accepted_remain = F_accepted(burnIn+1:end); % extract the values after the burn-in 

sum_F = F_accepted_remain{1}; 
for h=2:length(F_accepted_remain) 
    sum_F = sum_F + F_accepted_remain{h}; 
end
burn_in_mean_F = sum_F/length(F_accepted_remain); 
elapsedtime = toc; 

% save all workspace variables in the file experiment1.mat 
save('experiment1') 
% load('experiment1') %% load all variables from the file experiment1.mat 


% Create a 3D scatter plot
figure;
scatter3(XY_coordinates(1,:), XY_coordinates(2,:), link_func(burn_in_mean_F), 50 ,link_func(burn_in_mean_F), 'filled'); % marker size = 50 
% Apply a colormap to represent Z values
colormap jet;  % Use the 'jet' colormap (can use 'parula', 'cool', etc.)
colorbar;  % Display a colorbar to indicate the Z values
title(['N = ' num2str(N) ' Elapsed time: ' num2str(elapsedtime) 'seconds']);
xlabel('x');
ylabel('y');
axis([-1 1 -1 1 0.95 1.13])
clim([0.95 1.13])

% Create a 3D scatter plot
figure;
n_true_plot = n_true_1(XY_coordinates(1,:), XY_coordinates(2,:)); 
scatter3(XY_coordinates(1,:), XY_coordinates(2,:),n_true_plot ,50 ,n_true_plot, 'filled'); % marker size = 50 
% Apply a colormap to represent Z values
colormap jet;  % Use the 'jet' colormap (can use 'parula', 'cool', etc.)
colorbar;  % Display a colorbar to indicate the Z values
title('true n'); 
xlabel('x');
ylabel('y');
axis([-1 1 -1 1 0.95 1.13])
clim([0.95 1.13])

% Plot log-likelihood of each iterations (log scale on y axis)  
figure; 
iter_array = 1:iter; 
semilogy(iter_array,array_log_likelihood,'-o')

% Plot log-likelihood of each iterations (log scale on x axis)  
figure; 
iter_array = 1:iter; 
semilogx(iter_array,array_log_likelihood,'-o')