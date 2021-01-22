%% Computation of the identity map for a torus from an initial noisy map.
% The code here is set up for a convergence study. For smaller dx values
% the Compute Canada WestGrid servers had to be used.
%
% Note that due to use of random noisy map code must be ran multiple times
% to obtain the convergence rates. 
restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Set a list of grid sizes for convergence study
dx_max = 0.2;
imax = 4;           % convergence study up to imax 5 done on WestGrid
dx = dx_max*2.^(-(1:imax)+1);     % half each grid size

%% Banding: do calculation in a narrow band around the torus
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Repeat convergence study for num_rep trials
num_rep = 1;           % number of random trials
err = zeros(length(dx),num_rep);
t_scheme = zeros(length(dx),num_rep);
for j = 1:num_rep
    for i = 1:length(dx)    
        %% Construct a grid in the embedding space
        % make vectors of x, y, positions of the grid
        x1d = (-2.0:dx(i):2.0)';
        y1d = x1d;
        z1d = x1d;

        %% Find closest points on the surface
        % For each point (x,y,z), store closest point (cpx,cpy,cpz)
        % meshgrid is only needed for finding the closest points
        [xx, yy, zz] = meshgrid(x1d, y1d, z1d);
        % function cpTorus for finding the closest points on torus
        [cpx, cpy, cpz, dist] = cpTorus(xx,yy,zz);
        % make into vectors
        cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

        % compute band
        band = find(abs(dist) <= bw*dx(i));
        % store closest points in the band
        cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);

        %% Create Laplacian and extension matrix for heat equation
        disp('Constructing Laplacian and extension matrix')
        E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);
        L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

        %% Add noise to closest points.
        sigma = 0.05;
        rng('shuffle')
        N = sigma*[randn(length(cpxg),1),...
                   randn(length(cpyg),1),...
                   randn(length(cpzg),1)];
        Ncp = [cpxg, cpyg, cpzg] + N;

        % map noisy data back onto torus and assign initial u1, u2, u3
        [u1, u2, u3] = cpTorus(Ncp(:,1), Ncp(:,2), Ncp(:,3));

        %% Time evolution for 2-harmonic map.
        Tf = 0.01;
        dt = 0.1*dx(i)^2;
        numtimesteps = ceil(Tf/dt);
        % adjust for integer number of steps
        dt = Tf / numtimesteps;
        tStart = tic;
        for kt = 1:numtimesteps
            % explicit Euler timestepping
            utilde1 = u1 + dt*L*u1;
            utilde2 = u2 + dt*L*u2;
            utilde3 = u3 + dt*L*u3;
    
            unew1 = E*utilde1;
            unew2 = E*utilde2;
            unew3 = E*utilde3;

            [u1, u2, u3] = cpTorus(unew1, unew2, unew3);
            
            t = kt*dt
        end
        t_scheme(i,j) = toc(tStart);
        % average 2-norm of the difference between point locations
        point_diff = [u1, u2, u3] - [cpxg, cpyg, cpzg];
%         err(i,j) = sum(sqrt(sum(abs(point_diff).^2,2)))/length(point_diff);
        err(i,j) = max(sqrt(sum(abs(point_diff).^2,2)));      
    end
end
tfin = toc;

% compute order of convergence and average errors
order = log2(err(1:end-1,:)./err(2:end,:));
order_avg = sum(order,2)/num_rep
err_avg = sum(err,2)/num_rep
t_scheme_avg = sum(t_scheme,2)/num_rep
