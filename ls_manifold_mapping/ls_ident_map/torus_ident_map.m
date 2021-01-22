%% Computation of identity map for a torus from an initial noisy map.
% The code here is for the level set method and is set up for a convergence
% study. For smaller dx values the Compute Canada WestGrid servers had to 
% be used.
%
% Note that due to use of random noisy map code must be ran multiple times
% to obtain the convergence rates. 
%
% For more information on the level set method see Facundo Memoli,
% Guillermo Sapiro and Stanley Osher, "Solving variational problems and 
% partial differential equations mapping into general target manifolds",
% Journal of Computational Physics, 2004.
restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% set a list of grid sizes for convergence study
dx_max = 0.2;
imax = 4;           % convergence study up to imax 5 done on WestGrid
dx = dx_max*2.^(-(1:imax)+1);     % half each grid size

%% banding: do calculation in a narrow band around the torus
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% repeat convergence study for num_rep trials
num_rep = 1;           % number of random trials
err = zeros(length(dx),num_rep);
t_scheme = zeros(length(dx),num_rep);
tExtension = zeros(length(dx),num_rep);
for j = 1:num_rep
    for i = 1:length(dx)    
        %% Construct a grid in the embedding space
        % make vectors of x, y, positions of the grid
        x1d = (-6.0:dx(i):6.0)';
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
        band = find(abs(dist) <= 3*bw*dx(i));
        errband = find(abs(dist(band)) <= 2*dx(i));
        % store closest points in the band
        cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
        % store the signed distance function in the band
        sd = dist(band);
        % store locations where signed distance function is given in band
        xg = xx(band); yg = yy(band); zg = zz(band);

        %% Create upwind schemes
        disp('Constructing upwind schemes in matrix form')
        [Dxb,Dxf, Dyb,Dyf, Dzb,Dzf] = firstderiv_upw1_3d_matrices(x1d,y1d,z1d, band);
        
        %% Create extension matrix
        E = interp3_matrix(x1d,y1d,z1d, cpxg, cpyg, cpzg, p, band);
        
        %% Create projection matrices
        % note: only projection matrix for M can be constructed before time
        % stepping loop. Projection onto N must be updated each time step.
        disp('Constructing center difference schemes')
        [Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d,y1d,z1d, band);
        
        grad_sd = [Dxc*sd, Dyc*sd, Dzc*sd];
        
        rep_rows = kron(grad_sd,ones(3,1));
        T = grad_sd';
        rep_cols = kron(T(:),ones(1,3));
        big_eye = repmat(eye(3),length(cpxg),1);
        proj_M1 = sparse(big_eye - rep_rows.*rep_cols);

        d1 = 3*ones(length(cpxg),1);
        proj_M2 = mat2cell(proj_M1,d1,3);
        proj_M = blkdiag(proj_M2{:});
        
        % build interpolation functions to construct projection onto N
        Fx = scatteredInterpolant(xg, yg, zg, sd, 'linear');
       
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
        
        proj_u1 = zeros(length(cpxg),3);
        proj_u2 = zeros(length(cpxg),3);
        proj_u3 = zeros(length(cpxg),3);
        proj_MN = zeros(length(cpxg),3);
        proj_N = zeros(3,3,length(cpxg));
        tExt = zeros(floor(numtimesteps/5),1);
        
        tStart = tic;
        ExtCount = 0;
        for kt = 1:numtimesteps
            
            if mod(kt,5) == 0
               tExtensionStart = tic;
               % perform "extension evolution"
               u1 = E*u1;
               u2 = E*u2;
               u3 = E*u3;
               ExtCount = ExtCount + 1;
               tExt(ExtCount) = toc(tExtensionStart);
            end
            
            % compute gradient of map
            grad_u1 = [Dxf*u1, Dyf*u1, Dzf*u1]';
            grad_u2 = [Dxf*u2, Dyf*u2, Dzf*u2]';
            grad_u3 = [Dxf*u3, Dyf*u3, Dzf*u3]';
            
            % project onto torus M
            proj_u1 = proj_M*grad_u1(:);
            proj_u2 = proj_M*grad_u2(:);
            proj_u3 = proj_M*grad_u3(:);
 
            proj_u1 = reshape(proj_u1,3,length(cpxg))';
            proj_u2 = reshape(proj_u2,3,length(cpxg))';
            proj_u3 = reshape(proj_u3,3,length(cpxg))';
               
            % compute divergence of gradient of map
            div1 = Dxb*proj_u1(:,1) + Dyb*proj_u1(:,2) + Dzb*proj_u1(:,3);
            div2 = Dxb*proj_u2(:,1) + Dyb*proj_u2(:,2) + Dzb*proj_u2(:,3);
            div3 = Dxb*proj_u3(:,1) + Dyb*proj_u3(:,2) + Dzb*proj_u3(:,3);
            
            % project onto torus N
%            [~,~,~,sd_new] = cpTorus(u1, u2, u3);
            sd_new = Fx(u1, u2, u3);
            grad_sd_new = [Dxc*sd_new, Dyc*sd_new, Dzc*sd_new];
            
            rep_rows = kron(grad_sd_new,ones(3,1));
            T = grad_sd_new';
            rep_cols = kron(T(:),ones(1,3));
            proj_N = big_eye - rep_rows.*rep_cols;

            div = [div1, div2, div3];
            div_rep_rows = kron(div,ones(3,1));
            proj_MN1 = div_rep_rows.*proj_N;
            proj_MN = sum(proj_MN1,2);
            proj_MN = reshape(proj_MN,3,length(cpxg))';
            
            % explicit Euler timestepping
            u1 = u1 + dt*proj_MN(:,1);
            u2 = u2 + dt*proj_MN(:,2);
            u3 = u3 + dt*proj_MN(:,3);

            [u1, u2, u3] = cpTorus(u1, u2, u3);
            t = kt*dt
        end
        tExtension(i,j) = sum(tExt)/ExtCount;
        ExtCount
        t_scheme(i,j) = toc(tStart);
        % average 2-norm of the difference between point locations
        point_diff = [u1, u2, u3] - [cpxg, cpyg, cpzg];
        %err(i,j) = sum(sqrt(sum(abs(point_diff).^2,2)))/length(point_diff);
        err(i,j) = max(sqrt(sum(abs(point_diff(errband)).^2,2)));    
    end
end
tfin = toc;

% compute order of convergence and average errors
order = log2(err(1:end-1,:)./err(2:end,:));
order_avg = sum(order,2)/num_rep
err_avg = sum(err,2)/num_rep
t_scheme_avg = sum(t_scheme,2)/num_rep
tExtension_avg = sum(tExtension,2)/num_rep
percentExtTime = (tExtension_avg./t_scheme_avg)*100