%% Add noise to plotting vertices.
% This can be triangulated and still gives accurate description of the
% amount of noise added to the computational points.
function noisy_texture_map_plot(Vertices, Faces, C, sigma)

N_plot = sigma*[randn(length(Vertices(:,1)),1),...
    randn(length(Vertices(:,2)),1), randn(length(Vertices(:,3)),1)];

V_noise = Vertices + N_plot;

%% Try to do CP_N(u_noise).
% Best idea is to do knnsearch onto plotting vertices.
[IDX] = knnsearch(Vertices,V_noise,'NSMethod','exhaustive');

%% plot triangulation.
figure(3);
patch('Vertices', Vertices, 'Faces', Faces,'FaceVertexCData', C(IDX,:),'FaceColor',...
    'flat','EdgeColor', 'none');
camlight(35,0); view([65,-10]); axis off; axis([-0.9 0.9 -0.9 0.9 0 0.8])
material dull; set(gcf,'PaperSize',[4.5 4]);
