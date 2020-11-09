clc; close all; clear;
addpath(genpath('utils/'));
addpath('func_main/');

mesh_dir = 'data/';
s1_name = 'cat10.off';
s2_name = 'wolf0.off';

%% read the mesh and compute the LB basis
S1 = MESH.MESH_IO.read_shape([mesh_dir, s1_name]);
S2 = MESH.MESH_IO.read_shape([mesh_dir, s2_name]);

meshOptions = {'IfComputeGeoDist',false,'IfComputeLB',true,'IfComputeNormals',true,'numEigs',100};
S1 = MESH.preprocess(S1,meshOptions{:});
S2 = MESH.preprocess(S2,meshOptions{:});

%% set the initial map
B1 = S1.evecs(:,1:4); B2 = S2.evecs(:,1:4);
Ev1 = S1.evals(1:4); Ev2 = S2.evals(1:4); 
numTimes = 100;

fct1 = fMAP.waveKernelSignature(B1, Ev1, S1.A, numTimes);
fct2 = fMAP.waveKernelSignature(B2, Ev2, S2.A, numTimes);

para.type_orient = 'direct';
para.a = 2e-1; para.b = 1e-2;
para.c = 8e-4; para.d = 1e-1;
para.fMap_size = [4, 4];

C12 = fMAP.compute_fMap_from_descriptors(S1, S2, fct1, fct2, para);

figure(1);
for i = 1:4
    subplot(2,4,i); plot_func_on_mesh(S1, B1(:,i)); title([s1_name, ': LB',num2str(i)])
    subplot(2,4,i+4); plot_func_on_mesh(S2, B2(:,i));title([s2_name, ': LB',num2str(i)])
end
% we can see that the first four LB basis of the cat and the wolf align
% with each other, therefore, we can use the identity as the initial
% functionl map
% also: if we start with only the first three LBs, we will get the
% symmetric map
C21_ini = diag([1,1,1,-1]);
T12_ini = fMAP.fMap2pMap(B2,B1,C21_ini);
T12_ini = fMAP.fMap2pMap(B2,B1,C12');

figure(2);
subplot(1,2,1); visualize_map_on_source(S1, S2, T12_ini); title('Source');
subplot(1,2,2); visualize_map_on_target(S1, S2, T12_ini); title('The initial p2p-map w.r.t. the identity fMap')
%% apply zoomOut
para.k_init = 3;
para.k_step = 1;
para.k_final = 100;
tic
[T12, C21, all_T12, all_C21] = zoomOut_refine(S1.evecs, S2.evecs, T12_ini, para);
t = toc;
fprintf('ZoomOut runtime: %.2f sec\n',t)

% fast version
para.num_samples = 200;
tic
T12_fast = zoomOut_refine_fast(S1, S2, T12_ini, para,1);
t = toc;
fprintf('ZoomOut with sampling runtime: %.2f sec\n',t)

%%
figure(3);
plot_id = [2,3,4,5,18,48];
subplot(2,7,1); visualize_map_on_source(S1,S2,T12); title('Source');
for i = 1:length(plot_id)
    subplot(2,7,i+1); imagesc(all_C21{plot_id(i)}); axis square; 
    title(['fMap size: ' num2str(size(all_C21{plot_id(i)},1))]);
    subplot(2,7,i+8); visualize_map_on_target(S1, S2, all_T12{plot_id(i)});
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.34, 1, 0.7]);