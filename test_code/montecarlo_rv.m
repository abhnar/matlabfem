%% Initialize library
clc;
clearvars;
lib = genpath('../efem_framework');
addpath(lib);
lib = genpath('../ssfem_framework');
addpath(lib);
fprintf('Currently running code: ''%s.m''\n',mfilename);
clear lib;
%%

mcs = MonteCarlo;

mcs.verbose = 0;
mcs.LoadMesh('./models/Multi_Dielectric');
mcs.meshStatistics();
Mat = [1, 9-0.018j, 6-0.012j, 6-0.012j, 6-0.012j, 1.2-0.0024j, 1.2-0.0024];

mcs.setMaterials(Mat);
freq = linspace(2,4,100);
f = 2.5;
mcs.init(100,'kle',false,'p_order',2);

mcs.setSeed(1);
mcs.assignRandomMaterialVariation(2, 0.8)
mcs.MCSimulation(f);

mcs.plot('current')
