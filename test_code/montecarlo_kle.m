%% Initialize library
clc;
clearvars;
lib = genpath('../efem_framework');
addpath(lib);
lib = genpath('../ssfem_framework');
addpath(lib);
fprintf('Currently running code: ''%s.m''\n',mfilename);
clear lib;


mcs = MonteCarlo;

mcs.verbose = 0;
mcs.LoadMesh('./models/Multi_Dielectric');
mcs.meshStatistics();
Mat = [1, 9-0.018j, 6-0.012j, 6-0.012j, 6-0.012j, 1.2-0.0024j, 1.2-0.0024];

mcs.setMaterials(Mat);
mcs.buildSystem();

f = 2.5;

kle = KLE3D;
kle.process(mcs,2)
kle.evaluate_Phi('nlambda',3);
kle.getKLEData(1);
kle.KLDATA.sd = 0.8;
mcs.init(100,'kle',true,'p_order',2, 'nklterms', 3);
mcs.setSeed(1);
mcs.assignSpatialMaterialVariation(kle);
mcs.MCSimulation(f);
mcs.plot('current');
