%% Initialize library
clc;
clearvars;
lib = genpath('../ssfem_framework');
addpath(lib);
fprintf('Currently running code: ''%s.m''\n',mfilename);
clear lib;


ssfem = SSFEM;

ssfem.verbose = 0;
ssfem.LoadMesh('./models/Multi_Dielectric');
ssfem.meshStatistics();
Mat = [1, 9-0.018j, 6-0.012j, 6-0.012j, 6-0.012j, 1.2-0.0024j, 1.2-0.0024];

ssfem.setMaterials(Mat);
ssfem.buildSystem();

f = 2.5;

%% KLE Setup
kle = KLE3D;
kle.process(ssfem,2)
kle.evaluate_Phi('nlambda',10);
kle.getKLEData(1);

%% Spectral Stochastic
ssfem.init(100,'kle',true,'p_order',2, 'nklterms', 3);
kle.KLDATA.sd = 0.8;
ssfem.setSeed(2);
ssfem.assignSpatialMaterialVariation(kle);
ssfem.ssfemkle(f);
ssfem.plot('current');