clc;
clearvars;
lib = genpath('../efem_framework');
addpath(lib);
lib = genpath('../ssfem_framework');
addpath(lib);
clear lib;

ssfem = SSFEM;
ssfem.execute('Multi_Dielectric.sto');