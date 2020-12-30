clear all
close all

%We add to Matlab path the place where the code is
addpath('../../Code')

%We load the two unfilter half maps and average them and write the output
%to disk.
vol1 = ReadMRC('emd_10418_half_map_1.map');
vol2 = ReadMRC('emd_10418_half_map_2.map');
vol = 0.5*(vol1+vol2);
clear vol1 vol2;
WriteMRC(vol,0.82,'emd_10418.mrc');

%The threshold value is obtained from Chimera
mask = vol > 0.0153; 
%We remove small obhects, smaller than 25px from the mask
mask = bwareaopen(mask,25,6);

%% LocSpiral
[map W] = locSpiral(vol,mask,0.82,25,2.96,0.9,4.5);
WriteMRC(map,0.82,'locSpiralMap.mrc');

%% local Bfactor sharpening correction 
[M W] = locBSharpen(vol,mask,0.82,25,2.96,0.9,4.5);
WriteMRC(M,0.82,'bfactorSharpenMap.mrc');

%% local Bfactor estimation
[AMap BMap noise Mask Mod resSquare] = locBFactor(vol,mask,0.82,15,2.96,10,0.9,4.8);
WriteMRC(BMap*4,0.82,'bfactorMap.mrc');
WriteMRC(AMap,0.82,'AMap.mrc');

%% local Occupancy
OMap = locOccupancy(vol,mask,0.82,25,10,0.25,5);
WriteMRC(OMap,0.82,'OMap.mrc');
