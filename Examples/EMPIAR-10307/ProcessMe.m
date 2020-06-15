clear all
close all

%We add to Matlab path the place where the code is
addpath('../../Code')

%We load the two unfilter half maps and average them and write the output
%to disk.
vol1 = ReadMRC('relion_half1_class001_unfil.mrc');
vol2 = ReadMRC('relion_half2_class001_unfil.mrc');
vol = 0.5*(vol1+vol2);
clear vol1 vol2;
WriteMRC(vol,1.699,'class001.mrc');

%The threshold value is obtained from Chimera
mask = vol > 0.014;

%We remove small obhects, smaller than 25px from the mask
mask = bwareaopen(mask,25,6);

%% LocSpiral. For visualization it is recommended in Chimera threshold 1.89
%  and hide dust 100 
[map W] = locSpiral(vol,mask,1.699,30,4.25,0.95,9);
WriteMRC(map,1.699,'locSpiralMap.mrc');

%% local Bfactor sharpening correction 
[map W] = locBSharpen(vol,mask,1.699,30,4.25,0.9,6);
WriteMRC(map,1.699,'bfactorSharpenMap.mrc');

%% local Bfactor estimation - 1
[AMap BMap noise Mod resSquare] = locBFactor(vol,mask,1.699,15,4.25,10,0.9,6);
WriteMRC(BMap*4,1.699,'bfactorMap-1.mrc');
WriteMRC(AMap,1.699,'AMap-1.mrc');

[AMap BMap noise Mod resSquare] = locBFactor(vol,mask,1.699,20,10,10,0.9,6);
WriteMRC(BMap*4,1.699,'bfactorMap-2.mrc');
WriteMRC(AMap,1.699,'AMap-2.mrc');

%% local Occupancy
OMap = locOccupancy(vol,mask,1.699,25,10,0.25,6);
WriteMRC(OMap,1.699,'OMap.mrc');
