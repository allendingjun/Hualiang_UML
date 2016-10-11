clc;
clear all;
close all;

addpath('E:\Work\HFSS\MatlabHFSS\HFSSApi-source-v0.4.0');

hfssIncludePaths('E:\Work\HFSS\MatlabHFSS\HFSSApi-source-v0.4.0');



tmpPrjFile    = [pwd, '\tmpDipole.hfss'];
tmpDataFile   = [pwd, '\tmpData.m'];
tmpScriptFile = [pwd, '\MS_5.vbs'];

hfssExePath = 'D:\Program Files\AnsysEM\AnsysEM17.1\Win64\ansysedt.exe';

fid=fopen(tmpScriptFile,'wt');

hfssNewProject(fid);
hfssInsertDesign(fid,'MS_5');


% size of the array
% m is the number of layer
% r is the radius of the whole array
% height is post height
% sub_thick is the substrate thickness
% lattice is the lattice constant
m=10;
lattice=450;

height=600;
sub_thick=1500;
r=(m-0.5)*lattice;

% draw substrate
hfssCylinder(fid,'sub','z', [0,0,-sub_thick] ,r,sub_thick,'nm');
% draw posts

hfssCylinder(fid,'Post1_1','z', [0,0,0], radius_post(1,1,m), height, 'nm');

for ii=2:1:m
    for jj=1:1:6*(ii-1)
        hfssCylinder(fid,['Post' num2str(ii) '_' num2str(jj)], 'z', [coordinate_x(ii,jj),coordinate_y(ii,jj),0], radius_post(ii,jj,m), height, 'nm');
    end
end

% draw some additional posts to fill up the whole circular substrate
% additional layer should be less than 2 (for m<10)
for ii=m+1:1:m+2
    for jj=1:1:6*(ii-1)
        if (radius_post(ii,jj,m)+sqrt((coordinate_x(ii,jj))^2+(coordinate_y(ii,jj))^2))<r
            hfssCylinder(fid,['Post' num2str(ii) '_' num2str(jj)], 'z', [coordinate_x(ii,jj),coordinate_y(ii,jj),0], radius_post(ii,jj,m), height, 'nm');
        end
    end
end


% Freq. parameters
fSolve = 3.5294e14;
wSolve = 3e8/fSolve;
wSolve = wSolve*1000;



%hfssRemovePaths('E:\Work\HFSS\MatlabHFSS\HFSSApi-source-v0.4.0');
%rmpath('E:\Work\HFSS\MatlabHFSS\HFSSApi-source-v0.4.0');