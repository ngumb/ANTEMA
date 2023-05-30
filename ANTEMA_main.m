% Main script for Segmentation Analysis of HRTEM images
% Run this script to run the automated nanoparticle transmission electron
% micrograph analyzer (ANTEMA)


% Copyright (C) 2022  University of Duisburg-Essen, 
%                       Nina Gumbiowski
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


clear
close all
addpath('SegmentationFunctions','Nets','dm3read')

%% Settings

netname = 'NP-net-weights.mat';
Convexthresh = 0.95;    % threshold for particle separation based on convexity
PSep = true;            % Particle separation enabled (true) or disabled (false)
saveas = 'xlsx';        % set to 'csv' or 'xlsx'
removecount = 0.5;      % minimum particle size in nm
minmarker = 0.5;        % minimum marker size for particle separtion routine in nm
% Colormap for the raw segmentation map
cmap = [
         0 114 189      % Background
         237 177 32     % Particle 
         ];
     cmap = cmap ./ 255;
%% load network for segmentation
load(netname, 'net')

%% Initialize as batch or single file
batch=input('Segment Batch or single file?\n Press:\n 0 : single file\n 1 : batch\n');
if batch==0 % single file option
    fprintf('Choose dm3 file to be segmented \n')
    [baseName, pathdir] = uigetfile('*.dm3');
    filename = {fullfile(pathdir, baseName)};
       
elseif batch==1 % batch option for all images in one folder
    fprintf('Choose folder that contains all .dm3 files to be segmented \n')
    pathdir=uigetdir();
    [stat,scan]=fileattrib([pathdir '\*.dm3']);
    filename={scan.Name};
    fprintf('%d files found in specified folder\n',length(scan)) 
    allinone = input('How should the data be saved?\n 0: one file for each image\n 1: all image data in one file\n 2: both\n'); %%IMPLEMENT SAVING OPTION AT THE END
    if batch == 1  && allinone > 0
        AllpartProp = [];
    end
else
    fprintf('ERROR: Irregular input for Segment Batch or single file\n')
    return
end

%% start timer for performance evaluation
tic

%% Otion to save elsewhere
store = input('Where should the results be stored?\n 1: Results folder (default)\n 2: Directory of the origin data\n 3: Other \n');
if store == 1
    resfolder = 'Results';
elseif store == 2
    resfolder = pathdir;
elseif store == 3
    resfolder = uigetdir();
else
    fprintf('ERROR: Irregular input \n')
    return
end
%% Processing images
for i=1:length(filename)
    close all
    % load singular dm3 file 
    filepath=filename{i};
    [pathstr,expName,ext] = fileparts(filepath); 
    [Image,pxsz,units]=ReadDMFile(filepath);
    if isempty(units)
        warning('Unit of dm3 file is empty, Skip evaluation of %s\n',filepath)
        continue
    elseif units~='nm'
        warning('Unit not in nm, Case currently not implemented\n Skip evaluation of %s\n',filepath)
        continue
    end
    % apply filter
    Image=imgaussfilt(Image,2);
    % convert to uint8
    Image = uint8(255 * mat2gray(Image));
    % Segmentation
    tSegStart=toc
    C = semanticseg(Image, net);
    tSegEnd=toc
    fprintf('Segmentation of file No. %i completed\n',i)
    C =uint8(C);
    
    % Show image
%     figure
%     CB=labeloverlay(Image,C);
%     imshowpair(Image,CB,'montage')
%  
    %% get properties
    [Eval,partProp,B,Mlines] = particlePropertiesEval(C,2,pxsz,Convexthresh,PSep,Image,minmarker,removecount);

    if isempty(partProp)
        warning('Property measurement is empty\n No data generated\n')
        
        continue
    end

    % Make image Overlays
    B2=labeloverlay(Image,C,'Colormap',cmap);
    

    %% Add Scalebar to Image
    Image=scalebar(Image,pxsz,units);

    %% Save Results
    % make directory for results
    mkdir(resfolder,expName) 
    
    % save overlay image and raw image     
    imwrite(B,[resfolder '\' expName '\' expName '_ImageOverlay.png'])
    imwrite(B2,[resfolder '\' expName '\' expName '_ImageOverlayRaw.png'])
    imwrite(Image,[resfolder '\' expName '\' expName '_Image.png'])

    if PSep == true
        Mlines = imdilate(Mlines,strel('square',3));
         Mlines=labeloverlay(Image,Mlines,'Colormap',[1,1,1]);
        imwrite(Mlines,[resfolder '\' expName '\' expName 'ImageOverlayLines.png'])
    end

    % save results as Matlab file
    save([resfolder '\' expName '\EvalData_' expName '.mat'],'Eval','partProp','C','pxsz')

    if saveas == 'xlsx'
        % Write to Excel file
        writetable(partProp,[resfolder '\' expName '\EvalData_' expName '.xlsx'])
        xlswrite([resfolder '\' expName '\EvalData_' expName '.xlsx'],Eval,'MeanValues','B2')
        xlswrite([resfolder '\' expName '\EvalData_' expName '.xlsx'],{'mean Value','standard deviation'},'MeanValues','B1')
        rows={'Area [nm^2]'; 'Equivalent Diameter [nm]';...
            'min. Feret Diameter [nm]'; 'max. Feret Diameter [nm]';'Circularity [-]'};
        xlswrite([resfolder '\' expName '\EvalData_' expName '.xlsx'],rows,'MeanValues','A2')
        tend=toc
    elseif saveas =='csv'
        writetable(partProp,[resfolder '\' expName '\EvalData_' expName '.csv'])
    end

if batch == 1  && allinone > 0
    AllpartProp = [AllpartProp; partProp];
end
end

if batch == 1 && allinone > 0
    writetable(AllpartProp,[resfolder '\' expName '_all_Images_Results.xlsx'])
end
