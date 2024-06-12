function [Eval,T_prop,RGB,Mlines] = particlePropertiesEval(map,Class,pxsz,Image, varargin)
% Function to give out particle properties for a HRTEM image segmentation 
% 
% INPUTS:
%   map:      Segmentation map, uint8 file
%   Class:    integer
%             Class to analyze
% Recognized Name-Value Pair arguments
%   pxsz:     pixel size in nm
%   Convexthresh:  threshold for particle separation, 
%                  default value 0.95  
%   PSep:     Logical Value, if true a separation algorithm will be used on
%               non convex particles, default true
%   minmarker: minimum marker size for convex markers to regrow in
%              separation routine, default 0.5
%   removecount: remove particles smaller than input in nm, default 0.5
%   Convexthresh_min: minimum Convexity needed to feed into separation
%                     routine, particles will be discarded if below,
%                     default 0.5
%   Mode : Mode for particle separation. Can either be 'UECS' or 'Watershed'
%           Default: UECS
% OUTPUTS:
%   Eval: Eval=[p_Area_mean, p_Area_StD;...
%        p_LFmin_mean, p_LFmin_StD;...
%        p_LFmax_mean, p_LFmax_StD;...
%        p_circ_mean, p_circ_StD];
%   T_prop: Particle Properties as a Table
%   RGB:    Processed Particle Map
%   Mlines: Outlines of the separated particles to visualize overlapping
%   regions
%
% Copyright (C) 2022  University of Duisburg-Essen
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

% Default Values
Convexthresh = 0.95;
PSep = true;
minmarker = 0.5;
removecount = 0.5;
Convexthresh_min = 0.5;
defaultMode = 'UECS';
expectedModes = {'UECS','Watershed'};


p = inputParser;
addParameter(p,'Convexthresh',Convexthresh,@(x) isnumeric(x)&& (x<=1)&&(x>=0));
addParameter(p,'PSep', PSep, @islogical);
addParameter(p,'minmarker',minmarker, @isnumeric);
addParameter(p,'removecount',removecount, @isnumeric);
addParameter(p,'Convexthresh_min',Convexthresh_min,@(x) isnumeric(x)&& (x<=1)&&(x>=0));
addParameter(p,'Mode',defaultMode,@(x) any(validatestring(x,expectedModes)));
parse(p,varargin{:})

Convexthresh = p.Results.Convexthresh;
PSep = p.Results.PSep;
minmarker = p.Results.minmarker;
removecount = p.Results.removecount;
Convexthresh_min = p.Results.Convexthresh_min;
Mode = p.Results.Mode;

    %get binary map for single Particle Class
    M1=ismember(map,Class);

    % remove all areas with an equivalence diameter of 0.5 nm
    removecount = removecount/pxsz;
    removecount = double(ceil(pi/4*(removecount)^2));
    minmarker = minmarker/pxsz;
    minmarker = double(ceil(pi/4*(minmarker)^2));
    M1 = bwareaopen(M1,removecount);

    BB=regionprops(M1,'BoundingBox');
    BB = struct2cell(BB);

    % get convexity 
    %% Particle Separation
    
    Mresolved = M1;
    Mlines = zeros(size(M1,1),size(M1,2));

    if PSep == true 
        stats = regionprops(M1,'Perimeter','ConvexImage','EulerNumber');
            for ii=1:length(stats)
               cp = regionprops(+stats(ii).ConvexImage,'Perimeter');
               stats(ii).Convexity = cp(1).Perimeter / stats(ii).Perimeter;
               
            end
        
        Convex = [stats.Convexity]';
        
        Mstack = regionprops('table', M1, 'Image');
        Mstack = Mstack.Image;

        for ii = 1:length(Mstack)
            BBcurr = BB{ii};
            BBcurr = ceil(BBcurr);
            if Convex(ii) <= Convexthresh && Convex(ii) >= Convexthresh_min
                
                MnonConv = Mstack{ii};
                Msplit = Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1);
                Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
                    ~(Msplit==MnonConv);
                if stats(ii).EulerNumber <= -10
                    warning('Too many holes. Particle will be removed from evaluation \n')
                    continue
                end
                if strcmp(Mode,'UECS') == true
                    [SepParticles,~] = ParticleSeparation(MnonConv,minmarker,Convexthresh);
    
                    for bi = 1:length(SepParticles)
                        if ~isempty(SepParticles{bi})
                            Msplit = zeros(size(M1,1),size(M1,2));
                            Msplit(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1) = SepParticles{bi};
                            Mresolved(:,:,end+1) = Msplit;
                            B = bwperim(Msplit,8);
                            Mlines(B) = 1; 
                        end
                    end
                elseif strcmp(Mode,'Watershed') == true
                    SepParticles = ParticleSeparationWatershed(MnonConv,Convexthresh);
                    Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
                    Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1)+SepParticles;
                    % OVERRIDE PARTICLES WITH SEPARATED TO SAVE MEMORY
                else
                    error('Input for mode is not a valid input. Mode can be UECS oder Watershed \n')
                end
            elseif Convex(ii) < Convexthresh_min %QUALITY CONTROL
                MnonConv = Mstack{ii};
                Msplit = Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1);
                Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
                    ~(Msplit==MnonConv);
                warning('Convexity too low. Particle will be excluded from evaluation.\n')
                               
            end
        end
    end
    
    partProp = [];
    SepCount = 0;
    for ii = 1:size(Mresolved,3)
            
             M = Mresolved(:,:,ii);  


             M = imclearborder(M);
             M = bwareaopen(M,removecount);
             M = imfill(M,'holes');
             Mresolved(:,:,ii) = M;
             
             if M == 0
                 continue
             end

             % Get regionproperties of all particles as table
             partPropS=regionprops('table',M,'Area','Circularity','EquivDiameter',...
                'MaxFeretProperties','MinFeretProperties','Perimeter','BoundingBox','Centroid');
             %remove unneccessary parameters
             partPropS=removevars(partPropS,{'MaxFeretAngle','MaxFeretCoordinates',...
                'MinFeretAngle','MinFeretCoordinates'});

             if ii > 1
                SepCount= SepCount + height(partPropS);
             end

             if PSep == false
                 stats = regionprops(M,'Perimeter','ConvexImage');
                    for ii=1:length(stats)
                       cp = regionprops(+stats(ii).ConvexImage,'Perimeter');
                       stats(ii).Convexity = cp(1).Perimeter / stats(ii).Perimeter;
                       
                    end
                
                Convex = [stats.Convexity]';
             end
             partProp =[partProp; partPropS];
            
    end
    Mresvis = logical(Mresolved);
    RGB = insertObjectMask(Image,Mresvis,'LineColor','white','LineWidth',2);
    
    
    %% Area
    P_Area=partProp.('Area');
    p_Area=P_Area.*(pxsz^2);
    p_Area_mean=mean(p_Area);
    p_Area_StD=std(p_Area);
    
    %% Equivalence Diameter
    P_Dequ=partProp.('EquivDiameter');
    p_Dequ=P_Dequ.*pxsz;
    p_Dequ_mean=mean(p_Dequ);
    p_Dequ_StD=std(p_Dequ);
    
    %% Ferret min
    P_LFmin=partProp.('MinFeretDiameter');
    p_LFmin=P_LFmin.*pxsz;
    p_LFmin_mean=mean(p_LFmin);
    p_LFmin_StD=std(p_LFmin);
    
     %% Ferret max
    P_LFmax=partProp.('MaxFeretDiameter');
    p_LFmax=P_LFmax.*pxsz;
    p_LFmax_mean=mean(p_LFmax);
    p_LFmax_StD=std(p_LFmax);

     %% Aspect Ratio
    AspR = p_LFmax./p_LFmin;

    %% Circularity
    p_circ=partProp.('Circularity');
    p_circ_mean=mean(p_circ);
    p_circ_StD=std(p_circ);

    %% Perimeter
    Perim = partProp.('Perimeter');
    Perim = Perim.*pxsz;
    %% Bounding Box
    BB = partProp.('BoundingBox');

    %% Centroid
    Centroid2 = partProp.('Centroid');
    Centroid = Centroid2.*pxsz;
    
    %% Separation indicator
    if PSep == true
        Sep = zeros(height(partProp),1);
        Sep(end-SepCount:end)=1;
    else
        NonConv = Convex <= Convexthresh;
    end

    %% Save all mean values in one matrix
    Eval=[p_Area_mean, p_Area_StD;...
        p_Dequ_mean, p_Dequ_StD;
        p_LFmin_mean, p_LFmin_StD;...
        p_LFmax_mean, p_LFmax_StD;...
        p_circ_mean, p_circ_StD];
    
    %Create Table
    if PSep == true
        T_prop=table(p_Area, p_Dequ, p_LFmin, p_LFmax, AspR, p_circ,Perim,BB,Sep,Centroid);
        T_prop.Properties.VariableNames = {'Area [nm^2]', 'Equivalent Diameter [nm]',...
            'min. Feret Diameter [nm]', 'max. Feret Diameter [nm]','Aspect Ratio (max/min)','Circularity [-]',...
            'Perimeter [nm]','Bounding Box parameters (xmin,ymin,xwidth,ywidth) [pixel]','Separation used (1 = true)','Centroid position [nm]'};
    else
        T_prop=table(p_Area, p_Dequ, p_LFmin, p_LFmax, AspR, p_circ,Perim,BB,Convex,NonConv,Centroid);
        T_prop.Properties.VariableNames = {'Area [nm^2]', 'Equivalent Diameter [nm]',...
            'min. Feret Diameter [nm]', 'max. Feret Diameter [nm]','Aspect Ratio (max/min)','Circularity [-]',...
            'Perimeter [nm]','Bounding Box parameters (xmin,ymin,xwidth,ywidth) [pixel]','Convexity [-]','Convexity below threshold (1 = true)','Centroid position [nm]'};
    end
   
   