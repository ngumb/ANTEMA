function [Eval,T_prop,RGB,Mlines] = particlePropertiesEval(map,Class,pxsz,Convexthresh,PSep,Image)
% Function to give out particle properties for a HRTEM image segmentation 
% 
% INPUTS:
%   map:      Segmentation map, uint8 file
%   Class:    integer
%             Class to analyze
%   pxsz:     pixel size in nm
%   Convexthresh:  threshold for particle separation, 
%                  particles with a convexity below 
%   PSep:     Logical Value, if true a separation algorithm will be used on
%               non convex particles
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


    %get binary map for single Particle Class
    M1=ismember(map,Class);

    % remove all areas with an equivalence diameter of 0.5 nm
    removecount = 0.5/pxsz;
    removecount = double(ceil(pi/4*(removecount)^2));
    minmarker = 0.5/pxsz;
    minmarker = double(ceil(pi/4*(minmarker)^2));
    M1 = bwareaopen(M1,removecount);

    BB=regionprops(M1,'BoundingBox');
    BB = struct2cell(BB);
  
    
    % get convexity 
    %% Particle Separation
    
    Mresolved = M1;
    Mlines = zeros(size(M1,1),size(M1,2));

    if PSep == true 
        stats = regionprops(M1,'Perimeter','ConvexImage');
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
            if Convex(ii) <= Convexthresh
                MnonConv = Mstack{ii};
                Msplit = Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1);
                Mresolved(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
                    ~(Msplit==MnonConv);
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
                'MaxFeretProperties','MinFeretProperties','Perimeter','BoundingBox');
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
    
    %% Circularity
    p_circ=partProp.('Circularity');
    p_circ_mean=mean(p_circ);
    p_circ_StD=std(p_circ);

    %% Perimeter
    Perim = partProp.('Perimeter');
    Perim = Perim.*pxsz;
    %% Bounding Box
    BB = partProp.('BoundingBox');
    
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
        T_prop=table(p_Area, p_Dequ, p_LFmin, p_LFmax, p_circ,Perim,BB,Sep);
        T_prop.Properties.VariableNames = {'Area [nm^2]', 'Equivalent Diameter [nm]',...
            'min. Feret Diameter [nm]', 'max. Feret Diameter [nm]','Circularity [-]',...
            'Perimeter [nm]','Bounding Box parameters (xmin,ymin,xwidth,ywidth) [pixel]','Separation used (1 = true)'};
    else
        T_prop=table(p_Area, p_Dequ, p_LFmin, p_LFmax, p_circ,Perim,BB,Convex,NonConv);
        T_prop.Properties.VariableNames = {'Area [nm^2]', 'Equivalent Diameter [nm]',...
            'min. Feret Diameter [nm]', 'max. Feret Diameter [nm]','Circularity [-]',...
            'Perimeter [nm]','Bounding Box parameters (xmin,ymin,xwidth,ywidth) [pixel]','Convexity [-]','Convexity below threshold (1 = true)'};
    end
   
   