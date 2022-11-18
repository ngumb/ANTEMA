function [SepParticles,markers] = ParticleSeparation(MnonConv,minMarkerSize,Convthresh)
% Ultimate Erosion of Convex Shapes (UECS) based particle separation of non convex particles
%   based on the discribtion of 
%   Park, C. and Ding, Y. (2021) Data Science for Nano Image Analysis, Springer Nature, Cham, Switzerland
%   and the  corresponding code available at https://aml.engr.tamu.edu/book-dsnia/
% INPUT
%   MnonConv:       binary map of a non convex particle shape
%   minMarkerSize:  minimum pixel area size a marker has to be to still be considered
%   Convexthresh:   threshold for particles to be considered convex
% OUTPUT
%   SepParticles:   cell array of separated particle maps
%   markers:        binary map of all generated markers


Temp_Ero = MnonConv;
markers = 0.*Temp_Ero;
Temp_Count = 0.*Temp_Ero;

se1    = strel('disk', 1 );
se2    = strel('arbitrary', ones(2,2));


minMarkerSize = max([minMarkerSize,30]);
updated = 1; cnt = 0;

while updated
    updated = 0;
    [B, L] = bwboundaries(Temp_Ero,8,'noholes');
    s  = regionprops(L, 'ConvexArea', 'Area', 'Perimeter','ConvexImage','MinorAxisLength');
    if mod(cnt, 2) == 1
        se = se1;
    else
        se = se2;
    end
    
    for j = 1:size(B,1)

        pls = logical(L==j);

        cp = regionprops(+s(j).ConvexImage,'Perimeter');
        Conv = cp(1).Perimeter / s(j).Perimeter;

% Maxmimum iterations can also be implemented if needed
        if ((Conv <= Convthresh || 1 - s(j).Area / s(j).ConvexArea >= 0.05) && s(j).Area > minMarkerSize)&& s(j).MinorAxisLength > floor(minMarkerSize / 10.0)
            tmp = imerode(pls, se);
            tmp = imopen(tmp, se1);
            updated = 1;
            Temp_Ero(pls) = tmp(pls);
        else
            if s(j).Area > minMarkerSize
                markers(pls) = 1;
            end
            Temp_Ero(pls) = 0;
            Temp_Count(pls) = cnt +1;
        end
    end
    cnt = cnt + 1;
end


%% markers are generated, regrow

SepParticles= {};
if nnz(markers) == 0
    return
end

[B, L] = bwboundaries(markers,8,'noholes');
num_markers = size(B,1);


for j = 1:num_markers
    pls = logical(L==j);
    cnt_ero = max(unique(Temp_Count(pls)));

    updated = 1; cnt = 0;
    while updated

        if cnt + 1 > cnt_ero
            updated = 0;
            continue;
        end

        if mod(cnt, 2) == 1
            se = se1;
        else
            se = se2;
             updated = 0;
        end
      
        tmp = imdilate(pls, se);
        tmp = imclose(tmp, se1);
        
        tmp2 = tmp & (MnonConv > 0) & ~pls;
        tmp3 = tmp & ~pls;

        if (sum(sum(tmp2)) > 0) && (nnz(tmp2)/nnz(tmp3) >= 0.25) 
            updated = 1;
        end

        pls = tmp2 | pls;
        
        
        cnt = cnt + 1;
    end
    s  = regionprops(pls,'ConvexImage','Perimeter','Area','ConvexArea');
    cp = regionprops(+s(1).ConvexImage,'Perimeter');
    Conv = cp(1).Perimeter / s(1).Perimeter;

    if Conv <= Convthresh || 1 - s.Area / s.ConvexArea >= 0.05
        continue
    end

    SepParticles{j} = pls;
    
    
end


end