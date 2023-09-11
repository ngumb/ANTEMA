function SepParticles = ParticleSeparationWatershed(MnonConv,Convexthresh,AspCutoff)
% Ultimate Erosion of Convex Shapes (UECS) based particle separation of non convex particles
%   based on the discribtion of 
%   Park, C. and Ding, Y. (2021) Data Science for Nano Image Analysis, Springer Nature, Cham, Switzerland
%   and the  corresponding code available at https://aml.engr.tamu.edu/book-dsnia/
% INPUT
%   MnonConv:       binary map of a non convex particle shape

% OUTPUT
%   SepParticles:   cell array of separated particle maps
%   markers:        binary map of all generated markers

bw = MnonConv;
D = -bwdist(~bw);
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;


BB=regionprops(bw3,'BoundingBox');
BB = struct2cell(BB);
stats = regionprops(bw3,'Perimeter','ConvexImage','MaxFeretProperties','MinFeretProperties');

for ii=1:length(stats)
   cp = regionprops(+stats(ii).ConvexImage,'Perimeter');
   stats(ii).Convexity = cp(1).Perimeter / stats(ii).Perimeter;
   
end
AspR = [stats.('MaxFeretDiameter')]./[stats.('MinFeretDiameter')];
Convex = [stats.Convexity]';
clear stats
Mstack = regionprops('table', bw3, 'Image');
Mstack = Mstack.Image;

for ii = 1:length(Mstack)
    BBcurr = BB{ii};
    BBcurr = ceil(BBcurr);
    if Convex(ii) <= Convexthresh 
        
        MnonConv = Mstack{ii};
        Msplit = bw3(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1);
        bw3(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
            ~(Msplit==MnonConv);
    elseif AspR(ii) < AspCutoff
        MnonConv = Mstack{ii};
        Msplit = bw3(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1);
        bw3(BBcurr(2):BBcurr(2)+BBcurr(4)-1,BBcurr(1):BBcurr(1)+BBcurr(3)-1,1) =...
            ~(Msplit==MnonConv);
    end
end

SepParticles = bw3;
%% TEST CONVEXITY AND REMOVE NONCONVEX FROM IMAGE
% SepParticles= {};
% [B, L] = bwboundaries(bw3,8,'noholes');
% num_markers = size(B,1);
% 
% for j = 1:num_markers
%     pls = logical(L==j);
% 
%     s  = regionprops(pls,'ConvexImage','Perimeter','Area','ConvexArea');
%     cp = regionprops(+s(1).ConvexImage,'Perimeter');
%     Conv = cp(1).Perimeter / s(1).Perimeter;
% 
%     if Conv <= Convthresh 
%         continue
%     end
% 
%     SepParticles{j} = pls;
% end

% 
end