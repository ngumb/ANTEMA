function [Image] = scalebar(I,pxsz,unit)
% writes a scalebar into the image
% To be used with full size HRTEM images of Size 1024x1024 - 2048x2048
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


% check if unit is in nm
if unit ~= 'nm'
    error('unit of pixel scale not in nm, Case currently not implemented\n')
end

% Transfer image to uint8 if not already uint8
I = im2uint8(I);

% Basic bar parameters
barThick=13; %Tickness of bar center
barFrame=15; %Thickness of additional bar frame
%colorCenter=[255,255,255]; %white
%colorFrame=[0,0,0]; %black

%%  Decide bar length
% select upper and lower boundary
upperN = pxsz*350;
lowerN = pxsz*150;
% round to integer and build range
upperN = round(upperN);
lowerN = round(lowerN);
range=[lowerN:1:upperN];
% check if multiple of 5 is within range
multi5=(mod(range,5)==0);
range5=range(multi5);

if isempty(range5)
   n=lowerN+floor((upperN-lowerN)/2);
else
    n=min(range5);
end
barlength = floor(n/pxsz);

% Position of scalebar in lower right corner
% d px from the right hand side and the bottom
d = 0+barFrame;

% Font Size
FS=40;

% image Size
imSzy = size(I,1);
imSzx = size(I,2);

% outer boundary box
ymin = imSzy-d-barlength-2*barFrame;
ymax = imSzy-d+barFrame;
xmin = imSzx-d-barThick-2*barFrame-ceil(1.333*FS);
xmax = imSzx-d+barFrame;

% outer boundary box black
I(xmin:xmax,ymin:ymax)=0;

% inner boundary
ymin = imSzy-d-barlength-barFrame;
ymax = imSzy-d;
xmin = imSzx-d-barThick;
xmax = imSzx-d;

% inner boundary white
I(xmin:xmax,ymin:ymax)=255;

% insert scaling text
xpos = xmin+floor((xmax-xmin)/2);
ypos = ymin+floor((ymax-ymin)/2);
text_str = [num2str(n) ' nm']; %case of other units not implemented
Image = insertText(I,[ypos,xpos],text_str,'FontSize',FS,...
    'TextColor','white','BoxOpacity',0,'AnchorPoint','CenterBottom',...
    'Font','Arial Black');

end

