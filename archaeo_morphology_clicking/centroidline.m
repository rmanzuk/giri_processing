function [centerpoints,directionVector] = centroidline(innersInput,outersInput,in_or_out,scale_ratio,plot)
% Figures out centroid points for each layer, plots them, then plots a line
% of best fit through them and gives the direction vector of that line.
%
% User can specify whether to use inner or outer rings to calculate
% centroid points (see code below). 
% 
% IN
% innersInput: the cell array containing inner circle data
% innersOutput: the cell array containing outer circle data
% in_or_out: logical flag to say whether to use inner circles or outer
% circles when calculating center points and direction vecor for the line
% of best fit. 1 for inners, 0 for outers.
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
% plot: logical flag if you want the function to plot out center points and
% centroid line. 1 for plot, 0 for do not
% 
% OUT
% centerpoints: matrix containing coordinates of the center points for
% each layer
% directionVector: a 3x1 direction vector of the line of best fit. 
% 
% Sample usage:
% [centroidMatrix,V1]=centroidline(innerArray,outerArray)
% 
% 
% Nishant '23 | singhal@princeton.edu | July 2020

inners=innersInput;
outers=outersInput;

% construct the 'confident' sub-arrays:
innersconfident = {};
for i=1:length(inners)
    if sum(size(inners{i}))>0
        innersconfident{i}=inners{i}(inners{i}(:,3)==1,1:2);
    end
end
outersconfident = {};
for i=1:length(outers)
    if sum(size(outers{i}))>0
        outersconfident{i}=outers{i}(outers{i}(:,3)==1,1:2);
    end
end
clear i


% specify which circle to use for centroids (inner or outer) - 
% either innersconfident or outersconfident
if in_or_out
    tempcell=innersconfident;
    tempcell2 = outersconfident;
else
    tempcell=outersconfident;
    tempcell2 = innersconfident;
end

layerseparation=-scale_ratio;
centerpoints=[];

warning('off','all');
% create the array of polyshapes + extracts their centroids
for i=1:length(tempcell)
    if sum(size(tempcell{i}))>0
        ellipse = fit_ellipse(tempcell{i}(:,1),tempcell{i}(:,2));
        if isempty(ellipse) || strcmp(ellipse.status, 'Hyperbola found')
            pgon=polyshape(tempcell{i}(:,1),tempcell{i}(:,2));
            [x,y]=centroid(pgon);
        else
            x = ellipse.X0_in;
            y = ellipse.Y0_in;
        end
        centerpoints(:,i)=[x;y;layerseparation*(i)];
        clear x y
    end
end
warning('on','all')
centerpoints(centerpoints == 0) = NaN;


% plot line of best fit
xyz=reshape(centerpoints(~isnan(centerpoints)),3,[])';
r0=mean(xyz);
xyz=bsxfun(@minus,xyz,r0);
[~,~,V]=svd(xyz,0);
r0=r0';
V=V(:,1);

directionVector=V;

% plot centerpoints
if plot
    scatter3(centerpoints(1,:),centerpoints(2,:),centerpoints(3,:),'bo')
    t=linspace(-1*r0(3)/V(3),(layerseparation*(length(tempcell)-1)-r0(3))/V(3));
    line(r0(1)+t*V(1),r0(2)+t*V(2),r0(3)+t*V(3),'Color','blue','Linewidth',2)
else
    % do nothing
end

    
end

