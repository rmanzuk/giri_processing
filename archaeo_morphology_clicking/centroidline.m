function [centerpoints,directionVector] = centroidline(innersInput,outersInput)
% Figures out centroid points for each layer, plots them, then plots a line
% of best fit through them and gives the direction vector of that line.
%
% User can specify whether to use inner or outer rings to calculate
% centroid points (see code below). 
% 
% IN
% innersInput: the cell array containing inner circle data
% innersOutput: the cell array containing outer circle data
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
for i=1:length(inners)
    if sum(size(inners{i}))>0
        innersconfident{i}=inners{i}(inners{i}(:,3)==1,1:2);
    end
end
for i=1:length(outers)
    if sum(size(outers{i}))>0
        outersconfident{i}=outers{i}(outers{i}(:,3)==1,1:2);
    end
end
clear i


% specify which circle to use for centroids (inner or outer) - 
% either innersconfident or outersconfident
tempcell=outersconfident;

layerseparation=-3.125;
centerpoints=[];

% create the array of polyshapes + extracts their centroids
for i=1:length(tempcell)
    if sum(size(tempcell{i}))>0
        pgon{i}=polyshape(tempcell{i}(:,1),tempcell{i}(:,2));
        [x,y]=centroid(pgon{i});
        centerpoints=[centerpoints [x;y;layerseparation*(i-1)]];
        clear x y
    end
end

% plot centerpoints
scatter3(centerpoints(1,:),centerpoints(2,:),centerpoints(3,:),'bo')

% plot line of best fit
xyz=centerpoints';
r0=mean(xyz);
xyz=bsxfun(@minus,xyz,r0);
[~,~,V]=svd(xyz,0);
r0=r0';
V=V(:,1);
t=linspace(-1*r0(3)/V(3),(layerseparation*(length(tempcell)-1)-r0(3))/V(3));
line(r0(1)+t*V(1),r0(2)+t*V(2),r0(3)+t*V(3),'Color','blue','Linewidth',2)

directionVector=V;

end

