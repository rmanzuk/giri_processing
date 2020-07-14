function [] = plottracing(innersInput,outersInput)
% Plots archeo trace data. This function has no actual outputs. 
% 
% IN
% innersInput: the cell array containing inner circle data
% innersOutput: the cell array containing outer circle data
% 
% Inner data is plotted as black rings. 
% Outer data is plotted as red rings. 
% Outlier data is plotted as red/black X's. 
% 
% Sample usage:
% plottracing(innerArray,outerArray)
% 
% 
% Nishant '23 | singhal@princeton.edu | July 2020



inners=innersInput;
outers=outersInput;

%construct the four sub-arrays:
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

for i=1:length(inners)
    if sum(size(inners{i}))>0
        innersguesses{i}=inners{i}(inners{i}(:,3)==120,1:2);
    end
end

for i=1:length(outers)
    if sum(size(outers{i}))>0
        outersguesses{i}=outers{i}(outers{i}(:,3)==120,1:2);
    end
end
clear i

layerseparation=-3.125;
hold on

% plot inner circle ring
for i=1:(sum(size(innersconfident))-1)
    [numpoints,~] = size(innersconfident{i});
    if numpoints>0
        zvalues=layerseparation*(i-1)*ones(numpoints+1,1);
        plot3([innersconfident{i}(:,1); innersconfident{i}(1,1)], [innersconfident{i}(:,2); innersconfident{i}(1,2)],zvalues,'k-')
    end
    clear numpoints zvalues
end

% plot outer circle ring
for i=1:(sum(size(outersconfident))-1)
    [numpoints,~] = size(outersconfident{i});
    if numpoints>0
        zvalues=layerseparation*(i-1)*ones(numpoints+1,1);
        plot3([outersconfident{i}(:,1); outersconfident{i}(1,1)], [outersconfident{i}(:,2); outersconfident{i}(1,2)],zvalues,'r-')
    end
    clear numpoints zvalues
end

% plot inner circle outliers
for i=(sum(size(innersguesses))-1)
    [numpoints,~] = size(innersguesses{i});
    if numpoints>0
        zvalues=layerseparation*(i-1)*ones(numpoints+1,1);
        plot3([innersguesses{i}(:,1); innersguesses{i}(1,1)], [innersguesses{i}(:,2); innersguesses{i}(1,2)],zvalues,'k.')
    end
    clear numpoints zvalues
end

% plot outer circle outliers
for i=(sum(size(outersguesses))-1)
   [numpoints,~] = size(outersguesses{i});
    if numpoints>0
        zvalues=layerseparation*(i-1)*ones(numpoints+1,1);
        plot3([outersguesses{i}(:,1); outersguesses{i}(1,1)], [outersguesses{i}(:,2); outersguesses{i}(1,2)],zvalues,'r.')
    end 
    clear numpoints zvalues
end

clear i
axis([0 748 0 538])
box on
axis ij
daspect([1 1 1])
end

