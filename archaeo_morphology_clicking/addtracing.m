function [innercells,outercells] = addtracing(imageDirectory,startNum,endNum,existingInners,existingOuters)
% This function lets you create archeo trace data using ginput,
% for a single archeo. 
% It will loop through a set of photos and ask you to input inner and
% outer circle points. 
% 
% Begin by changing the code below to point matlab to the folder containing
% your images in .tif format. 
% 
% IN
% imageDirectory: fullfile directory with the folder where the images are
% stored and their extension. example; imageDirectory = 
% dir(fullfile('/Users/Nishant/Desktop/PEI/radius/', '*.tif'));
% startNum: what image # you want to start at. Eg if you only wanna
% input tracing data for a subset of the images in your stack
% endNum: what image # in your stack you want to end at. 
% existingInners: (optional) specify the cell array that already contains 
% inner circle data for the archeo, so that you can add to it
% rather than starting fresh. This allows you to trace an archeo across
% various sittings rather than all at once. 
% existingOuters: (optional) specify the cell array containing any existing
% outer circle data for the archeo you're currently working on
% 
% OUT
% innercells: create a cell array containing the traced inner circle data
% outercells: create a cell array containing the traced outer circle data
% 
% Matlab will ask for inner circle data, then outer circle data, then
% cycle to the next image, etc. Use the enter key to move through
% sections. Use the left-click to input points (in circular order), and 
% use the 'x' key to input outlier or "guess" points. 
% 
% There is some commented code at the bottom that you can use to export
% traced data to .mat files. Use this to save your work between sittings. 
% 
% Sample usage:
% [innernew,outernew]=addtracing(3,8,innerold,outerold);
%   This will allow you to input trace data for images 3 through 8 in
%   your folder. innerold and outerold  would be cell arrays containing
%   existing inner and outer data for the archeo (e.g. perhaps for images
%   1?5). 
%   The function will create new cell arrays, innernew and outernew, 
%   conrtaining the updated trace data. Each cell in the array represents
%   a layer. For example, innernew{4} would give me the inner points for 
%   the fourth image. Rows with a '1' in the 3rd column are points in the
%   circle, and rows with a '120' in the 3rd  column are outlier points. 
% 
%
% Nishant '23 | singhal@princeton.edu | July 2020


% CHANGE THIS to point to the folder containing your images
%filenames = dir(fullfile('/Users/Nishant/Desktop/PEI/radius/', '*.tif'));
fulln=length(imageDirectory);


% throws error if number of inputs isn't 2 or 4
if nargin ~= 3 && nargin ~= 5
  error('number of inputs must be 2 or 4')
end

% makes sure that the tracings are added to existing data rather than
% entirely replacing it
if nargin == 5
    inners=existingInners;
    outers=existingOuters;
    if startNum > numel(inners)
        add_empty = cell(1,endNum - numel(inners));
        inners(numel(inners)+1:endNum) = add_empty;
        outers(numel(outers)+1:endNum) = add_empty;
    else
        %do nothing
    end
else
    inners=cell(1,endNum);
    outers=cell(1,endNum);
end

% error messages to ensure that startNum ? endNum ? total no. of images
if startNum > endNum
    error('startNum cannot be bigger than endNum')
end
if endNum > fulln
    error('endNum cannot be greater than your total number of images')
end

% goes through each image in range, and uses ginput to get inners + outers
count = 0;
for i=startNum:endNum
    count = count+1;
    
    workingimage=imread(fullfile(imageDirectory(i).folder,imageDirectory(i).name));
    
    imshow(workingimage,'InitialMagnification',400);
    ax = gca;
    ax.Toolbar.Visible = 'off';
    hold on
    if i>1  && ~isempty(inners{i-1})
        plot(inners{i-1}(:,1),inners{i-1}(:,2),':','LineWidth',1)
    else
        % do nothing
    end
    title("inner for image #" + (i) + " (" + (i-startNum+1) + " of " + (endNum-startNum+1) + " for this tracing)"); 
    n=0;
    while true
        [x_i,y_i,button_i] = ginput(1);
        if isempty(x_i) ; break; end
        n = n+1;
        x(n) = x_i(1);
        y(n) = y_i(1);
        button(n) = button_i(1);
        plot(x,y,'r','LineWidth',1)
        drawnow
    end
    
    if exist('x')
        inners{i}=[x' y' button'];
    end
    clear x y button
    hold off
    
    imshow(workingimage,'InitialMagnification',400); 
    ax = gca;
    ax.Toolbar.Visible = 'off';
    hold on
    if i>1 && ~isempty(outers{i-1})
        plot(outers{i-1}(:,1),outers{i-1}(:,2),':','LineWidth',1)
    else
        % do nothing
    end
    title("outer for image #" + (i) + " (" + (i-startNum+1) + " of " + (endNum-startNum+1) + " for this tracing)"); 
    n=0;
    while true
        [x_i,y_i,button_i] = ginput(1);
        if isempty(x_i) ; break; end
        n = n+1;
        x(n) = x_i(1);
        y(n) = y_i(1);
        button(n) = button_i(1);
        plot(x,y,'r','LineWidth',1)
        drawnow
    end
    if exist('x')
        outers{i}=[x' y' button'];
    end
    clear x y button

end

innercells=inners;
outercells=outers;

% optionally use the following code to save the results to files:
% save('innersfile.mat', 'innercells')
% save('outersfile.mat', 'outercells')
end

