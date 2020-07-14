function [innercells,outercells] = importtracing(innerfilename,outerfilename)
% An optional function to import any existing trace data into the workspace. 
% files must be in the current folder directory and in a .mat format
% 
% IN
% innerfilename: file containing inner circle data, e.g. 'inners.mat'
% outerfilename: file containing outer circle data, e.g. 'outers.mat'
%
% OUT
% innercells: a cell array containing the imported inner circle data
% outercells: a cell array containing the imported outer circle data
%
% Sample usage:
% [innerArray,outerArray]=importtracing('innersfile.mat','outersfile.mat');
%
%
% Nishant '23 | singhal@princeton.edu | July 2020


% if specified files exist, import them
if isfile(innerfilename)
    innerstemp=load(innerfilename);
    inners=struct2cell(innerstemp);
    inners=inners{1};
end
if isfile(outerfilename)
    outerstemp=load(outerfilename);
    outers=struct2cell(outerstemp);
    outers=outers{1};
end

innercells=inners;
outercells=outers;


end

