function [myAngle] = branchangle(vector1,vector2)
% Receives two direction vectors and figures out the angle between them,
% after first "normalising" by making both vectors point upwards. 
% 
% IN
% vector1: the first direction vector. A 3x1 matrix
% vector2: the second direction vector. A 3x1 matrix
% 
% OUT
% myAngle: the angle in degress
% 
% Sample usage:
% archeoAngle=branchangle(V1,V2)
% 
% 
% Nishant '23 | singhal@princeton.edu | July 2022


tempV1=vector1;
tempV2=vector2;


% Make both vectors point in the positive z direction
% i.e. reverse their direction if they point downards
if tempV1(3)<0
    tempV1=-1*tempV1;
end
if tempV2(3)<0
    tempV2=-1*tempV2;
end

mycos = dot(tempV1,tempV2)/(norm(tempV1)*norm(tempV2));
myAngle = acosd(mycos);
end

