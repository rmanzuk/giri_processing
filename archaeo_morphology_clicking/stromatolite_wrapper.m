% script to run through and make a few stromatolites, coalesc them and
% measure their params
%
% R. A. Manzuk 03/29/2021

%% load in database stuff
ages = xlsread('PreC_out_20191012.xlsx','Sheet1','O3:O16321');

bed_base_top = xlsread('PreC_out_20191012.xlsx','Sheet1','R3:S16321');
thicknesses = bed_base_top(:,2) - bed_base_top(:,1);

[~,strom_notes] =  xlsread('PreC_out_20191012.xlsx','Sheet1','Y3:Y16321');
%% do some analysis

% we only really need to worry about ~700ma to 500ma where there is a strom
% note
precam_cam = ages < 700 & ages > 500 & ~strcmp(strom_notes,'none');

final_thicknesses = thicknesses(precam_cam);
final_ages = ages(precam_cam);
final_stroms = strom_notes(precam_cam);

% need to split up the stromatolites when multiple are listed
split_stroms = cellfun(@(x) strsplit(x, ', '), final_stroms, 'UniformOutput', false);

% id the classes of stroms that could build a reef
stroms_reef = {'thromb';'lds';'colbras';'constrom';'gistm'};

% and identify the occurrences of each kind
strom_occs = false(numel(final_ages),numel(stroms_reef));
for i = 1:numel(stroms_reef)
    for j = 1:numel(split_stroms)
       strom_occs(j,i) = ismember(stroms_reef(i),split_stroms{j}(:));
    end
end

% what are the unique ages
unique_ages = unique(final_ages);

% and what proportion of that thickness is made by each type of strom
strom_thicknesses = [];
for i = 1:numel(unique_ages)
    for j = 1:numel(stroms_reef)
        strom_thicknesses(i,j) = sum(final_thicknesses(final_ages == unique_ages(i) & strom_occs(:,j))); 
    end
end

proportional_thicknesses = strom_thicknesses./sum(strom_thicknesses,2);

%% now, calculate sa and enclosing volume for each type of strom
% start with thrombolites and lds (as same)

% these params hold for all stromatolites
height = 10;
vert_sampling_rate = 10;
circle_sampling_rate = 10;
bump_size = 1;
stromat_scale_ratio = 1/vert_sampling_rate;

% set of radii
rads_lds = [2,4,6,8,9,10,10,10,9,8];

stromat1 = make_stromatolite(rads_lds,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat2 = make_stromatolite(rads_lds,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat3 = make_stromatolite(rads_lds,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat4 = make_stromatolite(rads_lds,height,vert_sampling_rate,circle_sampling_rate,bump_size);

% now translate the stromatolites so they can be in different positions and we can coalesc them
translation_length = 12;

% move stromat2 in positive x direction
for i = 1:numel(stromat2)
    points_here = stromat2{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    stromat2{i} = points_here;
end

% move stromat3 in positive y direction
for i = 1:numel(stromat3)
    points_here = stromat3{i};
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat3{i} = points_here;
end

% move stromat4 in positive x and y direction
for i = 1:numel(stromat4)
    points_here = stromat4{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat4{i} = points_here;
end

% coalesc
all_stromats_lds = {stromat1,stromat2,stromat3,stromat4};

% identify branch points
[lds_branched_flags,lds_branch_points_3d] = process_branched_network(all_stromats_lds,stromat_scale_ratio);

% make those stromats 3d
lds_3d = make_clicking_3d(all_stromats_lds,stromat_scale_ratio);

[lds_surface_area,lds_volume] = sa_and_vol(all_stromats_lds,stromat_scale_ratio,1,lds_branched_flags);
[lds_convhull_points, lds_enclosing_volume] = get_enclosing_volume(lds_3d, 1);
%% columnar, branhcing

% set of radii
rads_colbras = [1,2,2.5,3,3.5,4,4,4,4,4]/2;

stromat1 = make_stromatolite(rads_colbras,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat2 = make_stromatolite(rads_colbras,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat3 = make_stromatolite(rads_colbras,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat4 = make_stromatolite(rads_colbras,height,vert_sampling_rate,circle_sampling_rate,bump_size);

% now translate the stromatolites so they can be in different positions and we can coalesc them
translation_length = 7;

% move stromat2 in positive x direction
for i = 1:numel(stromat2)
    points_here = stromat2{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    stromat2{i} = points_here;
end

% move stromat3 in positive y direction
for i = 1:numel(stromat3)
    points_here = stromat3{i};
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat3{i} = points_here;
end

% move stromat4 in positive x and y direction
for i = 1:numel(stromat4)
    points_here = stromat4{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat4{i} = points_here;
end

% coalesc
all_stromats_colbras = {stromat1,stromat2,stromat3,stromat4};

% identify branch points
[colbras_branched_flags,colbras_branch_points_3d] = process_branched_network(all_stromats_colbras,stromat_scale_ratio);

% make those stromats 3d
colbras_3d = make_clicking_3d(all_stromats_colbras,stromat_scale_ratio);

[colbras_surface_area,colbras_volume] = sa_and_vol(all_stromats_colbras,stromat_scale_ratio,1,colbras_branched_flags);
[colbras_convhull_points, colbras_enclosing_volume] = get_enclosing_volume(colbras_3d, 1);
%% conical

height = 7;
% set of radii
rads_constrom = [0:0.5:4];

stromat1 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat2 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat3 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat4 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);

% now translate the stromatolites so they can be in different positions and we can coalesc them
translation_length = 4;

% move stromat2 in positive x direction
for i = 1:numel(stromat2)
    points_here = stromat2{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    stromat2{i} = points_here;
end

% move stromat3 in positive y direction
for i = 1:numel(stromat3)
    points_here = stromat3{i};
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat3{i} = points_here;
end

% move stromat4 in positive x and y direction
for i = 1:numel(stromat4)
    points_here = stromat4{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat4{i} = points_here;
end

% coalesc
all_stromats_constrom = {stromat1,stromat2,stromat3,stromat4};

% identify branch points
[constrom_branched_flags,constrom_branch_points_3d] = process_branched_network(all_stromats_constrom,stromat_scale_ratio);

% make those stromats 3d
constrom_3d = make_clicking_3d(all_stromats_constrom,stromat_scale_ratio);

[constrom_surface_area,constrom_volume] = sa_and_vol(all_stromats_constrom,stromat_scale_ratio,1,constrom_branched_flags);
[constrom_convhull_points, constrom_enclosing_volume] = get_enclosing_volume(constrom_3d, 1);

%% giant

height = 50;
% set of radii
rads_gistm = [2,4,6,8,9,10,10,10,9,8]*4;

stromat1 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat2 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat3 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);
stromat4 = make_stromatolite(rads_constrom,height,vert_sampling_rate,circle_sampling_rate,bump_size);

% now translate the stromatolites so they can be in different positions and we can coalesc them
translation_length = 55;

% move stromat2 in positive x direction
for i = 1:numel(stromat2)
    points_here = stromat2{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    stromat2{i} = points_here;
end

% move stromat3 in positive y direction
for i = 1:numel(stromat3)
    points_here = stromat3{i};
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat3{i} = points_here;
end

% move stromat4 in positive x and y direction
for i = 1:numel(stromat4)
    points_here = stromat4{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat4{i} = points_here;
end

% coalesc
all_stromats_gistm = {stromat1,stromat2,stromat3,stromat4};

% identify branch points
[gistm_branched_flags,gistm_branch_points_3d] = process_branched_network(all_stromats_gistm,stromat_scale_ratio);

% make those stromats 3d
gistm_3d = make_clicking_3d(all_stromats_gistm,stromat_scale_ratio);

[gistm_surface_area,gistm_volume] = sa_and_vol(all_stromats_gistm,stromat_scale_ratio,1,gistm_branched_flags);
[gistm_convhull_points, gistm_enclosing_volume] = get_enclosing_volume(gistm_3d, 1);
%% finally just sa to enclosing volume over time given proportion of occurrences
sa_vol_mat = [sum(lds_surface_area)/lds_enclosing_volume;...
    sum(lds_surface_area)/lds_enclosing_volume;...
    sum(colbras_surface_area)/colbras_enclosing_volume;
    sum(constrom_surface_area)/constrom_enclosing_volume;...
    sum(gistm_surface_area)/gistm_enclosing_volume];

sa_vol_bins = sa_vol_mat' * proportional_thicknesses';