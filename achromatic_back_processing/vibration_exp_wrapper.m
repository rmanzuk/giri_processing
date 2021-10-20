% script to run through images of the target from vibration experiments and
% track relative dot positions
%
% R. A. Manzuk 03/19/2021
%% set up directories, files

raw_ext = '*.iiq';
tif_ext = '*.tif';

raw_dir = '/Users/ryan/Desktop/giri_vibration_experiments/iiqs';
tif_dir = '/Users/ryan/Desktop/giri_vibration_experiments/tifs';


tif_files = fullfile(tif_dir, tif_ext);
tifs = dir(tif_files);

raw_files = fullfile(raw_dir, raw_ext);
raws = dir(raw_files);
%%
im_times = NaT(1,numel(raws),'TimeZone','Atlantic/Cape_Verde','Format','yyyy:MM:dd HH:mm:ssZ');
names = {};
for i = 1:numel(raws)
    names{i} = raws(i).name;
    exif_chars = getexif(fullfile(raw_dir,raws(i).name));
   
    exif_struct=struct; 
    lines=split(exif_chars,char(10)); 
    for j=1:length(lines)-1 
    parts=split(lines{j},' :'); 
    exif_struct.(deblank(parts{1}))=deblank(parts{2}); 
    end
    im_times(i) = datetime(exif_struct.FileModifyDate,'TimeZone','Atlantic/Cape_Verde','Format','yyyy:MM:dd HH:mm:ssZ');
end

%%
blobs = {};
for i = 1:numel(tifs)
    img = im2double(imread(fullfile(tif_dir,tifs(i).name)));
    local_thresh = 3;
    % size of the blobs
    rad = 16;
    % how strong of a response defines a blob
    global_thresh = 0.05;

    these_blobs = detect_blobs(img,rad,local_thresh,global_thresh);
    % and bring blobs from image space to xy space
    blobs{i} = [these_blobs(:,2),-these_blobs(:,1)];
end

%%
movements = {};
for i = 2:numel(blobs)
    all_distances = pdist2(blobs{i},blobs{i-1});
    [movements{i},closest_inds] = min(all_distances);
    clear all_distances
end
%%
movements_mat = zeros(max(cellfun(@numel,movements)),numel(movements));

for i = 1:numel(movements)
    movements_mat(1:length(movements{i}),i) = movements{i};
end

%% read in the vibration data
vib_data_folder = '/Users/ryan/Desktop/giri_vibration_experiments/Princeton U Rock Grinder Vibration Data';
vib_file_ext = '*.xlsx';

vib_files = fullfile(vib_data_folder, vib_file_ext);
vib_spreadsheets = dir(vib_files);

vib_data = {};
vib_data_names = {};
start_times = NaT(1,numel(vib_spreadsheets));

for i = 1:numel(vib_spreadsheets)
    opts = detectImportOptions(fullfile(vib_data_folder,vib_spreadsheets(i).name));
    these_data = readmatrix(fullfile(vib_data_folder,vib_spreadsheets(i).name),opts);
    vib_data{i} = these_data(~any(isnan(these_data),2),:);
    vib_data_names{i} = vib_spreadsheets(i).name;
    
    [~,~,date] = xlsread(fullfile(vib_data_folder,vib_spreadsheets(i).name),1,'C81');
    [~,time] = xlsread(fullfile(vib_data_folder,vib_spreadsheets(i).name),1,'C82');
    start_times(i) = datetime(date{1},'ConvertFrom','excel')+ timeofday(datetime(time{1},'InputFormat', 'hh:mm:ss:SSS')) + hours(2);
end  

start_times.TimeZone = im_times.TimeZone;
%%
figure()
count = 1;
for i = [2,6,4,3,1,5]
    plot_times = start_times(i)+seconds(vib_data{i}(:,2));
    subplot(3,2,count)

    p4 = plot(plot_times,vib_data{i}(:,9),'DisplayName','Channel 4');
    p4.Color(4) = 0.8;
    hold on
    p1 = plot(plot_times,vib_data{i}(:,3),'DisplayName','Channel 1');
    p1.Color(4) = 0.8;
    p2 = plot(plot_times,vib_data{i}(:,5),'DisplayName','Channel 2');
    p2.Color(4) = 0.8;
    p3 = plot(plot_times,vib_data{i}(:,7),'DisplayName','Channel 3');
    p3.Color(4) = 0.8;
    axis tight
    title([string(count) + '. ' + vib_data_names{i}])
    ylabel('acceleration(m/s^{2})')
    data_here = isbetween(im_times,plot_times(1),plot_times(end));
    
    if i == 4
        legend()
    end

    
    count = count +1;
end
%%
figure()
scatter(im_times,mean(movements_mat),80,'filled')
ylabel('mean dot movement (pixels)')
text(start_times(2),0.5,'1','FontSize',14)
text(start_times(6),0.5,'2','FontSize',14)
text(start_times(4),0.5,'3','FontSize',14)
text(start_times(3),0.5,'4','FontSize',14)
text(start_times(1),0.5,'5','FontSize',14)
text(start_times(5),0.5,'6','FontSize',14)
xlim([datetime(2021,3,12,7,45,0,'TimeZone',im_times.TimeZone) datetime(2021,3,12,10,45,0,'TimeZone',im_times.TimeZone)])