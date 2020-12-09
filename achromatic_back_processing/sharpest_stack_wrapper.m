% Script to wrap through functions to select each channel's sharpest image
% and attempt to stack them into a single image
%
%
% R. A. Manzuk 11/20/2020
% last updated 12/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set everything up
n_windows = 5;
window_size = 300;
n_zsteps = 15;
z_step_size = 20;

nir_wavelength = 850;
red_wavelength = 625;
green_wavelength = 530;
blue_wavelength = 470;
uv_wavelength = 365;

blue_dir = '/Users/rmanzuk/Desktop/achromatic_tests/labrador_120mm/blue';
green_dir = '/Users/rmanzuk/Desktop/achromatic_tests/labrador_120mm/green';
red_dir = '/Users/rmanzuk/Desktop/achromatic_tests/labrador_120mm/red';
ir_dir = '/Users/rmanzuk/Desktop/achromatic_tests/labrador_120mm/ir';
uv_dir = '/Users/rmanzuk/Desktop/achromatic_tests/labrador_120mm/uv';

%% read all image stacks and get sharpness metrics

[blue_global_sharpness,blue_window_sharpness,edge_coords] = stack_sharpness(blue_dir,n_windows,window_size,n_zsteps,z_step_size);
[green_global_sharpness,green_window_sharpness,~] = stack_sharpness(green_dir,n_windows,window_size,n_zsteps,z_step_size,edge_coords);
[red_global_sharpness,red_window_sharpness,~] = stack_sharpness(red_dir,n_windows,window_size,n_zsteps,z_step_size,edge_coords);
[ir_global_sharpness,ir_window_sharpness,~] = stack_sharpness(ir_dir,n_windows,window_size,n_zsteps,z_step_size,edge_coords);
[uv_global_sharpness,uv_window_sharpness,~] = stack_sharpness(uv_dir,n_windows,window_size,n_zsteps,z_step_size,edge_coords);

%% find and read in the sharpest image for each stack...whole vs windows

[blue_sharpest,blue_sharpest_loc] = select_sharpest(blue_dir, blue_global_sharpness, blue_window_sharpness, 'whole',n_zsteps,z_step_size);
[green_sharpest,green_sharpest_loc] = select_sharpest(green_dir, green_global_sharpness, green_window_sharpness, 'whole',n_zsteps,z_step_size);
[red_sharpest,red_sharpest_loc] = select_sharpest(red_dir, red_global_sharpness, red_window_sharpness, 'whole',n_zsteps,z_step_size);
[ir_sharpest,ir_sharpest_loc] = select_sharpest(ir_dir, ir_global_sharpness, ir_window_sharpness, 'whole',n_zsteps,z_step_size);
[uv_sharpest,uv_sharpest_loc] = select_sharpest(uv_dir, uv_global_sharpness, uv_window_sharpness, 'whole',n_zsteps,z_step_size);

[blue_sharpest_windows,blue_sharpest_loc_windows] = select_sharpest(blue_dir, blue_global_sharpness, blue_window_sharpness, 'windows',n_zsteps,z_step_size);
[green_sharpest_windows,green_sharpest_loc_windows] = select_sharpest(green_dir, green_global_sharpness, green_window_sharpness, 'windows',n_zsteps,z_step_size);
[red_sharpest_windows,red_sharpest_loc_windows] = select_sharpest(red_dir, red_global_sharpness, red_window_sharpness, 'windows',n_zsteps,z_step_size);
[ir_sharpest_windows,ir_sharpest_loc_windows] = select_sharpest(ir_dir, ir_global_sharpness, ir_window_sharpness, 'windows',n_zsteps,z_step_size);
[uv_sharpest_windows,uv_sharpest_loc_windows] = select_sharpest(uv_dir, uv_global_sharpness, uv_window_sharpness, 'windows',n_zsteps,z_step_size);
%% normalize the sharpness curves

[blue_global_normalized,blue_window_normalized] = normalize_sharpness(blue_global_sharpness,blue_window_sharpness);
[green_global_normalized,green_window_normalized] = normalize_sharpness(green_global_sharpness,green_window_sharpness);
[red_global_normalized,red_window_normalized] = normalize_sharpness(red_global_sharpness,red_window_sharpness);
[ir_global_normalized,ir_window_normalized] = normalize_sharpness(ir_global_sharpness,ir_window_sharpness);
[uv_global_normalized,uv_window_normalized] = normalize_sharpness(uv_global_sharpness,uv_window_sharpness);
%%
subplot(5,1,1)
plot(([1:31]-n_zsteps-1).*z_step_size,uv_global_normalized)
hold on
plot(([1:31]-n_zsteps-1).*z_step_size,uv_window_normalized)
ylabel('Sharpness')
title('UV')
axis tight

subplot(5,1,2)
plot(([1:31]-n_zsteps-1).*z_step_size,blue_global_normalized)
hold on
plot(([1:31]-n_zsteps-1).*z_step_size,blue_window_normalized)
ylabel('Sharpness')
title('blue')
axis tight

subplot(5,1,3)
plot(([1:31]-n_zsteps-1).*z_step_size,green_global_normalized)
hold on
plot(([1:31]-n_zsteps-1).*z_step_size,green_window_normalized)
ylabel('Sharpness')
title('green')
axis tight

subplot(5,1,4)
plot(([1:31]-n_zsteps-1).*z_step_size,red_global_normalized)
hold on
plot(([1:31]-n_zsteps-1).*z_step_size,red_window_normalized)
ylabel('Sharpness')
title('red')
axis tight

subplot(5,1,5)
plot(([1:31]-n_zsteps-1).*z_step_size,ir_global_normalized)
hold on
plot(([1:31]-n_zsteps-1).*z_step_size,ir_window_normalized)
ylabel('Sharpness')
title('IR')
axis tight



%% for the table height vs wavelength plot
sharpness_stretch = 40;


xx = [350:900];
yy = spline([uv_wavelength,blue_wavelength,green_wavelength,red_wavelength,nir_wavelength],[uv_sharpest_loc,blue_sharpest_loc,green_sharpest_loc,red_sharpest_loc,ir_sharpest_loc],xx);
subplot(2,1,1)
plot([uv_wavelength,blue_wavelength,green_wavelength,red_wavelength,nir_wavelength],[uv_sharpest_loc,blue_sharpest_loc,green_sharpest_loc,red_sharpest_loc,ir_sharpest_loc],'o',xx,yy)
hold on
plot((uv_global_normalized.*sharpness_stretch)+uv_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((blue_global_normalized.*sharpness_stretch)+blue_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((green_global_normalized.*sharpness_stretch)+green_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((red_global_normalized.*sharpness_stretch)+red_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((ir_global_normalized.*sharpness_stretch)+nir_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
legend('actual data','spline fit','sharpness curves')
xlabel('wavelength [nm]')
ylabel('focus table height [um]')
title('whole image')
axis tight

xx = [350:900];
yy = spline([uv_wavelength,blue_wavelength,green_wavelength,red_wavelength,nir_wavelength],[uv_sharpest_loc_windows,blue_sharpest_loc_windows,green_sharpest_loc_windows,red_sharpest_loc_windows,ir_sharpest_loc_windows],xx);
subplot(2,1,2)
plot([uv_wavelength,blue_wavelength,green_wavelength,red_wavelength,nir_wavelength],[uv_sharpest_loc_windows,blue_sharpest_loc_windows,green_sharpest_loc_windows,red_sharpest_loc_windows,ir_sharpest_loc_windows],'o',xx,yy)
hold on
plot((uv_window_normalized.*sharpness_stretch)+uv_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((blue_window_normalized.*sharpness_stretch)+blue_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((green_window_normalized.*sharpness_stretch)+green_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((red_window_normalized.*sharpness_stretch)+red_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((ir_window_normalized.*sharpness_stretch)+nir_wavelength,([1:31]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
legend('actual data','spline fit', 'sharpness curves')
xlabel('wavelength [nm]')
ylabel('focus table height [um]')
title('edge windows')
axis tight