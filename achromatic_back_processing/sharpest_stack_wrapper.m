% Script to wrap through functions to select each channel's sharpest image
% and attempt to stack them into a single image
%
%
% R. A. Manzuk 11/20/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set everything up
n_windows = 10;
window_size = 1000;
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

[blue_global_sharpness,blue_window_sharpness] = stack_sharpness(blue_dir,n_windows,window_size,n_zsteps,z_step_size);
[green_global_sharpness,green_window_sharpness] = stack_sharpness(green_dir,n_windows,window_size,n_zsteps,z_step_size);
[red_global_sharpness,red_window_sharpness] = stack_sharpness(red_dir,n_windows,window_size,n_zsteps,z_step_size);
[ir_global_sharpness,ir_window_sharpness] = stack_sharpness(ir_dir,n_windows,window_size,n_zsteps,z_step_size);
[uv_global_sharpness,uv_window_sharpness] = stack_sharpness(uv_dir,n_windows,window_size,n_zsteps,z_step_size);

%% find and read in the sharpest image for each stack

[blue_sharpest,blue_sharpest_loc] = select_sharpest(blue_dir, blue_global_sharpness, blue_window_sharpness, 'both',n_zsteps,z_step_size);
[green_sharpest,green_sharpest_loc] = select_sharpest(green_dir, green_global_sharpness, green_window_sharpness, 'both',n_zsteps,z_step_size);
[red_sharpest,red_sharpest_loc] = select_sharpest(red_dir, red_global_sharpness, red_window_sharpness, 'both',n_zsteps,z_step_size);
[ir_sharpest,ir_sharpest_loc] = select_sharpest(ir_dir, ir_global_sharpness, ir_window_sharpness, 'both',n_zsteps,z_step_size);
[uv_sharpest,uv_sharpest_loc] = select_sharpest(uv_dir, uv_global_sharpness, uv_window_sharpness, 'both',n_zsteps,z_step_size);


%%
whole_im = cat(3,red_sharpest,green_sharpest,blue_sharpest);
%%
subplot(2,1,1)
plot(([1:31]-n_zsteps).*z_step_size,uv_global_sharpness,'DisplayName','global sharpness')
ylabel('Sharpness')
title('UV whole image')
axis tight
subplot(2,1,2)
plot(([1:31]-n_zsteps).*z_step_size,uv_window_sharpness,'LineStyle',':')
xlabel('location [um]')
ylabel('Sharpness')
title('UV windows')
axis tight

