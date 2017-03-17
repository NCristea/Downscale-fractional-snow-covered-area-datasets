%this script is reading the orthophoto-derived 500-m fractional snow cover
%area data for August 16, 2011 across Mt. Rainier, WA, to generate binary
%(presence/absence)snow data at 30-m resolution

%load previously derived DAH and TPI grids (see how in the document
%attached) 
load dah30m_30_mt_rainier;
load tpi_30m_2_mt_rainier; % this used 60m TPI searching distance

%load binary (snow/no snow)validation data 
load naip30msnow;

%load fSCA tiff grid and get the coarse grid X and Y coords
[scene_SCA, cmap, R, ~] = geotiffread('fsca_rec_228_11_res4');
rows = size(scene_SCA,1);
cols = size(scene_SCA,2);
[scene_500_X, scene_500_Y]= refmat2meshgrid(R, rows, cols);
fsca_500 = double(scene_SCA);

%plot the fSCA scene 
pcolor(scene_500_X, scene_500_Y, fsca_500);colorbar;

%load DEM tiff grid and get the high resolution grid X and Y coords
[dem, cmap_dem, R_dem, ~] = geotiffread('dem30_mt_rainier');
rows_dem = size(dem,1);
cols_dem = size(dem,2);
[dem_X, dem_Y]= refmat2meshgrid(R_dem, rows_dem, cols_dem);
%%
%call the downscale function
%inputs: X(coarse scale grid), Y(coarse scale grid), fSCA(coarse scale
%grid in percetange 0-100), X(fine scale grid), Y(fine scale grid), DAH
%grid (fine scale), TPI grid(fine scale), weight value

r_composite =  downscale_composite_index(scene_500_X, scene_500_Y, fsca_500.*100, dem_X, dem_Y, scene_dah, scene_tpi, 0.2);

%plot the downscaled map

figure; h = pcolor(dem_X, dem_Y, r_composite);colorbar;
set(h, 'EdgeColor', 'none');

%calculate statistics 

stats = stats30(naip30msnow, r_composite)