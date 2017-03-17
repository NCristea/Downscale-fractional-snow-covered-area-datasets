function [ result ] = downscale_composite_index( low_res_x, low_res_y, low_res, hi_res_x, hi_res_y, hi_res, c_res, wei)
%downscale 
%   This function is used to downscale coarse resoltion fractional snow
%   cover (fSCA) products to binary (presence/absence) high resoltion snow data
% To call the function: 
% r = downscale_composite_index( low_res_x, low_res_y, low_res, hi_res_x, hi_res_y, hi_res, c_res, wei)
%  low_res_x = coarse resolution X fSCA grid in UTM coordinates
%  low_res_y = coarse resolution Y fSCA grid in UTM coordinates
%  low_res = coarse resolution fSCA grid 
%  hi_res_x = high resolution X grid in UTM coordinates
%  hi_res_y = high resolution Y grid in UTM coordinates
%  hi_res =  high resolution diurnal anisotropic heating index (DAH)
%  c_res = high resolution topographic position index (TPI)


[low_res_x_size, low_res_y_size] = size(low_res);
% x step for low resolution
low_res_dx = low_res_x(low_res_x_size, 2) - low_res_x(low_res_x_size, 1);
% y step for low resolution
low_res_dy = low_res_y(1,2) - low_res_y(2,2);

[high_res_x_length high_res_y_length] = size(hi_res);
% x step for high resolution
high_res_dx = hi_res_x(high_res_x_length, 2) - hi_res_x(high_res_x_length, 1);
% y step for high resolution
high_res_dy = hi_res_y(1,2) - hi_res_y(2,2);
result = zeros(size(hi_res));

% map factors (one for x and one for y) for from low resolution to high resolution
x_map_factor = low_res_dx / high_res_dx;
y_map_factor = low_res_dy / high_res_dy;
tmp_array = [];

for i = 1:low_res_x_size
    for j = 1:low_res_y_size
        snow_percentage = double(low_res(i, j))/100;
        if snow_percentage <= 0
            %fprintf('No snow for cell <%d, %d>, skipping...\n', i, j);
            continue
        end
        % compute the range of high resolution cells that map to current
        % low resolution cell 
        x_range = high_res_range(i, x_map_factor);
        y_range = high_res_range(j, y_map_factor);
        composite_index = compute_composite_index(double(hi_res(x_range, y_range)), wei, c_res(x_range, y_range), (1-wei));
        % for each high resolution cell store in tmp_array
        % the cell indexes (x, y) and the high resolution data
        for hires_x = 1:(1 + max(x_range) - min(x_range))
            for hires_y = 1:(1 + max(y_range) - min(y_range))
                hi = composite_index(hires_x, hires_y);
                p = [hi (min(x_range)+hires_x-1) (min(y_range)+hires_y-1)];
                tmp_array = [tmp_array; p];
            end
        end
        
        % sort the high resolution cell by the composite index
        sorted_array = sortrows(tmp_array, 1);
        
        % tmp_array content erased here since it's no longer needed
        % until next next iteration (cell that has snow_percentace > 0)
        tmp_array = []; 
        
        area_factor = 1 / (x_map_factor * y_map_factor);
        q = 0;
        r = (1 + max(x_range) - min(x_range))*(1 + max(y_range) - min(y_range));
        
        % allocate snow: while there's still snow remaining (based on
        % snow_percentage) then add snow to the next high resolution cell
        % Here 'next' is given by the sorting done above.
        
%         
        last_index = round(r * snow_percentage);
        last_index1 = int32(last_index);
%         
%         if snow_percentage > 0.95
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.95 && snow_percentage > 0.8
%             r = 220; %r = 220;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.8 && snow_percentage > 0.7
%             r = 170;%r = 220;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.7 && snow_percentage > 0.6
%             r = 120;%r = 220;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.6 && snow_percentage > 0.5
%             r = 110;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.5 && snow_percentage > 0.4
%             r = 100;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         elseif snow_percentage <= 0.4 && snow_percentage > 0.3
%             r = 89;%89;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         else 
%             r = 64;%64;
%             last_index = round(r * snow_percentage);
%             last_index1 = int32(last_index);
%         end

        a = size(sorted_array);
        if isnan(last_index)
            for idx = 1:a(1, 1)
                hires_x = sorted_array(idx, 2);
                hires_y = sorted_array(idx, 3);
                result(hires_x, hires_y) = NaN;
            end
        end
             
        for idx = 1:last_index1
            hires_x = sorted_array(idx, 2);
            hires_y = sorted_array(idx, 3);
            result(hires_x, hires_y) = 1;
        end
        sorted_array = [];
    end
end

end

function composite_index = compute_composite_index(heat_index, heat_weight, curvature_index, curvature_weight)
    % nomalize hi and cv matrices
    mmin_hi = min(heat_index(:));
    mmax_hi = max(heat_index(:));
    % first subtract mmin to have [0; (mmax-mmin)], then normalize by highest value
    hi30norm = (heat_index-mmin_hi) ./ (mmax_hi - mmin_hi); 
    
    mmin_cv = min(curvature_index(:));
    mmax_cv = max(curvature_index(:));
    % first subtract mmin to have [0; (mmax-mmin)], then normalize by highest value
    cv30norm = double((curvature_index-mmin_cv) ./ (mmax_cv-mmin_cv));
    
    %composite_index = heat_index * curvature_index;
    %composite_index = hi30norm .* cv30norm;
    composite_index = heat_weight * hi30norm + curvature_weight * cv30norm;
    %composite_index = heat_weight * heat_index + curvature_weight * curvature_index;
end


