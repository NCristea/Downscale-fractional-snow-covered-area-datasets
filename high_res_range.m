function [ output_args ] = high_res_range(i, map_factor)
    %high_res_range Returns the range of high resolution indexes given a low resolution cell
    %distance   = floor(map_factor);
    distance   = floor(map_factor);
    range_low  = 1 + floor((i - 1) * map_factor);
    range_high = range_low + distance;
    %fprintf('High res range (%d %g): [%d, %d]\n', i, map_factor, range_low, range_high);
    output_args = [ range_low : range_high ];
end