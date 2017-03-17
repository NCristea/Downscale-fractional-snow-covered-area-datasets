function stats = stats30(naip30msnow, rcomp)

% this function computes binary statistics precision recall and F score for
% binary grids comparisons 
%function call: stats = stats30(binary_map1, binary_map2)

no_pixels_r = size(find(rcomp > 0), 1);


%prcentage of snow cover over the entire domain for different snow maps
%pr_naip = no_sn_pixels_naip/no_pixels;

%estimate common areas
%
com_r = rcomp.* naip30msnow;


%estimate pixels covered by snow in common areas

no_sn_pixels_com = size(find(com_r > 0), 1);

% 

naip_umpapped_r = naip30msnow - com_r;%false negatives 

no_sn_pixels_naip_umpapped_r = size(find(naip_umpapped_r > 0), 1);

dwn_unmapped_r = rcomp - com_r;%false positive

no_sn_pixels_dwn_unmapped_r = size(find(dwn_unmapped_r > 0), 1);%false positive

%calculate precision = true_pos/(true_pos + false_pos)
pr = no_sn_pixels_com/(no_sn_pixels_com + no_sn_pixels_dwn_unmapped_r);
%calculate recall  = true_pos/(true_pos + false_neg)
recall = no_sn_pixels_com/(no_sn_pixels_com + no_sn_pixels_naip_umpapped_r);
% fscore - both
Fscore = (2 * no_sn_pixels_com)/(2 * no_sn_pixels_com + no_sn_pixels_dwn_unmapped_r + no_sn_pixels_naip_umpapped_r);

stats = [pr recall Fscore];
