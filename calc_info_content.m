function [norm_info_content, info_content] = calc_info_content(cd,x_sorted,q)

if nargin < 3 || isempty(q)
    q = .9;
end

q_lo = (1-q)/2;
q_hi = 1-q_lo;

[~,idx_lo] = min(abs(q_lo - cd),[],2);
idx_lo = squeeze(idx_lo);

[~,idx_hi] = min(abs(q_hi - cd),[],2);
idx_hi = squeeze(idx_hi);

idx1 = reshape(repmat(1:size(x_sorted,1),1,size(x_sorted,3)), size(idx_lo));
idx3 = reshape(repelem(1:size(x_sorted,3),1,size(x_sorted,1)),size(idx_lo));

lower = x_sorted(sub2ind(size(x_sorted),idx1,idx_lo,idx3));
upper = x_sorted(sub2ind(size(x_sorted),idx1,idx_hi,idx3));

CI = upper-lower;

ranges = squeeze(max(x_sorted,[],[1,2]) - min(x_sorted,[],[1,2]));

info_content = 1 - CI./ranges';
norm_info_content = info_content./max(info_content,[],1);