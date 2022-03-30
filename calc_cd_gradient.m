function gradients = calc_cd_gradient(cd,x_sorted,bins)

if nargin < 3 || isempty(bins)
    bins = min(size(x_sorted,2)/2, 20);
end

min_x = squeeze(min(x_sorted,[],[1,2]));
max_x = squeeze(max(x_sorted,[],[1,2]));

dx = (max_x - min_x)/bins;

breaks = min_x + (dx .* repmat(0:bins,size(dx,1),1));

gradients = zeros(size(cd,1), bins, size(cd,3));

for p=1:size(dx,1)
    this_p_breaks = breaks(p,:);
    this_p_x = x_sorted(:,:,p);
    this_p_cd_breaks = zeros(size(x_sorted,1), bins+1);
    for b=1:bins+1
        [~,idx_break] = min(abs(this_p_x - this_p_breaks(b)),[],2);
        this_p_cd_breaks(:,b) = cd(sub2ind(size(cd),1:size(cd,1),idx_break',repelem(p,size(cd,1))));
    end
    gradients(:,:,p) = diff(this_p_cd_breaks,1,2)/dx(p);
end
