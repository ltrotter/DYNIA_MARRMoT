function sample = lhs_sample_par(model, n)

ranges = model.parRanges;
lhs = lhsdesign(n,size(ranges, 1));
sample = lhs'.*diff(ranges,1,2) + ranges(:,1);