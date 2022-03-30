function sample = unif_sample_par(model, n)

ranges = model.parRanges;
rd = rand(size(ranges, 1),n);
sample = rd.*diff(ranges,1,2) + ranges(:,1);