function [r,p] = correl(data1,data2,type)
% function [r,p] = correl(data1,data2)
% just a wrapper for corr.m, handles NaNs (excludes them from both
% vectors it has to correlate) and returns single R and P values, not
% matrices. Works with pairs of columns on matrices.

if size(data1) ~= size(data2)
    error('2 inputs not same size!')
end

if find([size(data1)]==1)
    data1 = data1(:);
    data2 = data2(:);
end

if ~exist('type','var')
    type = 'Pearson';
end

for c = 1:size(data1,2)
    idx = intersect(find(~isnan(data1(:,c))),find(~isnan(data2(:,c))));
%     idx = ~isnan(data1(:,c)) & ~isnan(data2(:,c));
    [r(c),p(c)] = corr(data1(idx,c),data2(idx,c),'type',type);
%     r(c) = temp1(1,2); % the correlation coefficient "R"
%     p(c) = temp2(1,2); % the associted p-value
end
