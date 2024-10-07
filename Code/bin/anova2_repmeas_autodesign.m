function [out,str] = anova2_repmeas_autodesign(data,factNames,printOutFlag)
% function out = anova2_repmeas_autodesign(data,factNames,printOutFlag)
% computes 2-way repeated measures ANOVA on data matrix: 1st dim = subjects,
% 2nd dim = factor 1, 3nd dim = factor 2.
%
% Johannes Schultz spring 2007

if ~exist('factNames','var')
  for i = 1:length(size(data))-1
    factNames{i} = ['Factor ' num2str(i)];
  end
end

if ~exist('printOutFlag','var')
  printOutFlag = 1;
end

dimSizes = size(data);

subjects = repmat([1:dimSizes(1)]',dimSizes(2)*dimSizes(3),1);
factor1 = kron(repmat([1:dimSizes(2)]',dimSizes(3),1),ones(dimSizes(1),1));
factor2 = kron([1:dimSizes(3)]',ones(dimSizes(1)*dimSizes(2),1));

table = anova2_repmeas(data(:),subjects,factor1,factor2,factNames);

out.table = table;
out.factorNames = factNames;
out.dataVector = data(:);
out.subjectVector = subjects;
out.factor1Vector = factor1;
out.factor2Vector = factor2;

if nargout == 0 | nargout == 2
  tests = factNames; tests{3} = 'Interaction';
  for c = 2:4
    str{c-1} = sprintf('%s: F(%d,%d)=%.2f, p=%.4f',tests{c-1},out.table{c,3},out.table{c+3,3},out.table{c,5},out.table{c,6});
  end
  str = char(str{:});
  if printOutFlag
    disp(str)
  end
end
