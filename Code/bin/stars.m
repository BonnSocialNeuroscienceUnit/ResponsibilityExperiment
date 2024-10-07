function out = stars(p)

for i = 1:length(p)
  if p(i) < 0.001; out{i} = '***';
  elseif p(i) < 0.01; out{i} = '**';
  elseif p(i) < 0.05; out{i} = '*';
  elseif p(i) < 0.1; out{i} = '^';
  else out{i} = ' ';
  end
end
