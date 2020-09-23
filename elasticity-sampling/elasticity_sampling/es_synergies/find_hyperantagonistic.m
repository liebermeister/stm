function find_hyperantagonistic(fx,Fxx)

% find_hyperantagonistic(fx,Fxx)
%
% check a given fitness derivative vector and curvature matrix
% for hyperantagonistic responses to gene decrease

nr = length(fx);

for i1 = 1:nr,
  for i2 = 1:i1-1,
    fa  = fx(i1);
    fb  = fx(i2);
    faa = Fxx(i1,i1);
    fab = Fxx(i1,i2);
    fbb = Fxx(i2,i2);
    [success,solution] = check_hyperantagonistic(fa,fb,faa,fab,fbb);
    display(sprintf('Enzyme pair %d/%d; success: %d',i1,i2,success));
  end
end
