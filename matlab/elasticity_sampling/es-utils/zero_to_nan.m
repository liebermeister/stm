function MM = zero_to_nan(M,epsilon)

% MM = zero_to_nan(M,epsilon)

MM=M;
if exist('epsilon','var'),
  MM(abs(MM)<=epsilon) = nan;
else
  MM(abs(MM)==0) = nan;
end