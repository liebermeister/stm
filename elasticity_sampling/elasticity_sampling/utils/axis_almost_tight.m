function axis_almost_tight()

axis tight
a = axis;

if a(3)<=0,
  axis([a(1)-0.05*(a(2)-a(1)), a(2)+0.05*(a(2)-a(1)),a(3)-0.05*(a(4)-a(3)), a(4)+0.05*(a(4)-a(3))])
else
  axis([a(1)-0.05*(a(2)-a(1)), a(2)+0.05*(a(2)-a(1)),exp(log(a(3))-0.05*(log(a(4))-log(a(3)))), exp(log(a(4))+0.05*(log(a(4))-log(a(3))))])
end