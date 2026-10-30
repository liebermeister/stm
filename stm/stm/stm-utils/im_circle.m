function im_circle(X,range,rownames,colnames,color)

% im_circle(X,range,rownames,colnames,color)
%
% Show the positive entries of a matrix X as filled circles (area proportional to the value)

if ~exist('range','var'),  range = []; end
if ~exist('rownames','var'), rownames=[]; end
if ~exist('colnames','var'), colnames=[]; end
if ~exist('color','var'), color='k'; end

if isempty(range), range = max(max(abs(X))); end

X   = fliplr(X');
sgn = sign(fliplr(X'));
nX  = 1/range;

for it1=1:size(X,1)
  for it2=1:size(X,2)
    if X(it1,it2)>0, circle(it1,it2, nX * X(it1,it2),color); hold on; end
    if X(it1,it2)<0, circle(it1,it2,-nX * X(it1,it2),'w'); hold on; end
  end
end
hold off;

axis equal;
axis([0.5 size(X,1)+0.5 0.5  size(X,2)+0.5])

if isempty(rownames), set(gca,'YTick',[]);
else, set(gca,'YTick',1:size(X,2)); set(gca,'YTickLabel',flipud(rownames));
end

if isempty(colnames), set(gca,'XTick',[]);
else, set(gca,'XTick',1:size(X,1)); set(gca,'XTickLabel',colnames);
end

% ---------------------------------------

function circle(x,y,l,col)

  % l: value to be shown -> corresponds to area of the circle

  fill(x + 0.5*sqrt(l) * cos(0:0.1:2*pi), y + 0.5*sqrt(l) * sin(0:0.1:2*pi),col)
