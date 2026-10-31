function matrix_circle_scale(X_values,X_sizes,rownames,colnames,colormap,fontsize)

% matrix_circle_scale(X_values,X_sizes,rownames,colnames,colormap,fontsize)
%
% plot values of matrix as dots: X_sizes are shown by area (default: equal sizes),
% X_values by colour
%
% X_values are assumed to be between -1 and 1; if there are bigger
% absolute values, all values are normalised to the maximal absolute value
%
% USAGE EXAMPLE: matrix_circle_scale(randn(2),rand(2),{'a1','a2'},{'b1','b2'},'jet')
%

if ~exist('rownames','var'), rownames=[]; end
if ~exist('colnames','var'), colnames=[]; end
if ~exist('colormap','var'), colormap='rb_colors'; end
if ~exist('fontsize','var'), fontsize = 10; end
if isempty(X_sizes), X_sizes = ones(size(X_values)); end

cm = eval(colormap);
X_sizes = sqrt(abs(fliplr(X_sizes')));  % area proportional to absolute value
nX = 0.5/max(max(X_sizes));
X_values = fliplr(X_values');

% normalise to maximum absolute value
if max(abs(X_values(:)))>1
  display('Normalising the values');
  X_values = X_values/max(abs(X_values(:)));
end

X_values = 0.5 * [X_values + max(abs(X_values(:)))];

for it1=1:size(X_sizes,1)
  for it2=1:size(X_sizes,2)
    ind_color = 1+ceil([size(cm,1)-1]*X_values(it1,it2));
    if isfinite(ind_color)
      circle(it1,it2,nX * X_sizes(it1,it2),cm(ind_color,:)); hold on;
    end
  end
end

hold off;
axis([0.5 size(X_sizes,1) + 0.5 0.5  size(X_sizes,2)+0.5])

if isempty(rownames), set(gca,'YTick',[]);
else, set(gca,'YTick',1:size(X_sizes,2)); set(gca,'YTickLabel',flipud(rownames),'Fontsize',fontsize);
end

if isempty(colnames), set(gca,'XTick',[]);
else, set(gca,'XTick',1:size(X_sizes,1)); set(gca,'XTickLabel',colnames);
  my_xticklabel([],size(X_values,2)+1,[],fontsize);
end

axis equal; axis square;  %axis tight;

function circle(x,y,l,col)

  fill(x + l * cos(0:0.2:2*pi), y + l * sin(0:0.2:2*pi),col,'EdgeColor',[.7 .7 .7])
