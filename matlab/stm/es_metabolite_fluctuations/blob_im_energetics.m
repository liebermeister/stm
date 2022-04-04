function blob_im_energetics(CJ_sc,my_xlabel,my_ylabel)

maxval = max(max(abs(CJ_sc)));
matrix_circle_scale(CJ_sc/maxval,CJ_sc/maxval,[],[],'rb_colors',20)
xlabel(my_xlabel); ylabel(my_ylabel); colorbar
axis equal; colormap(rb_colors); 
set(gca,'Clim',maxval*[-1,1]); set(gca,'FontSize',24);
colorbar
