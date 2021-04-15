function es_synergy_ensemble_comparison_graphics(network_CoHid, res, condition1, condition2, print_flag, filenames, result_file)

% es_synergy_ensemble_comparison_graphics(network_CoHid, res, condition1, condition2, print_flag, filenames, result_file)
%
% Graphics: Compare synergy values between two model ensembles - show statistically significant differences
  
% --------------------------------------------------------------
% graphics 

% cd(filenames.data_dir); load(result_file);

col = [0.8 0.8 0.8; rb_colors];

figure(1);
subplot(2,2,1); im(res.p_interaction.mean_a,[]); colormap(col); axis equal; axis off;  title(['Mean ' condition1 ]);
subplot(2,2,3); im(res.p_interaction.mean_a_significant,2); colormap(col);  axis equal; axis off;  colorbar off; title(['Significance ' condition1])
subplot(2,2,2); im(res.p_interaction.mean_b); colormap(col);  axis equal; axis off;  title(['Mean ' condition2 ])
subplot(2,2,4); im(res.p_interaction.mean_b_significant,2); colormap(col);  axis equal; axis off;   colorbar off; title(['Significance ' condition2])

figure(2);
subplot(2,2,1); im(res.p_interaction.mean_total); colormap(col); axis equal; axis off;  title('Mean')
subplot(2,2,3); im(res.p_interaction.mean_total_significant,2); colormap(col); axis equal; axis off;  colorbar off; title('Significance: Mean > Overall mean')
subplot(2,2,2); im(res.p_interaction.mean_delta); colormap(col); axis equal; axis off;  title('Difference')
subplot(2,2,4); im(res.p_interaction.mean_delta_significant,2); colormap(col); axis equal; axis off;  colorbar off; title('Significance: Difference > 0')

% ----------------------------------------------------------
% network graphics

col = rb_colors(20); col = col(4:17,:);

% mean values (all values are shown)

figure(3); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_total,res.p_interaction.mean_total, col,struct('relative_threshold',0.1)); 

figure(4); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_delta,res.p_interaction.mean_delta, col,struct('relative_threshold',0.1)); 

figure(5); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_a,res.p_interaction.mean_a, col,struct('relative_threshold',0.1)); 

figure(6); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_b,res.p_interaction.mean_b, col,struct('relative_threshold',0.1)); 


% mean values (only significant values are shown)

figure(7); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_total_significant,res.p_interaction.mean_total_significant, col, struct('relative_threshold',0.1)); 

figure(8); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_delta_significant,res.p_interaction.mean_delta_significant, col, struct('relative_threshold',0.1)); 

figure(9); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_a_significant,res.p_interaction.mean_a_significant, col, struct('relative_threshold',0.1)); 

figure(10); clf; 
interaction_network_plot(network_CoHid, res.p_influence.mean_b_significant,res.p_interaction.mean_b_significant, col, struct('relative_threshold',0.1)); 


if print_flag,
  if strcmp(result_file(end-3:end),'.mat'), result_file = result_file(1:end-4); end
  cd(filenames.psfile_dir);
  print([result_file, '_significance_conditions.eps'],'-f1','-depsc');
  print([result_file, '_significance_total_delta.eps'],'-f2','-depsc');
  print([result_file, '_mean_total.eps'],'-f3','-depsc');
  print([result_file, '_mean_delta.eps'],'-f4','-depsc');
  print([result_file, '_mean_' condition1 '.eps'],'-f5','-depsc');
  print([result_file, '_mean_' condition2 '.eps'],'-f6','-depsc');
  print([result_file, '_mean_significant_total.eps'],'-f7','-depsc');
  print([result_file, '_mean_significant_delta.eps'],'-f8','-depsc');
  print([result_file, '_mean_significant_' condition1 '.eps'],'-f9','-depsc');
  print([result_file, '_mean_significant_' condition2 '.eps'],'-f10','-depsc');
end

figure(3); title('Mean influences / interactions total')
figure(4); title(['Difference influences / interactions'])
figure(5); title(['Mean influences / interactions ' condition1])
figure(6); title(['Mean influences / interactions ' condition2])
figure(7); title('Mean influences / interactions total')
figure(8); title(['Difference influences / interactions'])
figure(9); title(['Mean influences / interactions ' condition1])
figure(10); title(['Mean influences / interactions ' condition2])
