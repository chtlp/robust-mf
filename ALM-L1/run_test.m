a1=0;a2=0;mu=10;

x = -5:0.2:5;
y = arrayfun(@(d) alm_l1_loss(d, a1, a2, mu), x);
plot(x, y);
desc = sprintf('a1=%.1f, a2=%.1f, mu=%.1f.fig', a1, a2, mu);
saveas(gcf, desc, 'fig');

%%
a1=0;a2=0;mu=1;

x = -5:0.2:5;
y = arrayfun(@(d) alm_l1_loss(d, a1, a2, mu), x);
plot(x, y);
desc = sprintf('a1=%.1f, a2=%.1f, mu=%.1f.fig', a1, a2, mu);
saveas(gcf, desc, 'fig');

%%
a1=0;a2=10;mu=10;

x = 0:0.2:10;
y = arrayfun(@(d) alm_l1_loss(d, a1, a2, mu), x);
plot(x, y);
desc = sprintf('a1=%.1f, a2=%.1f, mu=%.1f.fig', a1, a2, mu);
saveas(gcf, desc, 'fig');

%%
a1=0;a2=10;mu=1;

x = 0:0.2:10;
y = arrayfun(@(d) alm_l1_loss(d, a1, a2, mu), x);
plot(x, y);
desc = sprintf('a1=%.1f, a2=%.1f, mu=%.1f.fig', a1, a2, mu);
saveas(gcf, desc, 'fig');



