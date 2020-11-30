%% plot gamma distribution
%
%

mu_set=20;
sd_set=0.1:0.1:5;
sd_r_set=abs(((sd_set-min(sd_set))/(max(sd_set)-min(sd_set)))*(.635-1)+1);
sd_g_set=abs(((sd_set-min(sd_set))/(max(sd_set)-min(sd_set)))*(.078-.893)+.893);
sd_b_set=abs(((sd_set-min(sd_set))/(max(sd_set)-min(sd_set)))*(.2-.85)+.85);
colormap0=[sd_r_set;sd_g_set;sd_b_set];
for mu_id=1:length(mu_set)
    for sd_id=1:length(sd_set)
        param_gamma=[(mu_set(mu_id)/sd_set(sd_id)).^2,mu_set(mu_id)/((sd_set(sd_id)).^2)];
        f=@(x)(gampdf(x,param_gamma(1),1/param_gamma(2)));
        fplot(f,[0,40],'color',[sd_r_set(sd_id),sd_g_set(sd_id),sd_b_set(sd_id)],"LineWidth",1.5)
        %fplot(f,[0,40])
        hold on
    end
end
colormap(colormap0');
colorbar;