%% Compare One-Sample Tests
% Tianye Ma, Nov 23
% 

%% Prologue
clear
close all
BMWdir='/home/tianye/Matlab_WS/BMW/Bayesian_Modeling_of_Working_Memory';
addpath(BMWdir)
BMW('AddPath');
savedir='/home/tianye/Matlab_WS/NHST_MC';
p_type=3;
% gaussian=1, von mises=2, gamma=3;
% Note that for gamma dist., the mean parameter should never be zero.
if isempty(gcp('nocreate'))
    parpool('local',16); % start parallel computing (use 16 threads)
end
%% Simulation
mu_set=[0, 0.5, 1, 1.5, 2]+20; % mean diff
sd_set=0.1:0.1:5; %
kappa_set=sd2kappa(sd_set);
N_set=4:4:200;
Nsim=100;
%dataset=cell(Nsim,length(N_set),length(sd_set));
dataset=cell(Nsim,1);
Nsim_left=Nsim;
fprintf('\nNow start simulation...\n')
sim_min_id=1;
for i=1:10
    Nsim_1=min(Nsim_left,ceil(Nsim/10));
    for sim_id=sim_min_id:sim_min_id+Nsim_1-1
        dataset{sim_id}=datagenerate(mu_set,sd_set,kappa_set,N_set,p_type);
    end
    Nsim_left=Nsim_left-Nsim_1;
    sim_min_id=sim_min_id+Nsim_1;
    fprintf('%d simulations left.\n',Nsim_left)
end
fprintf('\nSimulation finished...\n')

%% NHSTs
% method_set={'signtest','mult_comp_perm_t1','ttest'};
method_set={'mult_comp_perm_t1','ttest'};
Nmethod=length(method_set);
method=cell(Nmethod,1);
for i=1:Nmethod
    method{i}=str2func(method_set{i});  
end
NHST=cell(length(mu_set),1);
fprintf('\nNow start doing NHSTs...\n')
% mean difference test
parfor mu_id=1:length(mu_set)
    NHST{mu_id}=dotest(N_set,sd_set,Nmethod,method_set,method,dataset,Nsim,mu_id);
end
fprintf('\nNHSTs finished...\n')
% normality test
NHST_norm=cell(length(mu_set),1);
fprintf('\nNow start doing normality tests...\n')
parfor mu_id=1:length(mu_set)
    NHST_norm{mu_id}=donormtest(dataset,length(N_set),length(sd_set),Nsim,mu_id);
end
fprintf('\nNormality tests finished...\n')

%% Comparison
fprintf('\nVisualization...\n')
testdiff=cell(1,length(mu_set));
if Nmethod==2
    Nplot=3;
else
    Nplot=2;
end
% set color scale
map_color=[[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(255,0,256-23)' ]]/255;
% heatmaps of rejection rate
for j=1:Nmethod
    for i=1:length(mu_set)
        subplot(Nplot,length(mu_set),(j-1)*length(mu_set)+i);
        testdata=NHST{i};
        imagesc(reshape(testdata.Prej(j,:,:),[length(N_set),length(sd_set)]))
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        caxis([-.1,1]);
        colormap(map_color);
    end
end
% heatmaps of rejection rate difference
if Nmethod==2
    for i=1:length(mu_set)
        testdata=NHST{i};
        testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(sd_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(sd_set)]);
        subplot(Nplot,length(mu_set),i+Nmethod*length(mu_set))
        imagesc(testdiff)
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        caxis([-.1,1]);
        colormap(map_color);
    end
end
% sd
figure(2)
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(sd_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(sd_set)]);
    plot(mean(testdiff))
    hold on
end
% ss
figure(3)
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(sd_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(sd_set)]);
    plot(mean(testdiff,2))
    hold on
end
% normality (mean)
figure(4)
norm_sum=zeros(length(N_set),length(sd_set));
for i=1:length(mu_set)
    norm_sum=norm_sum+NHST_norm{i}.Prej(:,:);
end
imagesc(norm_sum/length(mu_set))
% kurtosis
figure(5)
kurtosis_sum=zeros(length(N_set),length(sd_set));
for i=1:length(mu_set)
    for j=1:length(N_set)
        for k=1:length(sd_set)
            for l=1:Nsim
                kurtosis_sum(j,k)=kurtosis_sum(j,k)+kurtosis(dataset{l}{j,i,k});
            end
        end
    end
end
imagesc(abs(kurtosis_sum/length(mu_set)/Nsim-3))
% %% Epilogue
% rmpath(BMWdir)
% fprintf('\nSave all data...\n')
% cd(savedir)
% savename=sprintf('type%dpermutation',p_type);
% save(savename)
% fprintf('\nAll done!\n')

%%
function dataset=datagenerate(mu_set,sd_set,kappa_set,N_set,p_type)
dataset=cell(length(N_set),length(mu_set),length(sd_set));
for mu_id=1:length(mu_set)
    for sd_id=1:length(sd_set)
        for N_id=1:length(N_set)
            switch p_type
                case 1
                    param_n=[mu_set(mu_id),sd_set(sd_id)];
                    dataset0=param_n(1)+param_n(2)*randn([N_set(N_id),1]);
                case 2
                    param_vm=[mu_set(mu_id),kappa_set(sd_id)];
                    dataset0=180/pi*vmrand(param_vm(1)/180*pi,param_vm(2),[N_set(N_id),1]);
                case 3
                    param_gamma=[(mu_set(mu_id)/sd_set(sd_id)).^2,mu_set(mu_id)/((sd_set(sd_id)).^2)];
                    dataset0=gamrnd(param_gamma(1),1/param_gamma(2),[N_set(N_id),1]);
            end
            dataset{N_id,mu_id,sd_id}=dataset0;
        end
    end
end
end

%%
function NHST=dotest(N_set,sd_set,Nmethod,method_set,method,dataset,Nsim,mu_id)
NHST.Prej=zeros(Nmethod, length(N_set), length(sd_set));
NHST.pvalue=zeros(Nmethod, length(N_set), length(sd_set));
for N_id=1:length(N_set)
    for sd_id=1:length(sd_set)
        for method_id=1:Nmethod
            h=zeros(Nsim,1);
            p=zeros(Nsim,1);
            for sim_id=1:Nsim
                data00=dataset{sim_id};
                data0=data00{N_id,mu_id,sd_id}-20;
                switch method_set{method_id}
                    case 'ttest'
                        [h(sim_id),p(sim_id),~,~]=method{method_id}(data0,0);
                    case 'signtest'
                        [p(sim_id),h(sim_id),~]=method{method_id}(data0,0);
                    case 'mult_comp_perm_t1'
                        [p(sim_id),~,~]=method{method_id}(data0);
                        h(sim_id)=double(p(sim_id)<0.05);
                end
            end
            NHST.Prej(method_id,N_id,sd_id)=sum(h)/Nsim;
            NHST.pvalue(method_id,N_id,sd_id)=mean(p);
        end
    end
end
end

%% 
function norm_test=donormtest(dataset,x,y,Nsim,mu_id)
norm_test.Prej=zeros(x, y);
norm_test.pvalue=zeros(x, y);
for i=1:x
    for j=1:y
        h=zeros(Nsim,1);
        p=zeros(Nsim,1);
        for s=1:Nsim
            testdata0=dataset{s};
            testdata=testdata0{i,mu_id,j};
            [h(s),p(s),~]=swtest(testdata);
        end
        norm_test.Prej(i,j)=sum(h)/Nsim;
        norm_test.pvalue(i,j)=mean(p);
    end
end
end

%%
% -> paulbays.com
function kappa=sd2kappa(sd)
sd0=sd*pi/180;
R = exp(-sd0.^2/2);
kappa = nan(size(R));
ix = R >= 0 & R < 0.53;
kappa(ix) = 2 * R(ix) + R(ix).^3 + (5 * R(ix).^5)/6;
ix = R >= 0.53 & R < 0.85;
kappa(ix) = -0.4 + 1.39 * R(ix) + 0.43./(1 - R(ix));
ix = R >= 0.85 & R <= 1;
kappa(ix) = 1./(R(ix).^3 - 4 * R(ix).^2 + 3 * R(ix));
end
