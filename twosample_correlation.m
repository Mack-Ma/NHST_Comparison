%% Compare two-sample ttests (correlation)
% Tianye Ma, Nov 25
%

%% Prologue
clear
close all
BMWdir='/home/tianye/Matlab_WS/BMW/Bayesian_Modeling_of_Working_Memory';
addpath(BMWdir)
BMW('AddPath');
savedir='/home/tianye/Matlab_WS/NHST_MC';
if isempty(gcp('nocreate'))
    parpool('local',16); % start parallel computing (use 16 threads)
end
%% Simulation
mu_set=[0, 0.5, 1, 1.5, 2]; % mean difference
sd=2.5;
corr_set=0:0.02:0.98;
N_set=4:4:200;
Nsim=100;
%dataset=cell(Nsim,length(N_set),length(corr_set));
dataset1=cell(Nsim,1);
dataset2=cell(Nsim,1);
Nsim_left=Nsim;
fprintf('\nNow start simulation...\n')
sim_min_id=1;
% group 1
fprintf('\nGroup 1...\n')
dataset=cell(Nsim,1);
for i=1:10
    Nsim_1=min(Nsim_left,ceil(Nsim/10));
    for sim_id=sim_min_id:sim_min_id+Nsim_1-1
        dataset{sim_id}=datagenerate(mu_set,corr_set,N_set,sd);
    end
    Nsim_left=Nsim_left-Nsim_1;
    sim_min_id=sim_min_id+Nsim_1;
    fprintf('%d simulations left.\n',Nsim_left)
end
fprintf('\nSimulation finished...\n')

%% NHSTs
method_set={'ttest2','ttest'};
Nmethod=length(method_set);
method=cell(Nmethod,1);
for i=1:Nmethod
    method{i}=str2func(method_set{i});  
end
NHST=cell(length(mu_set),1);
fprintf('\nNow start doing NHSTs...\n')
for mu_id=1:length(mu_set)
    NHST{mu_id}=dotest(N_set,corr_set,method_set,method,dataset,Nsim,mu_id);
%     NHST1.Prej=zeros(length(method_set), length(N_set), length(corr_set));
% NHST1.pvalue=zeros(length(method_set), length(N_set), length(corr_set));
% for N_id=1:length(N_set)
%     N_id
%     for sd_id=1:length(corr_set)
%         sd_id
%         for method_id=1:length(method_set)
%             h=zeros(Nsim,1);
%             p=zeros(Nsim,1);
%             for sim_id=1:Nsim
%                 data00_1=dataset1{sim_id};
%                 data0_1=data00_1{N_id,mu_id,sd_id};
%                 data00_2=dataset2{sim_id};
%                 data0_2=data00_2{N_id,mu_id,sd_id};
%                 [h(sim_id),p(sim_id),~,~]=ttest2(data0_1,data0_2,'method',method_set{method_id});
%             end
%             NHST1.Prej(method_id,N_id,sd_id)=sum(h)/Nsim;
%             NHST1.pvalue(method_id,N_id,sd_id)=mean(p);
%         end
%     end
% end
% NHST{mu_id}=NHST1;
end
fprintf('\nNHSTs finished...\n')

%% Comparison
% set color scale
map_color=[[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(255,0,256-23)' ]]/255;
% heatmaps
figure(1)
% heatmaps of rejection rate
for j=1:length(method_set)
    for i=1:length(mu_set)
        subplot(length(method_set)+1,length(mu_set),(j-1)*length(mu_set)+i);
        testdata=NHST{i};
        imagesc(reshape(testdata.Prej(j,:,:),[length(N_set),length(corr_set)]))
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        caxis([-.1,1]);
        colormap(map_color);
    end
end
% heatmaps of rejection rate difference
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(corr_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(corr_set)]);
    subplot(length(method_set)+1,length(mu_set),i+length(method_set)*length(mu_set))
    imagesc(testdiff)
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    caxis([-.1,1]);
    colormap(map_color);
end
%
figure(2)
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(corr_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(corr_set)]);
    plot(mean(testdiff))
    hold on
end
figure(3)
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(corr_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(corr_set)]);
    plot(mean(testdiff,2))
    hold on
end


%% Epilogue
rmpath(BMWdir)
fprintf('\nSave all data...\n')
cd(savedir)
savename=sprintf('twosample_correlation');
save(savename)
fprintf('\nAll done!\n')

%%
function dataset=datagenerate(mu_set,corr_set,N_set,sd)
dataset=cell(length(N_set),length(mu_set),length(corr_set));
for mu_id=1:length(mu_set)
    for corr_id=1:length(corr_set)
        for N_id=1:length(N_set)
            mu0=[mu_set(mu_id),0];
            vc_m=[sd.^2,corr_set(corr_id)*sd.^2;corr_set(corr_id)*sd.^2,sd.^2];
            dataset0=mvnrnd(mu0,vc_m,N_set(N_id));
            dataset{N_id,mu_id,corr_id}=dataset0;
        end
    end
end
end

%%
function NHST=dotest(N_set,corr_set,method_set,method,dataset,Nsim,mu_id)
NHST.Prej=zeros(length(method_set), length(N_set), length(corr_set));
NHST.pvalue=zeros(length(method_set), length(N_set), length(corr_set));
for N_id=1:length(N_set)
    for corr_id=1:length(corr_set)
        for method_id=1:length(method_set)
            h=zeros(Nsim,1);
            p=zeros(Nsim,1);
            for sim_id=1:Nsim
                data00_1=dataset{sim_id};
                data0_1=data00_1{N_id,mu_id,corr_id};
                switch method_set{method_id}
                    case 'ttest'
                        [h(sim_id),p(sim_id),~,~]=method{method_id}(data0_1(:,1),data0_1(:,2));
                    case 'ttest2'
                        [h(sim_id),p(sim_id),~,~]=method{method_id}(data0_1(:,1),data0_1(:,2));
                end
            end
            NHST.Prej(method_id,N_id,corr_id)=sum(h)/Nsim;
            NHST.pvalue(method_id,N_id,corr_id)=mean(p);
        end
    end
end
end
