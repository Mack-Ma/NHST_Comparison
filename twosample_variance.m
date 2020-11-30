%% Compare two-sample ttests (variance)
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
    parpool('local',12); % start parallel computing (use 16 threads)
end
%% Simulation
mu_set=[0, 0.5, 1, 1.5, 2]; % mean difference
sd_set=-2.4:0.1:2.5; % sd difference
sd_control=2.5;
N_set=4:4:200;
Nsim=100;
%dataset=cell(Nsim,length(N_set),length(sd_set));
dataset1=cell(Nsim,1);
dataset2=cell(Nsim,1);
Nsim_left=Nsim;
fprintf('\nNow start simulation...\n')
sim_min_id=1;
% group 1
fprintf('\nGroup 1...\n')
for i=1:10
    Nsim_1=min(Nsim_left,ceil(Nsim/10));
    for sim_id=sim_min_id:sim_min_id+Nsim_1-1
        dataset1{sim_id}=datagenerate(zeros(size(mu_set)),sd_control*ones(size(sd_set)),N_set);
    end
    Nsim_left=Nsim_left-Nsim_1;
    sim_min_id=sim_min_id+Nsim_1;
    fprintf('%d simulations left.\n',Nsim_left)
end
fprintf('\nGroup 1 finished...\n')
fprintf('\nGroup 2...\n')
% group 2
Nsim_left=Nsim;
sim_min_id=1;
for i=1:10
    Nsim_1=min(Nsim_left,ceil(Nsim/10));
    for sim_id=sim_min_id:sim_min_id+Nsim_1-1
        dataset2{sim_id}=datagenerate(mu_set,sd_set+sd_control,N_set);
    end
    Nsim_left=Nsim_left-Nsim_1;
    sim_min_id=sim_min_id+Nsim_1;
    fprintf('%d simulations left.\n',Nsim_left)
end
fprintf('\nGroup 2 finished...\n')

%% NHSTs
vartype_set={'equal','unequal'};
% NHST=cell(length(mu_set),1);
fprintf('\nNow start doing NHSTs...\n')
% for mu_id=1:length(mu_set)
%     NHST{mu_id}=dotest(N_set,sd_set,vartype_set,dataset1,dataset2,Nsim,mu_id);
%     NHST1.Prej=zeros(length(vartype_set), length(N_set), length(sd_set));
% NHST1.pvalue=zeros(length(vartype_set), length(N_set), length(sd_set));
% for N_id=1:length(N_set)
%     N_id
%     for sd_id=1:length(sd_set)
%         sd_id
%         for vartype_id=1:length(vartype_set)
%             h=zeros(Nsim,1);
%             p=zeros(Nsim,1);
%             for sim_id=1:Nsim
%                 data00_1=dataset1{sim_id};
%                 data0_1=data00_1{N_id,mu_id,sd_id};
%                 data00_2=dataset2{sim_id};
%                 data0_2=data00_2{N_id,mu_id,sd_id};
%                 [h(sim_id),p(sim_id),~,~]=ttest2(data0_1,data0_2,'vartype',vartype_set{vartype_id});
%             end
%             NHST1.Prej(vartype_id,N_id,sd_id)=sum(h)/Nsim;
%             NHST1.pvalue(vartype_id,N_id,sd_id)=mean(p);
%         end
%     end
% end
% NHST{mu_id}=NHST1;
% end
fprintf('\nNHSTs finished...\n')
% variance test
NHST_var=cell(length(mu_set),1);
fprintf('\nNow start doing variance tests...\n')
parfor mu_id=1:length(mu_set)
    NHST_var{mu_id}=dovariancetest(dataset1,dataset2,length(N_set),length(sd_set),Nsim,mu_id);
end
fprintf('\nVariance tests finished...\n')

%% Comparison
% set color scale
map_color=[[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(0,255,256-23)' ],[linspace(0,255,23)';linspace(255,0,256-23)' ]]/255;
% heatmaps of rejection rate
for j=1:length(vartype_set)
    for i=1:length(mu_set)
        subplot(length(vartype_set)+1,length(mu_set),(j-1)*length(mu_set)+i);
        testdata=NHST{i};
        imagesc(reshape(testdata.Prej(j,:,:),[length(N_set),length(sd_set)]))
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        caxis([-.1,1]);
        colormap(map_color);
    end
end
% heatmaps of rejection rate difference
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(sd_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(sd_set)]);
    subplot(length(vartype_set)+1,length(mu_set),i+length(vartype_set)*length(mu_set))
    imagesc(testdiff)
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    caxis([-.1,1]);
    colormap(map_color);
end
% variance (mean)
figure(4)
var_sum=zeros(length(N_set),length(sd_set));
for i=1:length(mu_set)
    var_sum=var_sum+NHST_var{i}.Prej(:,:);
end
imagesc(var_sum/length(mu_set))
threshold_var=.5;
mask=var_sum/length(mu_set)>threshold_var;
hold on
plot((1-mean(mask))*length(N_set)+.5,'LineWidth',2,'color','r')
% 
figure(5)
diff_mask=zeros(length(mu_set),1);
for i=1:length(mu_set)
    testdata=NHST{i};
    testdiff=reshape(testdata.Prej(2,:,:),[length(N_set),length(sd_set)])-reshape(testdata.Prej(1,:,:),[length(N_set),length(sd_set)]);
    diff_mask(i)=sum(sum(testdiff.*double(mask)))/sum(sum(mask));
    if mu_set(i)==0
        diff_mask(i)=-diff_mask(i);
    end
end
bar(diff_mask);

%% Epilogue
rmpath(BMWdir)
fprintf('\nSave all data...\n')
cd(savedir)
savename=sprintf('twosample_variance');
save(savename)
fprintf('\nAll done!\n')

%%
function dataset=datagenerate(mu_set,sd_set,N_set)
dataset=cell(length(N_set),length(mu_set),length(sd_set));
for mu_id=1:length(mu_set)
    for sd_id=1:length(sd_set)
        for N_id=1:length(N_set)
            param_n=[mu_set(mu_id),sd_set(sd_id)];
            dataset0=param_n(1)+param_n(2)*randn([N_set(N_id),1]);
            dataset{N_id,mu_id,sd_id}=dataset0;
        end
    end
end
end

%%
function NHST=dotest(N_set,sd_set,vartype_set,dataset1,dataset2,Nsim,mu_id)
NHST.Prej=zeros(length(vartype_set), length(N_set), length(sd_set));
NHST.pvalue=zeros(length(vartype_set), length(N_set), length(sd_set));
for N_id=1:length(N_set)
    for sd_id=1:length(sd_set)
        for vartype_id=1:length(vartype_set)
            h=zeros(Nsim,1);
            p=zeros(Nsim,1);
            for sim_id=1:Nsim
                data00_1=dataset1{sim_id};
                data00_2=dataset2{sim_id};
                [h(sim_id),p(sim_id),~,~]=ttest2(data00_1{N_id,mu_id,sd_id},data00_2{N_id,mu_id,sd_id},'vartype',vartype_set{vartype_id});
            end
            NHST.Prej(vartype_id,N_id,sd_id)=sum(h)/Nsim;
            NHST.pvalue(vartype_id,N_id,sd_id)=mean(p);
        end
    end
end
end

%%
function var_test=dovariancetest(dataset1,dataset2,x,y,Nsim,mu_id)
var_test.Prej=zeros(x, y);
var_test.pvalue=zeros(x, y);
for i=1:x
    for j=1:y
        h=zeros(Nsim,1);
        p=zeros(Nsim,1);
        for s=1:Nsim
            data00_1=dataset1{s};
            data00_2=dataset2{s};
            [p(s),~]=vartestn([data00_1{i,mu_id,j},data00_2{i,mu_id,j}],'testtype','LeveneQuadratic','displayopt','off');
            h(s)=double(p(s)<.05);
        end
        var_test.Prej(i,j)=sum(h)/Nsim;
        var_test.pvalue(i,j)=mean(p);
    end
end
end