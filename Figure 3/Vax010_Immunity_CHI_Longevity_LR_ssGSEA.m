%10/14/18
%Thomas Hagan
%Emory/CHI H5N1+AS03/TIV Ab Durability Analysis, linear regression prediction

clear
load('AffyU133plus2_Mar2016_entrez_symbol_desc.mat')
%Immunity (TIV Dataset)
load('Immunity_D180_D30_MaxFC_submatch.mat')
load('Immunity_FC_ssGSEA.mat')
probelist=BTM_labels;
Immunity_exp_FC=Immunity_ssGSEA_FC;
%Vax010 (Emory H5N1+AS03 Dataset)
load('Vax010_HAI_MN_SPR_FC.mat')
load('Vax010_FC_ssGSEA.mat')
%Split into adj/nonadj
load ('Vax010_Adj_Groups.mat')
adj_sub_labels_FC=cell(size(Vax010_ssGSEA_FC));
nonadj_sub_labels_FC=cell(size(Vax010_ssGSEA_FC));
adj_exp_FC=cell(size(Vax010_ssGSEA_FC));
nonadj_exp_FC=cell(size(Vax010_ssGSEA_FC));
for i=1:length(Vax010_ssGSEA_FC)
    [~,~,ind_adj]=intersect(Adj,Vax010_sub_labels_FC{i});
    adj_sub_labels_FC{i}=Vax010_sub_labels_FC{i}(ind_adj);
    adj_exp_FC{i}=Vax010_ssGSEA_FC{i}(:,ind_adj);
    [~,~,ind_nonadj]=intersect(Nonadj,Vax010_sub_labels_FC{i});
    nonadj_sub_labels_FC{i}=Vax010_sub_labels_FC{i}(ind_nonadj);
    nonadj_exp_FC{i}=Vax010_ssGSEA_FC{i}(:,ind_nonadj);
end
%CHI data (CHI H5N1+AS03 Dataset)
load('CHI_FC_ssGSEA.mat')
CHI_exp_FC=CHI_ssGSEA_FC;


%% Process Immunity Data

%Regression to define D180/D30 residual
p=polyfit(D30_D180_Max_FC(1,:),D30_D180_Max_FC(2,:),1);
Immunity_HAI_resid=D30_D180_Max_FC(2,:)-polyval(p,D30_D180_Max_FC(1,:));
Immunity_HAI_resid_sub_labels=D180_sub_labels;


%% Process Vax010 Data

%Convert titers/SPR FC to log2
Adj_MN=cellfun(@log2,Adj_MN,'UniformOutput',0);
Adj_MN_FC=cellfun(@log2,Adj_MN_FC,'UniformOutput',0);
Adj_HAI=cellfun(@log2,Adj_HAI,'UniformOutput',0);
Adj_HAI_FC=cellfun(@log2,Adj_HAI_FC,'UniformOutput',0);
Adj_SPR=cellfun(@log2,Adj_SPR,'UniformOutput',0);
Adj_SPR_FC=cellfun(@log2,Adj_SPR_FC,'UniformOutput',0);
Nonadj_MN=cellfun(@log2,Nonadj_MN,'UniformOutput',0);
Nonadj_MN_FC=cellfun(@log2,Nonadj_MN_FC,'UniformOutput',0);
Nonadj_HAI=cellfun(@log2,Nonadj_HAI,'UniformOutput',0);
Nonadj_HAI_FC=cellfun(@log2,Nonadj_HAI_FC,'UniformOutput',0);
Nonadj_SPR=cellfun(@log2,Nonadj_SPR,'UniformOutput',0);
Nonadj_SPR_FC=cellfun(@log2,Nonadj_SPR_FC,'UniformOutput',0);

%Calculate D42D21 FC for MN/HAI
[Adj_D42D21_MN_FC_sub,ia,ib]=intersect(Adj_MN_SPR_sub_labels{2},Adj_MN_SPR_sub_labels{3});
Adj_D42D21_MN_FC=Adj_MN{3}(ib)-Adj_MN{2}(ia);
[Adj_D42D21_HAI_FC_sub,ia,ib]=intersect(Adj_HAI_sub_labels{2},Adj_HAI_sub_labels{3});
Adj_D42D21_HAI_FC=Adj_HAI{3}(ib)-Adj_HAI{2}(ia);


%Calculate D100 residuals for MN and HAI titer
%MN
[Adj_D100_MN_resid_sub,ia,ib]=intersect(Adj_MN_SPR_FC_sub_labels{2},Adj_MN_SPR_FC_sub_labels{3});
p=polyfit(Adj_MN_FC{2}(ia),Adj_MN_FC{3}(ib),1);
Adj_D100_MN_resid=Adj_MN_FC{3}(ib)-polyval(p,Adj_MN_FC{2}(ia));
Adj_D100_MN_SSE=sum(Adj_D100_MN_resid.^2);
%Calculate D100D42 FC
Adj_D100D42_MN=Adj_MN_FC{3}(ib)-Adj_MN_FC{2}(ia);
%HAI
[Adj_D100_HAI_resid_sub,ia,ib]=intersect(Adj_HAI_FC_sub_labels{2},Adj_HAI_FC_sub_labels{3});
p=polyfit(Adj_HAI_FC{2}(ia),Adj_HAI_FC{3}(ib),1);
Adj_D100_HAI_resid=Adj_HAI_FC{3}(ib)-polyval(p,Adj_HAI_FC{2}(ia));
Adj_D100_HAI_SSE=sum(Adj_D100_HAI_resid.^2);
Adj_D100_HAI_resid_sub_labels=Adj_HAI_FC_sub_labels{3}(ib);
%Calculate D100D42 FC
Adj_D100D42_HAI=Adj_HAI_FC{3}(ib)-Adj_HAI_FC{2}(ia);

%% Correlate BTM/gene expression with residual

%Immunity
R_p_Immunity_HAI_resid=cell(size(Immunity_exp_FC));
R_p_Immunity_HAI_resid_table=cell(size(Immunity_exp_FC));
for i=1:size(Immunity_exp_FC,2)
    %Correlate with residual
    [~,ia,ib]=intersect(Immunity_sub_labels_FC{i},Immunity_HAI_resid_sub_labels);
    [R_p_Immunity_HAI_resid{i}(:,1),R_p_Immunity_HAI_resid{i}(:,2)]=...
        corr(Immunity_HAI_resid(ib)',Immunity_exp_FC{i}(:,ia)');
    %Sort by correlation
    [~,ind]=sort(R_p_Immunity_HAI_resid{i}(:,1),'descend');
    R_p_Immunity_HAI_resid_table{i}=[probelist(ind) num2cell(R_p_Immunity_HAI_resid{i}(ind,:))];
end

R_p_adj_HAI_resid=cell(size(adj_exp_FC));
R_p_adj_HAI_resid_table=cell(size(adj_exp_FC));
for i=1:size(adj_exp_FC,2)
    %Correlate with residual
    [~,ia,ib]=intersect(adj_sub_labels_FC{i},Adj_D100_HAI_resid_sub_labels);
    [R_p_adj_HAI_resid{i}(:,1),R_p_adj_HAI_resid{i}(:,2)]=...
        corr(Adj_D100_HAI_resid(ib)',adj_exp_FC{i}(:,ia)');
    %Sort by correlation
    [~,ind]=sort(R_p_adj_HAI_resid{i}(:,1),'descend');
    R_p_adj_HAI_resid_table{i}=[probelist(ind) num2cell(R_p_adj_HAI_resid{i}(ind,:))];
end

%Check intersection of top genes/probelist in Immunity/Vax010
p_cut=0.25;
common_BTM=cell(size(adj_exp_FC));
common_gene=cell(size(adj_exp_FC));
for i=1:size(Immunity_exp_FC,2)
    %probelist
    %Find significant probelist (immunity)
    ind_immunity_pos=find(R_p_Immunity_HAI_resid{i}(:,1)>0);
    ind_immunity_neg=find(R_p_Immunity_HAI_resid{i}(:,1)<0);
    ind_immunity_sig=find(R_p_Immunity_HAI_resid{i}(:,2)<p_cut);
    ind_immunity_pos_sig=intersect(ind_immunity_sig,ind_immunity_pos);
    ind_immunity_neg_sig=intersect(ind_immunity_sig,ind_immunity_neg);
    %Find significant probelist (vax010-prime)
    ind_adj_pos=find(R_p_adj_HAI_resid{i}(:,1)>0);
    ind_adj_neg=find(R_p_adj_HAI_resid{i}(:,1)<0);
    ind_adj_sig=find(R_p_adj_HAI_resid{i}(:,2)<p_cut);
    ind_adj_pos_sig=intersect(ind_adj_sig,ind_adj_pos);
    ind_adj_neg_sig=intersect(ind_adj_sig,ind_adj_neg);
    %Common
    ind_pos=intersect(ind_immunity_pos_sig,ind_adj_pos_sig);
    ind_neg=intersect(ind_immunity_neg_sig,ind_adj_neg_sig);
    ind=[ind_pos; ind_neg];
    common_BTM{i}=[probelist(ind) num2cell([R_p_Immunity_HAI_resid{i}(ind,:) R_p_adj_HAI_resid{i}(ind,:)])];
    %Find significant probelist (vax010-boost)
    ind_adj_pos=find(R_p_adj_HAI_resid{i+3}(:,1)>0);
    ind_adj_neg=find(R_p_adj_HAI_resid{i+3}(:,1)<0);
    ind_adj_sig=find(R_p_adj_HAI_resid{i+3}(:,2)<p_cut);
    ind_adj_pos_sig=intersect(ind_adj_sig,ind_adj_pos);
    ind_adj_neg_sig=intersect(ind_adj_sig,ind_adj_neg);
    %Common
    ind_pos=intersect(ind_immunity_pos_sig,ind_adj_pos_sig);
    ind_neg=intersect(ind_immunity_neg_sig,ind_adj_neg_sig);
    ind=[ind_pos; ind_neg];
    common_BTM{i+3}=[probelist(ind) num2cell([R_p_Immunity_HAI_resid{i}(ind,:) R_p_adj_HAI_resid{i+3}(ind,:)])];
end

%% Forward selection

%Select timepoint
day=3; %1=Day 1, 3=Day 7 %CHI has no Day 3, we are using boost data only
day_CHI=2; %1=Day 1, 2=Day 7 %CHI has no Day 3

probes=common_BTM{day+3}(:,1); %Boost

%Setup inputs/targets
%Vax010
%Setup targets
[~,ia,ib]=intersect(adj_sub_labels_FC{day+3},Adj_D100_HAI_resid_sub_labels); %Boost
targets_train=Adj_D100_HAI_resid(ib);
%Setup inputs
[~,ind]=ismember(probes,probelist);
probes_train=adj_exp_FC{day+3}(ind,ia); %Boost
%Immunity
%Setup targets
[~,ia,ib]=intersect(Immunity_sub_labels_FC{day},Immunity_HAI_resid_sub_labels);
targets_test=Immunity_HAI_resid(ib);
%Setup inputs
[~,ind]=ismember(probes,probelist);
probes_test=Immunity_exp_FC{day}(ind,ia);

%CHI (no targets as it is blind prediction)
%Setup inputs
[~,ind]=ismember(probes,probelist);
probes_CHI=CHI_exp_FC{day_CHI+2}(ind,:); %Boost

tic
%Total signatures to select
n_sig=1;
best_inputs_probes_total=cell(n_sig,1);
Means_v_inputs_total=cell(n_sig,1);
Means_final_total=cell(n_sig,1);
Std_final_total=cell(n_sig,1);
CHI_outputs_total=cell(n_sig,1);
R_p_total=cell(n_sig,1);
R_p_mean_std=cell(n_sig,1);
for q=1:n_sig
    remaining_inputs_train=probes_train;
    remaining_inputs_test=probes_test;
    remaining_inputs_CHI=probes_CHI;
    remaining_inputs_probes=probes;
    Means_v_inputs=[];
    accuracynew=1e9;
    bestaccuracy=1e9;
    best_inputs_train=[];
    best_inputs_test=[];
    best_inputs_CHI=[];
    best_inputs_probes={};
    %While there are remaining inputs to add, continue forward selection
    while isempty(remaining_inputs_train)==0

        Means=zeros(size(remaining_inputs_train,1),2);
        Means_sort=zeros(size(remaining_inputs_train,1),3);
        for k=1:size(remaining_inputs_train,1)
            %Build updated input variables
            inputs_train=[best_inputs_train; remaining_inputs_train(k,:)];
            inputs_test=[best_inputs_test; remaining_inputs_test(k,:)];

            %Set # of trials to repeat network training and testing
            nTrials = 1;
            
            percCorrect=zeros(nTrials,2);
            for i=1:nTrials
                
                %Compute linear regression
                LR=[ones(1,size(inputs_train,2)); inputs_train]'\targets_train';
                
                %Compute outputs
                outputs_test=([ones(1,size(inputs_test,2)); inputs_test]'*LR)';

                %Recalculate Training Output

                %Compute network outputs for training and validation points
                outputs_train=([ones(1,size(inputs_train,2)); inputs_train]'*LR)';
                
                %Compute performance
                percCorrect(i,1) = mean((outputs_train-targets_train).^2);
                percCorrect(i,2) = mean((outputs_test-targets_test).^2);

            end
            %Store train mean accuracy
            Means(k,1)=mean(percCorrect(:,1));
            %Store test mean accuracy
            Means(k,2)=mean(percCorrect(:,2));
            %Store Combined Variable for sorting
            Means_sort(k,:)=[Means(k,:) k];
        end
        %Sort by ascending mse in continuous mode
        %Means_sort=sortrows(Means_sort,1); %Sort on training
        Means_sort=sortrows(Means_sort,2); %Sort on test
        
        %Store best results for each input addition
        Means_v_inputs=[Means_v_inputs; Means_sort(1,1:end-1)];

        accuracyold=accuracynew;
        %accuracynew=Means_sort(1,1); %Sort on training
        accuracynew=Means_sort(1,2); %Sort on test
        %1-round check
        %If adding an input does not lower MSE, quit loop
        if accuracyold<=accuracynew
            break
        end
        
        %Set best inputs based on highest scoring combination
        best_inputs_train=[best_inputs_train; remaining_inputs_train(Means_sort(1,end),:)];
        best_inputs_test=[best_inputs_test; remaining_inputs_test(Means_sort(1,end),:)];
        best_inputs_CHI=[best_inputs_CHI; remaining_inputs_CHI(Means_sort(1,end),:)];
        %Store list of best inputs
        best_inputs_probes=[best_inputs_probes; remaining_inputs_probes(Means_sort(1,end),:)];
        %Set remaining inputs by removing input added to best inputs
        remaining_inputs_train(Means_sort(1,end),:)=[];
        remaining_inputs_test(Means_sort(1,end),:)=[];
        remaining_inputs_CHI(Means_sort(1,end),:)=[];
        remaining_inputs_probes(Means_sort(1,end),:)=[];
    end

    % Run full committee on best set of inputs

    inputs_train=best_inputs_train;
    inputs_test=best_inputs_test;
    inputs_CHI=best_inputs_CHI;

    nTrials = 10;
    %Perform cross-validation
    cv=cvpartition(size(inputs_train,2),'KFold',nTrials);

    percCorrect=zeros(nTrials,2);
    R_p=zeros(nTrials,4);
    outputs_CHI_all=zeros(nTrials,size(inputs_CHI,2));

    for i=1:nTrials
        %Compute linear regression
        LR=[ones(1,size(inputs_train(:,~cv.test(i)),2)); inputs_train(:,~cv.test(i))]'\targets_train(~cv.test(i))';

        %Compute outputs
        outputs_test=([ones(1,size(inputs_test,2)); inputs_test]'*LR)';
        outputs_CHI=([ones(1,size(inputs_CHI,2)); inputs_CHI]'*LR)';
        %Store CHI outputs for each trial
        outputs_CHI_all(i,:)=outputs_CHI;

        %Recalculate Training Output

        %Compute network outputs for training and validation points
        outputs_train=([ones(1,size(inputs_train,2)); inputs_train]'*LR)';

        %Compute performance
        %Use categorical outputs
        percCorrect(i,1) = 100*length(find((outputs_train(~cv.test(i))>0)==(targets_train(~cv.test(i))>0)))/length(targets_train(~cv.test(i)));
        percCorrect(i,2) = 100*length(find((outputs_test>0)==(targets_test>0)))/length(targets_test);
        %Compute correlations with residual
        [R_p(i,1),R_p(i,2)]=corr(outputs_train(~cv.test(i))',targets_train(~cv.test(i))');
        [R_p(i,3),R_p(i,4)]=corr(outputs_test',targets_test');
    end
    
    Means_final=zeros(1,size(percCorrect,2));
    Std_final=zeros(1,size(percCorrect,2));
    R_mean_final=zeros(1,2);
    R_std_final=zeros(1,2);
    for i=1:size(percCorrect,2)
        if i==1 %For training set use weighted mean/std
            Means_final(i)=sum(percCorrect(:,i).*cv.TrainSize')/sum(cv.TrainSize);
            Std_final(i)=sqrt(var(percCorrect(:,i),cv.TrainSize'));
            R_mean_final(i)=sum(R_p(:,i).*cv.TrainSize')/sum(cv.TrainSize);
            R_std_final(i)=sqrt(var(R_p(:,i),cv.TrainSize'));
            R_mean_final(i+1)=sum(R_p(:,i+2).*cv.TrainSize')/sum(cv.TrainSize);
            R_std_final(i+1)=sqrt(var(R_p(:,i+2),cv.TrainSize'));
        else
            Means_final(i)=mean(percCorrect(:,i));
            Std_final(i)=std(percCorrect(:,i));
        end
    end
    
    %Store results for each signature
    best_inputs_probes_total{q,1}=best_inputs_probes;
    Means_v_inputs_total{q,1}=Means_v_inputs;
    Means_final_total{q,1}=Means_final;
    Std_final_total{q,1}=Std_final;
    %Store output votes for CHI
    CHI_outputs_total{q,1}=outputs_CHI_all;
    %Store correlations (for continuous mode)
    R_p_total{q,1}=R_p;
    R_p_mean_std{q,1}=[R_mean_final R_std_final];
end
toc