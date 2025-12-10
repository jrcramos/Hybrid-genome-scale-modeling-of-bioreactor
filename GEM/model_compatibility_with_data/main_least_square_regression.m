%---------------------------------------------------------------------------
% This codes loads the exchange fluxes in phases P1*, P2*, P3*, P4*
% and P4* calculated from the DoE data and
% performs linear least square regression of the fluxes predicted by the
% GEM. This ensures that the GEM is capable of being used as the last layer
% in the neural network to predict exchange fluxes. An average flux from
% the 9 DoE cultivations is used for this purpose (data(10).val), but this
% test can be performed for an individual fed-batch by changing "ib" in
% line 34
% *(quasi)steady-state phases: P1-Exponential cell growth (approximately 0-36 h),
% P2-Second Exponential cell growth (approximately 38-84 h),
% P3-Early stationary (approximately 100-144 h), 
% P4-Late stationary (approximately 146-240 h)
%---------------------------------------------------------------------------

close all
load('data.mat')
load('model.mat')

batchs=1:9;
ratesp1=[];ratesp2=[];ratesp3=[];ratesp4=[];
for i=batchs % get experimental exchnage fluxes
    ratesp1=[ratesp1;data(i).val(:,1)'];
    ratesp2=[ratesp2;data(i).val(:,2)'];
    ratesp3=[ratesp3;data(i).val(:,3)'];
    ratesp4=[ratesp4;data(i).val(:,4)'];
end
data(10).val(:,1)=mean(ratesp1);data(10).val(:,2)=mean(ratesp2);data(10).val(:,3)=mean(ratesp3);data(10).val(:,4)=mean(ratesp4);
data(10).std(:,1)=std(ratesp1);data(10).std(:,2)=std(ratesp2);
data(10).std(:,3)=std(ratesp3);data(10).std(:,4)=std(ratesp4);data(10).batchid={'mean'}; % an average batch exchange rates


for ib=10

    for phase=1:4 %(quasi)steady-state phases
        pred=[];exp_data=[];exp_data_sd=[];
        [feasCheckStatus,vopt,vext,d,sd,fval]=LsqLin_data(data,model,phase,ib);
        if feasCheckStatus>=0
            pred=[pred vext];
            exp_data=[exp_data d'];
            exp_data_sd=[exp_data_sd sd'];
        end

        figure('units','normalized','outerposition',[0 0.1 0.35 0.35])

        errorbar(exp_data,pred,exp_data_sd,'horizontal','MarkerSize',13,'LineWidth',1.5,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
        mod=fitlm(reshape(exp_data,numel(exp_data),1),reshape(pred,numel(pred),1),'linear','RobustOpts','off');
        R2=mod.Rsquared;
        mod.plot('Marker','none') 
        title(['R^{2} ' num2str(R2.Adjusted,'%.3f') ' : P' num2str(phase)],'FontName','Times','fontsize',18,'FontWeight','Bold','interpreter','tex')
        ylabel('Prediction')
        xlabel('Experimental')
        y1=min(min(exp_data-exp_data_sd))-0.1;y2=max(max(exp_data-exp_data_sd))+0.2;
        ylim([y1 y2]);xlim([y1 y2]);
        axis('square')
        line([y1 y2],[y1 y2],'color','k')
        legend('off')
        % set(gca,'FontWeight','Bold','FontSize',12);
        % saveas(gcf,['lsqlin_mean' '_P'  num2str(phase)  '.tiff'])
    end
end


