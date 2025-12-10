function [feasCheckStatus,vopt,vext,d,sd,fval]=LsqLin_data(data,model,phase,batch)

Aint=model.Sint;
Aext=model.Sext;
metstr_ext = model.mets(~model.isintr);
mets_data=data(1).names;
ids=[];
for i=1:numel(mets_data)
    ids=[ids strmatch(mets_data{i},metstr_ext)];
end
Aext=Aext(ids,:);

ub=ones(size(Aint,2),1)*100;
lb=zeros(size(Aint,2),1);lb(model.isrev==1)=-100;

exp_data=[data(batch).val(:,phase)']; % flux
exp_data_sd=[data(batch).std(:,phase)']; % flux std
% exp_data_sd=[abs(data(batch).val(:,phase))'*0.3]; % flux std
idbio=strmatch('Xv[e]',data(1).names);
exp_data(idbio)=max([exp_data(idbio) 0]); % if mu<0 then their flux is 0
ids=find(~isnan(exp_data)); % find for which metabolite exp flux data is available




ub_ex=exp_data(ids)+exp_data_sd(ids);
lb_ex=exp_data(ids)-exp_data_sd(ids);

A=[Aext(ids,:);-Aext(ids,:)];
b=[ub_ex';-lb_ex'];


options = optimoptions('lsqlin','MaxIter',30500,'Display','none','TolCon',1e-6,'TolX',1e-6,'OptimalityTolerance',5e-6);%);%,'Algorithm','active-set','OptimalityTolerance',1e-3,'ConstraintTolerance', 1e-7);%'OptimalityTolerance',1e-12);%,'TolCon',1e-9) ;
options = optimoptions('lsqlin','Algorithm','active-set','ConstraintTolerance',1e-9,...
        'MaxIter',20000,'Display','iter');
options = optimoptions('lsqlin','Algorithm','active-set','ConstraintTolerance',1e-7,...
        'MaxIter',20000,'Display','iter');%,'ObjectiveLimit',0, 'OptimalityTolerance',1e-20,'StepTolerance',1e-20 );
sd=exp_data_sd(ids);
exp_data_sd(exp_data_sd==0)=1;
C=Aext(ids,:)./repmat(exp_data_sd(ids)',1,size(Aext(ids,:),2));
d=exp_data(ids)./exp_data_sd(ids);d(isnan(d))=1;
% [vref,SSE] = lsqlin(C,d,A,b,Aeq,beq,lb,ub,zeros(size(lb)),options);
[vopt,~,~,feasCheckStatus] = lsqlin(C,d,A,b,Aint,zeros(size(Aint,1),1),lb,ub,zeros(size(lb)),options);

d=exp_data(ids);
if feasCheckStatus>=0
    vext=Aext(ids,:)*vopt;
    exp_data_sd(exp_data_sd==0)=1e-6;
    fval=(vext'-exp_data(ids))./exp_data_sd(ids);
    fval=fval'*fval;
else
    fprintf('Problem is unfeasible for "%s" cultivation phase %d \n',data(batch).batchid{:},phase);
    vopt=[];fval=inf;vext=NaN(numel(ub_ex),1);
end
end

