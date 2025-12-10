function [control_to_subsys]=plot_control_nodes_paper(simfile,batchs)

hek_subsystems={1:12	'Glycolysis'
                13:20	'PPP'
                21:31	'TCA'
                32:32	'Urea cycle'
                33:36	'Glutaminolysis'
                37:38	'Glycogen synthesis'
                39:127	'AA metabolism'
                128:128	'Biomass synthesis'
                129:135	'Folate/C1 metabolism'
                136:157	'Nucleotide interconversion'
                158:177	'Nucleotide synthesis'
                178:179	'Inositol metabolism'
                180:186	'Sphingolipids'
                187:204	'Glico-phospholipid metabolism'
                205:212	'Fatty acids metabolism'
                213:219	'Carotenoids and isoprenoids'
                220:230	'Steroid metabolism/cholesterol synthesis'
                231:239	'Reticulum to cytoplasm'
                240:243	'Nucleus to cytoplasm'
                244:285	'Mitochondria to cytoplasm'
                286:292	'Oxidative phosphorylation/ maintenance'
                293:329	'Exchange Reactions'};


load([simfile '.mat']);
data=runs.data_iter{1, ind};
% get control nodes
[~,data]=hnetsimul(hnet,data,0,1);
train=[];test=[];
for ib=batchs
    if data.batch(ib).istrain==1
        train=[train ib];
    else
        test=[test ib];
    end
end

% for phase=1:2
%     if phase==1
%         idphase=data.age_h<=120;
%     else
%         idphase=data.age_h>=170;
%     end
    idphase=1:numel(data.age_h);
    control_to_subsys=zeros(size(data.batch(2).control_nodes,2),size(hek_subsystems,1));
    for inodes=1:size(data.batch(2).control_nodes,2)
        for isubsys=1:size(hek_subsystems,1)
            temp=0;temp2=0;
            for id=hek_subsystems{isubsys,1}
                for ib=batchs
                    if data.batch(ib).istrain==1
                        y=data.batch(ib).rates_full(idphase,id);
                        x=data.batch(ib).control_nodes(idphase,inodes);
                        mod=fitlm(x,y,'linear','RobustOpts','off');
                        R_sq=mod.Rsquared.Ordinary;temp=temp+R_sq;temp2=temp2+1;
                    end
                end
            end
            control_to_subsys(inodes,isubsys)=temp/temp2;
        end
    end
%     if phase==1
%         control_to_subsys_growth=control_to_subsys;
%     else
%         control_to_subsys_death=control_to_subsys;
%     end
% end


