function plot_predicted_intra_rates_paper(simfile,nrow,batchs)

rev=[false,true,false,true,true,true,true,true,true,false,true,false,true,false,false,...
    true,true,true,true,true,false,false,true,false,false,true,true,true,true,true,false,...
    false,true,true,false,false,true,true,true,false,false,false,true,true,false,false,true,...
    false,true,true,true,true,false,false,false,true,true,false,true,true,true,false,true,true,...
    false,false,false,true,true,true,true,false,true,false,false,true,true,true,false,true,false,true,...
    false,false,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,...
    false,true,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,...
    false,false,false,false,false,false,false,false,false,false,true,true,true,true,true,false,false,true,true,...
    true,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,...
    true,false,true,false,true,true,true,true,false,true,true,true,false,true,false,true,false,false,false,false,false,...
    false,false,false,false,true,false,false,false,false,false,false,false,false,true,false,false,true,false,true,false,...
    false,false,false,false,false,false,false,true,true,false,false,false,false,false,true,false,false,false,true,false,...
    false,false,true,true,true,true,true,true,false,false,false,false,true,true,true,true,true,true,false,false,true,true,...
    true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,...
    true,true,true,true,true,true,true,false,true,true,true,false,false,true,false,false,false,false,false,false,false,false,...
    false,true,true,false,false,false,false,true,true,false,false,true,false,true,false,true,true,false,true,true,true,true,...
    true,true,true,false,true,true,true,true,true,true,true,true,true,true,false,true,false,false,true,true, true];
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
                293:329	'Exchange Reactions'
                find(rev==true) 'Reversible'
                find(rev==false) 'Non-reversible'};


load([simfile '.mat']);
ci=1;
names=hek_subsystems(:,2);names2=names;

divi=1;divie=numel(names);

Names=names;Names2=names2;
Names=names;Names2=names2;
for i=2:numel(Names)
    Names{i}=['  ' Names{i}];
end
lis='A':'Z';
fwidth=20.5*0.85;
fheight=26*0.85;
xlabelheight=0.7;
ylabelwidth=2;
xspaceheight=0.5;
yspacewidth=0.5;
pheight=(fheight-xlabelheight-5.5*xspaceheight)/9;
pwidth=(fwidth-ylabelwidth-yspacewidth*3.2)/5.3;
% plot_rates_pos
Pos1=[  16-(16-1)/4*4,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*1.0275,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];

Pos2=[  16-(16-1)/4*3,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*1.0275,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];

Pos3=[  16-(16-1)/4*2,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*1.0275,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];
    
Pos4=[  16-(16-1)/4*1,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*1.0275,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];
  
  Pos5=[16-(16-1)/4*0,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*1.0275,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];


Pos10=Pos1;Pos20=Pos2;Pos30=Pos3;Pos40=Pos4;Pos50=Pos5;
for ic=1:ci
    Pos1=Pos10;Pos2=Pos20;Pos3=Pos30;Pos4=Pos40;Pos5=Pos50;
    names=Names(divi(ic):divie(ic));names2=Names2(divi(ic):divie(ic));
    max_c=7-ceil(numel(names)/nrow)+2;
    Pos1=Pos1(max_c:end,:);Pos2=Pos2(max_c:end,:);Pos3=Pos3(max_c:end,:);Pos4=Pos4(max_c:end,:);Pos5=Pos5(max_c:end,:);
    hFig = figure;
    set(hFig, 'Color', 'w',...
        'paperpositionmode', 'auto',...
        'paperunits', 'centimeters',...
        'paperposition', [0.5 1. 31*0.86 40*0.85],...
        'units', 'centimeters',...
        'position', [0.5 1.5 17.5*0.86 23.7*0.85]);
    %[0,0.64,1]=[0,0.8,0.8];
    
    fwidth=20.5*0.85;
    fheight=26*0.85;
    xlabelheight=0.7;
    ylabelwidth=2;
    xspaceheight=0.5;
    yspacewidth=0.5;
    pheight=(fheight-xlabelheight-5.5*xspaceheight)/9;
    pwidth=(fwidth-ylabelwidth-yspacewidth*3.2)/5.3;
    
    
    id=0;
    data=runs.data_iter{1, ind};
    np=numel(names);
    train=[];test=[];
    for ib=batchs
        if data.batch(ib).istrain==1
            train=[train ib];
        else
            test=[test ib];
        end
    end
            
        
    end
    try
        for i=1:7 % rows
            maxys=[];maxxs=[];
            for row=1:nrow %coolumns
                if row==1; position=Pos1(i,:);elseif row==2;position=Pos2(i,:);elseif row==3;position=Pos3(i,:);elseif row==4;position=Pos4(i,:);else; position=Pos5(i,:);end
                id=id+1;if id>np;error;end
                hAx1t=axes('Parent', hFig,'XColor',[1 1 1],'YColor',[1 1 1 ]);
                
                
                set(hAx1t,'Units','centimeters','position',position,'YTickLabel',[],'XTickLabel',[]);
                hAx1=axes('Parent', hFig,'Color','none');
                hold(hAx1,'on');hold on;
                infm=ones(size(data.batch(ib).rates_full,1),1);
                maxy1=-inf*infm;miny1=infm*inf; maxx1=-inf*infm;minx1=inf*infm;
                for ib=batchs
                    x = data.age_h';
                    y = data.batch(ib).rates_full(:,hek_subsystems{id,1});
                    maxy1=max(max([maxy1.*infm y]));miny1=min(min([miny1.*infm y]));
                    maxx1=max(max([maxx1.*infm x]));minx1=min(min([minx1.*infm x]));
                    if data.batch(ib).istrain==1
                        trainingPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#0080ff','Marker','none','LineStyle','-','Linewidth',1.2,'HandleVisibility','off');
                        
                    else
                    end
                end
                maxys(row)=maxy1;maxxs(row)=maxx1;
                for ib=batchs
                    x = data.age_h';
                    y = data.batch(ib).rates_full(:,hek_subsystems{id,1});
                    if data.batch(ib).istrain==1

                    else
                        testPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#2bbf5c','Marker','none','LineStyle','-','Linewidth',0.25,'HandleVisibility','off');
                        
                    end
                end
                
%                 if id==np
%                     il=[];
%                     if ~isempty(train)
%                         il=[il 1];
%                         x = data.age_h';
%                         y = data.batch(train(1)).rates_full(:,hek_subsystems{id,1});
%                         trainingPlot = plot(x,y(:,1),...
%                             'MarkerSize',5,'Color','#0080ff','Marker','none','LineStyle','-','Linewidth',1.2,'HandleVisibility','on');
%                     end
%                     if ~isempty(test)
%                         il=[il 2];
%                         x = data.age_h';
%                         y = data.batch(test(1)).rates_full(:,hek_subsystems{id,1});
%                         testPlot = plot(x,y(:,1),...
%                             'MarkerSize',5,'Color','#2bbf5c','Marker','none','LineStyle','-','Linewidth',1.2,'HandleVisibility','on');
%                     end
%                     lg={'train' 'test'};
%                     legend(lg(il),'Location', 'north','FontSize',9,'FontName','Arial')
%                 end

                
                
                %             pl=line([48.1667 48.1667],[min([0  maxy1]) max([  maxy1] )]*1.3,'LineStyle','-','LineWidth',1.1);
                %             pl=line([48.1667+12 48.1667+12],[min([0  maxy1]) max([ maxy1] )]*1.3,'LineStyle','-','LineWidth',1.1);
                
                miny1=miny1-abs(mean([miny1 maxy1]))*1.45;maxy1=maxy1+abs(mean([miny1 maxy1]))*1.45;

                ylim(hAx1,[ min([miny1 ]) max([maxy1])*1.45])
                xlim(hAx1,[ 0 max([maxx1])*1.15])
                if row==1
                    if max([abs(miny1) maxy1])>0.015 && max([abs(miny1) maxy1])<0.1
                        hAx1.YAxis.Exponent=-2;
                    elseif max([abs(miny1) maxy1])>0.1
                        hAx1.YAxis.Exponent=-1;
                    end
                end
                axis square
%                 if np-id>nrow-1
                    set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
%                 else
%                     set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
%                     
%                 end
                
                hAx = gca;             % handle to current axes
                
                
                if ~strcmp(names(id),'  Biomass synthesis')
                    ylabel({'(mmol/gDWh)'},'FontSize',7,'interpreter','tex','FontName','Arial','FontWeight','bold');
                    if id>23
                        ylabel({'(mmol/gDWh or h^{-1})'},'FontSize',7,'interpreter','tex','FontName','Arial','FontWeight','bold');
                    end
                else
                    ylabel({'(h^{-1})'}','FontSize',7,'interpreter','tex','FontName','Arial','FontWeight','bold');
                end
                if id>22
                    xlabel('Time (h)','FontSize',7,'FontName','Arial','FontWeight','bold');
                end

                
                n=names2{id};
                
                set(gca,'FontSize',8);
                yrange=max([maxy1])-min([miny1 ]);
                xrange=max([ maxx1])-min([ minx1]);
%                 if id<3 || id==20
%                     text(min([ minx1])+xrange*0.75,min([miny1 ])+yrange*1.13,lis(id),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','bold')
%                 else
                    text(min([ minx1])+xrange*1,min([miny1 ])+yrange*1.225,lis(id),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',11.5,'FontName','Arial','FontWeight','bold')
%                 end
                
            end

            if ic>1 && i==5
                6;
            end
            ids=(i-1)*row+1:(i-1)*row+row;

                       
        end
    catch kl
    end


end
