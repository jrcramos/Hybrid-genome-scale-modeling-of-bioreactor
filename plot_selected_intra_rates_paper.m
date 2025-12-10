function plot_selected_intra_rates_paper(simfile,batchs,varargin)
if nargin==2
    ids=[1 13 10 40 21 11 42 50 24 26 35 28 291 286 297 305 304];
elseif nargin>2
    ids=varargin{1};
end

load([simfile '.mat']);

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
Pos1=[16-(16-1)/4*4.3,24-(24-1.24)/7*7.425,2.60849056603774,2.07222222222222];
id=0;
data=runs.data_iter{1, ind};
train=[];test=[];
for ib=batchs
    if data.batch(ib).istrain==1
        train=[train ib];
    else
        test=[test ib];
    end
end
ysim=[];
for ib=batchs
    data.batch(ib).rates_full(:,42)=data.batch(ib).rates_full(:,42)+data.batch(ib).rates_full(:,41); % r_ME1 + r_ME2
    data.batch(ib).rates_full(:,50)=data.batch(ib).rates_full(:,50)+data.batch(ib).rates_full(:,49); % r_AspTa1 + r_AspTa2
    ysim=[ysim;data.batch(ib).rates_full(:,ids)];
end

Names={};
for i=1:329;Names{i,1}=['R' num2str(i)];end
Names_model={'R1' 'R13' 'R10' 'R40' 'R21' 'R11' 'R42' 'R50' 'R24' 'R26' 'R35' 'R28' 'R291' 'R286' 'R297' 'R305' 'R304'};
Names={'r_{HK}' 'r_{G6PDH}' 'r_{PK}' 'r_{PC}' 'r_{PDH}' 'r_{LDH}' 'r_{ME}' 'r_{AspTA}', ...
                      'r_{ICDH}' 'r_{KDH}' 'r_{GS}' 'r_{SDH}' 'r_{ATPase}' 'r_{mATP}' 'r_{GLCT}' 'r_{GLNT}' 'r_{GLUT}'};

for ic=1:numel(ids)
    hFig = figure;
    set(hFig, 'Color', 'w',...
        'paperpositionmode', 'auto',...
        'paperunits', 'centimeters',...
        'paperposition', [0.5 1. 31*0.86 40*0.85],...
        'units', 'centimeters',...
        'position', [0.5 1.5 17.5*0.86*0.287 23.7*0.85*0.086]);
    %[0,0.64,1]=[0,0.8,0.8];
    
    fwidth=20.5*0.85;
    fheight=26*0.85;
    xlabelheight=0.7;
    ylabelwidth=2;
    xspaceheight=0.5;
    yspacewidth=0.5;
    pheight=(fheight-xlabelheight-5.5*xspaceheight)/9;
    pwidth=(fwidth-ylabelwidth-yspacewidth*3.2)/5.3;
    
    
    
    
    
    
    maxys=[];maxxs=[];
    
    position=Pos1(end,:);
    id=id+1;
    hAx1t=axes('Parent', hFig,'XColor',[1 1 1],'YColor',[1 1 1 ]);
    
    
    set(hAx1t,'Units','centimeters','position',position,'YTickLabel',[],'XTickLabel',[]);
    hAx1=axes('Parent', hFig,'Color','none');
    hold(hAx1,'on');hold on;
    infm=ones(size(data.batch(ib).rates_full,1),1);
    x = data.age_h';
    maxy1=-inf*infm;miny1=infm*inf; maxx1=-inf*infm;minx1=inf*infm;
    y=[];y1=[];y2=[];y3=[];id1=x<=75;id2=x<=120;id2(id1)=false;id3=x<=520;id3(id1)=false;id3(id2)=false;
    for ib=batchs
        if data.batch(ib).istrain==1
            y = [y data.batch(ib).rates_full(:,ids(id))];
            y1=[y1; data.batch(ib).rates_full(id1,ids(id))];
            y2=[y2; data.batch(ib).rates_full(id2,ids(id))];
            y3=[y3; data.batch(ib).rates_full(id3,ids(id))];
            maxy1=max(max([maxy1.*infm y]));miny1=min(min([miny1.*infm y]));
            maxx1=max(max([maxx1.*infm x]));minx1=min(min([minx1.*infm x]));
        else
        end
    end
    d=[mean(y1) mean(y2) mean(y3)];
    barh(1:3,ones(1,3),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Linewidth',1.2,'HandleVisibility','off');
    c=[0 0.5 0; 
        0 0.2 0.8;
        1 0.4 0; ];c=c([ 2 1 3],:);ibb=1;
    for ibar=[3 2 1]
        barh(ibb,d(ibar)/max(max(ysim)),'FaceColor',c(ibar,:),'EdgeColor',c(ibar,:),'HandleVisibility','off');ibb=ibb+1;
    end
    ax = gca;
    ax.XColor = [1 1 1];
    ax.YColor = [1 1 1];

    
%     miny1=miny1-abs(mean([miny1 maxy1]))*1.45;maxy1=maxy1+abs(mean([miny1 maxy1]))*1.45;
%     miny1=1;
%     
%     ylim(hAx1,[ min([miny1 ]) max([maxy1])*1.45])
%     xlim(hAx1,[ 0 max([maxx1])*1.15])
    
    axis square
    %                 if np-id>nrow-1
    set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
    %                 else
    %                     set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
    %
    %                 end
    ibb=3;
    for itxt=1:3
        n=d(ibb);ibb=ibb-1;
        if abs(n) > 0.1
            strn = sprintf('%.2f', n); % Convert to string with one decimal place
        elseif abs(n) < 0.1
            strn = sprintf('%1.1e', n); % Convert to scientific notation with one decimal place
        end
        text(1+0.1,itxt+0.45,strn,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',14,'FontName','Arial','FontWeight','bold')
    end
    hAx = gca;             % handle to current axes
    
    
    
%     ylabel({[Names{id} ' (mmol/gDWh)']},'FontSize',7,'interpreter','tex','FontName','Arial','FontWeight','bold');
        
    
    
    set(gca,'FontSize',8);
%     yrange=max([maxy1])-min([miny1 ]);
%     xrange=max([ maxx1])-min([ minx1]);
    %                 if id<3 || id==20
    %                     text(min([ minx1])+xrange*0.75,min([miny1 ])+yrange*1.13,lis(id),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','bold')
    %                 else
%     text(min([ minx1])+xrange*1,min([miny1 ])+yrange*1.225,lis(id),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',11.5,'FontName','Arial','FontWeight','bold')
    %                 end
    
    
   export_fig(['C:\Users\joao3\OneDrive\Documents\SilicoMedia\Publications\HybridGEM\Figs\selected_rates\' Names{id} '.png'],'-q101','-transparent','-nocrop','-m4')
    close all
end