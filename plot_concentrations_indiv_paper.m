function plot_concentrations_indiv_paper(simfile,nrow,batchs)
load([simfile '.mat']);
ci=1;
names=data.ylabel;names2=names;

divi=1;divie=numel(names);

Names=names;Names2=names2;Names=replace(Names,{'Cells' 'Nh4' 'val' 'Gly'},{'Xv  ' 'NH4' 'Val' 'Glyc'});

lis='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
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
        16-(16-1)/4*4,24-(24-1.24)/7*1,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*4,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];

Pos2=[  16-(16-1)/4*3,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*1,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*3,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];

Pos3=[  16-(16-1)/4*2,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*1,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*2,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];
    
Pos4=[  16-(16-1)/4*1,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*1,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*1,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];
  
  Pos5=[16-(16-1)/4*0,24-(24-1.24)/7*0,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*1,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*2,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*3,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*4,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*5,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*6,2.60849056603774,2.07222222222222;
        16-(16-1)/4*0,24-(24-1.24)/7*7,2.60849056603774,2.07222222222222];


Pos10=Pos1;Pos20=Pos2;Pos30=Pos3;Pos40=Pos4;Pos50=Pos5;
% normalize conc
y=[];ysim=[];
for ib=batchs
    y=[y;data.batch(ib).y];ysim=[ysim;data.batch(ib).ysim];
end
for ib=batchs
    data.batch(ib).y=data.batch(ib).y./max(y);
    data.batch(ib).ysim=data.batch(ib).ysim./max(ysim);
end

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
        'position', [0.5 1.5 17.5*0.86 27.2*0.85]);
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
    np=numel(names);
    train=[];test=[];
    for ib=batchs
        if data.batch(ib).istrain==1
            train=[train ib];
        else
            test=[test ib];
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
                
                maxy1=-inf;miny1=inf; maxx1=-inf;minx1=inf;
                x1=[];x2=[];y1=[];y2=[];
                for ib=batchs
                    x = data.batch(ib).y(:,id);
                    y = data.batch(ib).ysim(:,id);
                    maxy1=max([maxy1;y]);miny1=min([miny1;y]);maxx1=max([maxx1;x]);minx1=min([minx1;x]);
                    if data.batch(ib).istrain==1
                        x1=[x1;x];y1=[y1;y];
                        trainingPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#0080ff','Marker','.','LineStyle','none','MarkerSize',12,'HandleVisibility','off');
                    elseif data.batch(ib).istrain==2
%                          testPlot = plot(x,y,...
%                             'MarkerSize',5,'Color','#2bbf5c','Marker','x','LineStyle','none','MarkerSize',10);
                    end
                end
                maxys(row)=maxy1;maxxs(row)=maxx1;
                
                if ~isempty (y1)
                    mod=fitlm(x1,y1,'linear','RobustOpts','off');
                    R2_tr=mod.Rsquared;
                end
                
%                 plot(x1(1),y1(2),'Marker','.','LineStyle','none','MarkerSize',12,'Color','#0080ff', 'HandleVisibility','on');
                                
                
                for ib=batchs
                    x = data.batch(ib).y(:,id);
                    y = data.batch(ib).ysim(:,id);
%                     maxy1=max([maxy1;y]);miny1=min([miny1;y]);maxx1=max([maxx1;x]);minx1=min([minx1;x]);
                    if data.batch(ib).istrain==1
%                         trainingPlot = plot(x,y,...
%                             'MarkerSize',5,'Color','#0080ff','Marker','.','LineStyle','none','MarkerSize',12);
                    elseif data.batch(ib).istrain==2
                        x2=[x2;x];y2=[y2;y];
                         testPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#2bbf5c','Marker','.','LineStyle','none','MarkerSize',6,'HandleVisibility','off');
                    end
                end
                if ~isempty(y2)
%                     plot(x2(1),y2(2),'Marker','.','LineStyle','none','MarkerSize',12,'Color','#2bbf5c', 'HandleVisibility','on');
                    mod=fitlm(x2,y2,'linear','RobustOpts','off');
                    R2_tt=mod.Rsquared;
                else
                    R2_tt.Ordinary=inf;
                end
                
                                
                if id==np
                    il=[];
                    if ~isempty(train)
                        il=[il 1];
                        x = data.batch(ib).y(:,id);
                        y = data.batch(ib).ysim(:,id);
                        trainingPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#0080ff','Marker','.','LineStyle','none','MarkerSize',12,'HandleVisibility','on');
                    end
                    if ~isempty(test)
                        il=[il 2];
                        x = data.batch(ib).y(:,id);
                        y = data.batch(ib).ysim(:,id);
                        testPlot = plot(x,y,...
                            'MarkerSize',5,'Color','#2bbf5c','Marker','.','LineStyle','none','MarkerSize',12,'HandleVisibility','on');
                    end
                    lg={'train' 'test'};
                    legend(lg(il),'Location', 'north','FontSize',7,'FontName','Arial')
                end
                
                
%                 miny1=miny1-0.1;maxy1=maxy1+0.1;
%                 minx1=minx1-0.1;maxx1=maxx1+0.1;
                miny1=0-0.15;maxy1=1+0.15;
                minx1=0-0.15;maxx1=1+0.15;
                if strcmp(Names(id),'Cys')
%                     1;
%                     miny1=-0.35;maxy1=1.16;
%                     minx1=-0.36;maxx1=1.2;
                elseif strcmp(Names(id),'Pro')
%                     1;
%                      miny1=1.5;maxy1=5.5;
%                     minx1=1.5;maxx1=5.5;
                end

                ylim(hAx1,[ min([miny1 minx1]) max([maxy1 maxx1])*1.])
                xlim(hAx1,[ min([miny1 minx1]) max([maxy1 maxx1])]*1.)
                line([ min([miny1 minx1]) max([maxy1 maxx1])]*1.,[ min([miny1 minx1]) max([maxy1 maxx1])*1.],'color',[.75 0.75 0.75],'HandleVisibility','off');
    
    
                axis square
%                 if np-id>nrow-1
                set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
%                 else
%                     set(hAx1,'Units','centimeters','position',position,'Box','on','YMinorGrid','off','XMinorGrid','off','FontSize',7,'FontName','Arial');hold on;
%                     
%                 end
                
                hAx = gca;             % handle to current axes
                n=names2{id};
                set(gca,'FontSize',7);
                yrange=max([maxy1 maxx1])-min([miny1 minx1]);
                xrange=max([maxy1 maxx1])-min([miny1 minx1]);
%                 text(min([miny1 minx1])+xrange*0.75,min([miny1 minx1])+yrange*1,[lis((i)) num2str(row) '  '],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8,'FontName','Arial','FontWeight','bold')
                xlabel('Experimental','FontSize',8,'FontName','Arial','FontWeight','bold');
                ylabel('Predicted','FontSize',8,'FontName','Arial','FontWeight','bold');
                text(min([miny1 minx1])*1.1+xrange*0.95,min([miny1 minx1])+yrange*0.2,Names(id),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',9,'FontName','Arial','FontWeight','bold')
                if ~isempty (y1)
                    plot(min([miny1 minx1])+xrange*0.06,min([miny1 minx1])+yrange*0.9, 'MarkerSize',5,'Color','#0080ff','Marker','.','LineStyle','none','MarkerSize',12,'HandleVisibility','off')
                    text(min([miny1 minx1])+xrange*0.55,min([miny1 minx1])+yrange*0.97,['R^{2}: ' num2str(R2_tr.Ordinary,'%.2f') ],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6.,'FontName','Arial','FontWeight','bold')
                end
                if ~isempty (y2)
                    plot(min([miny1 minx1])+xrange*0.06,min([miny1 minx1])+yrange*0.75,'MarkerSize',5,'Color','#2bbf5c','Marker','.','LineStyle','none','MarkerSize',12,'HandleVisibility','off')
                    text(min([miny1 minx1])+xrange*0.55,min([miny1 minx1])+yrange*0.82,['R^{2}: ' num2str(R2_tt.Ordinary,'%.2f') ],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6.,'FontName','Arial','FontWeight','bold')
                end
               try; patchline([min([miny1 minx1]) max([maxy1 maxx1]) ], [min([miny1 minx1]) max([maxy1 maxx1])],'edgecolor','k','linewidth',1,'edgealpha',0.1,'HandleVisibility','off');catch;end
                set(gca,'FontWeight','bold')
            end
%             if ic>1 && i==5
%                 6;
%             end
%             ids=(i-1)*row+1:(i-1)*row+row;
                       
        end
    catch kl
        1
    end
    
end
    
