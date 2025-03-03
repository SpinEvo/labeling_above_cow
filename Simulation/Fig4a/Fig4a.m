load('E_array');
fig = figure('Position', [0, 0, 500, 250]);

method_color = [1 0 0;
            0 0 1;
            0.6 0.9 0.2];
trans = 0.5;

subplot(131);
marker = {'o','^','v','s','d','p','h'};
x = [0 1];
x = categorical(x);
%RF_sep = [1560,1380,1200];
for i = 1
    a = [E_array(:,1) E_array(:,2)];
    for subID = 1:7
        
        scatter(x(1),a(subID,1),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(1,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans);hold on;
        scatter(x(2),a(subID,2),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(2,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans);hold on;
        
        line([x(1),x(2)],[a(subID,1),a(subID,2)],'Color',[0 0 0],'LineWidth',1.5);
    end
    box off
    yticks(0.3:0.6:0.9); yticklabels({'0.3','0.9'});ylim([0.3 0.92]);
    xticklabels({'OES / Default','IOES / Default'});
    %yticks(0.1:0.8:0.9); yticklabels({'0.1','0.9'});ylim([0.1 0.92]);
    
    set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');

end

subplot(132);
marker = {'o','^','v','s','d','p','h'};
x = [0 1];
x = categorical(x);
%RF_sep = [1560,1380,1200];
for i = 1
    a = [E_array(:,2) E_array(:,3)];
    for subID = 1:7
        scatter(x(1),a(subID,1),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(2,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans);hold on;
        scatter(x(2),a(subID,2),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(3,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans+0.2);hold on;
        line([x(1),x(2)],[a(subID,1),a(subID,2)],'Color',[0 0 0],'LineWidth',1.5);
    end
    box off
    yticks(0.3:0.6:0.9); yticklabels({'0.3','0.9'});ylim([0.3 0.92]);
    xticklabels({'IOES / Default','IOES / Optimal'});
    %yticks(0.1:0.8:0.9); yticklabels({'0.1','0.9'});ylim([0.1 0.92]);
    
    set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');

end

subplot(133);
marker = {'o','^','v','s','d','p','h'};
x = [0 1];
x = categorical(x);
%RF_sep = [1560,1380,1200];
for i = 1
    a = [E_array(:,1) E_array(:,3)];
    for subID = 1:7
        scatter(x(1),a(subID,1),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(1,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans);hold on;
        scatter(x(2),a(subID,2),'x','Marker',marker{subID},'MarkerFaceColor',[method_color(3,:)],'MarkerEdgeColor','none','MarkerFaceAlpha', trans+0.2);hold on;
        line([x(1),x(2)],[a(subID,1),a(subID,2)],'Color',[0 0 0],'LineWidth',1.5);
    end
    box off
    yticks(0.3:0.6:0.9); yticklabels({'0.3','0.9'});ylim([0.3 0.92]);
    xticklabels({'OES / Default','IOES / Optimal'});
    %yticks(0.1:0.8:0.9); yticklabels({'0.1','0.9'});ylim([0.1 0.92]);
    
    set(gca,'linewidth',1.5,'FontSize',12,'FontWeight','bold');
end