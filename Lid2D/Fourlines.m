% Create figure with multiple axes
% Source : https://www.mathworks.com/help/matlab/creating_plots/graph-with-multiple-x-axes-and-y-axes.html
clc
clear
close all

name_ghia = 'Ghia(1982)';
nnn       = 16; % resolution for no-model

for Re = [7500,10000]
    for p = [1,2,3]
        
        %--Ghia's----%
        uG = csvread('upaper.csv',1,0);
        uG_y    = uG(:,2);
        if Re ==7500
            u_ghia = uG(:,8);
            paper_y = uG_y;
            p_u = u_ghia;
        elseif Re == 10000
            u_ghia = uG(:,9);
            
            % wrong point in Ghia's
            paper_y = uG_y;
            paper_y(9) = [];
            p_u = u_ghia;
            p_u(9) = [];
            
        end
        
        %--------------u----------------%
        % stabilized 128x128 
        refu1 = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(3),'_m=',num2str(128),'/u_y.csv'];
        uref1 = csvread(refu1,1,0);
        ref_y1 = uref1(:,end-1);
        ref_u1 = uref1(:,1);
        name_refu1 = 'Well-resolved';
        
        % no model - u 16x16
        refu2 = ['NM_results/nm_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(nnn),'/u_y.csv'];
        uref2 = csvread(refu2,1,0);
        ref_y2 = uref2(:,end-1);
        ref_u2 = uref2(:,1);
        name_refu2 = 'Unstabilized';
        
           
        nsu = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(nnn),'/u_y.csv'];
        uns = csvread(nsu,1,0);
        ns_y = uns(:,end-1);
        ns_u = uns(:,1);
        name_uns = 'Stabilized';
        
        %-----------v-------------%
        
        
         %--Ghia's----%
        vG = csvread('vpaper.csv',1,0);
        vG_x    = vG(:,2);
        if Re ==7500
            v_ghia = vG(:,8);
        elseif Re == 10000
            v_ghia = vG(:,9);
        end
        
        
        %stabilized 128x128 
        refv1 = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(3),'_m=',num2str(128),'/v_x.csv'];
        vref1 = csvread(refv1,1,0);
        ref_x1 = vref1(:,end-2);
        ref_v1 = vref1(:,1);
        name_refv1 = 'Well-resolved';
       
        % no model - v
        refv2 = ['NM_results/nm_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(nnn),'/v_x.csv'];
        vref2 = csvread(refv2,1,0);
        ref_x2 = vref2(:,end-2);
        ref_v2 = vref2(:,1);
        name_refv2 = 'Unstabilized';
        
        % stabilized - v
        nsv = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(nnn),'/v_x.csv'];
        vns = csvread(nsv,1,0);
        ns_x = vns(:,end-2);
        ns_v = vns(:,1);
        name_vns = 'Stabilized';
        

        A =  figure('Renderer', 'painters', 'Position', [0 0 1350 1350]);
        
        line(vG_x,v_ghia,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'LineWidth',3, 'DisplayName',name_ghia)
        line(ref_x1,ref_v1,'Color','b','LineStyle',':','LineWidth',3,'DisplayName',name_refv1)
        line(ref_x2,ref_v2,'Color','b','LineStyle','-.','LineWidth',3,'DisplayName',name_refv2)
        line(ns_x,ns_v,'Color','b','LineStyle','-','LineWidth',3,'DisplayName',name_vns)
        
        ax1 = gca; % current axes
        ax1.XColor = 'b';
        ax1.YColor = 'b';
        xlabel(ax1,'$x$','interpreter','latex')
        ylabel(ax1,'$v$','interpreter','latex')
        xlim(ax1,[0 1])
        ylim(ax1,[-0.7 0.5])
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times')
        set(gca,'FontSize',34)
        legend('FontSize',40,'Location','north','interpreter','latex')
       

        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','None');
        
        line(paper_y,p_u,'Color','k','LineStyle','none','Marker','o','MarkerSize',12,'LineWidth',3,'DisplayName',name_ghia)
        line(ref_y1,ref_u1,'Parent',ax2,'Color','k','LineStyle',':','LineWidth',3,'DisplayName',name_refu1)
        line(ref_y2,ref_u2,'Parent',ax2,'Color','k','LineStyle','-.','LineWidth',3,'DisplayName',name_refu2)
        line(ns_y,ns_u,'Parent',ax2,'Color','k','LineStyle','-','LineWidth',3,'DisplayName',name_uns)
        xlabel(ax2,'$y$','interpreter','latex')
        ylabel(ax2,'$u$','interpreter','latex')
        xlim(ax2,[0 1])
        ylim(ax2,[-.5 1])
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times')
        set(gca,'FontSize',34)
        legend('FontSize',40,'Location','south','interpreter','latex')
        
        savename = ['fig/Re=',num2str(Re),'_p=',num2str(p),'_n=',num2str(nnn),'.png'];
        saveas(A,savename)
    end
end