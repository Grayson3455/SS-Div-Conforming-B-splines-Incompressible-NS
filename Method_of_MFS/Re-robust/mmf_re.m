clc
clear
close all

% integrate files
for p = [1,2,3]
    ss_p = zeros(4,2);
    count = 1;
    for Re = [1,10,100,1000]
        filename = ['error/ss_p=',num2str(p),'_Re=',num2str(Re),'_m=16.csv'];
        error = csvread(filename, 0,0);
        ss_p(count,1) = error(1);
        ss_p(count,2) = error(2);
        count = count + 1;
    end 
    savename = ['ss_p',num2str(p),'.csv'];
    csvwrite(savename,ss_p)
end

% get figures
for p = [1,2,3]
    A =  figure('Renderer', 'painters', 'Position', [200 100 800 800]);
    error = csvread(['ss_p',num2str(p),'.csv'], 0,0);
    namel2 = '$L^2$';
    nameh1 = '$\mathbf{H}^1$';
    loglog([1,10,100,1000],error(:,1),'-ko','MarkerSize',20,'DisplayName',namel2,'LineWidth',3);
    hold on
    loglog([1,10,100,1000],error(:,2),'-kx','MarkerSize',20,'DisplayName',nameh1,'LineWidth',3);
    ylim([1e-8 1e-1])
    xlabel('\it{Re}','interpreter','latex')
    ylabel('Error', 'FontName', 'Times')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times')
    set(gca,'FontSize',40)
    legend('FontSize',40,'Location','southeast','interpreter','latex')
    savefig = ['fig/MFS-Re-p',num2str(p),'.png'];
    saveas(A,savefig)
end
    
