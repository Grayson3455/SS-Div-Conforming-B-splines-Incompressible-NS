clc
clear
close all

dis = 0.02;
stepsize = 0.05;
max_number_vertices = 3000;
        
% well-resolved 
n = 128 + 1; % number of sample points

for Re = [7500, 10000]
    for p = [3]
        
        %no model - u
        ufile = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/u-xy.csv'];
        vfile = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/v-xy.csv'];
        uarray = csvread(ufile,1,0);
        varray = csvread(vfile,1,0);
        
        x = uarray(:,2);
        y = uarray(:,3);
        u = uarray(:,1);
        v = varray(:,1);
        
        % TRANSPOSE!!!!!
        X = reshape(x,[n,n]);
        X = X';
        Y = reshape(y,[n,n]);
        Y = Y';
        U = reshape(u,[n,n]);
        U = U';
        V = reshape(v,[n,n]);
        V = V';
        
        A =  figure('Renderer', 'painters', 'Position', [200 200 1000 1000]);
        hold on
        for loc = [0.005, 0.05]
        
            % line 1
            startx = 0:dis:1;
            starty = zeros(size(startx))+loc ;

            h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
            set(h,'color','black')
            set(gca,'FontSize',40)

            % line 2
            starty = 0:dis:1;
            startx = zeros(size(starty))+loc ;

            h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
            set(h,'color','black')
            set(gca,'FontSize',40)

            % line 3
            starty = 0:dis:1;
            startx = zeros(size(starty))+(1-loc) ;

            h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
            set(h,'color','black')
            set(gca,'FontSize',40)
            set(gca,'XTick',[], 'YTick', [])
            set(gca,'XColor', 'none','YColor','none')
            
            hold on
        end

        savename = ['fig/Best-',num2str(Re),'.png'];
        saveas(A,savename)
        
    end
end

% change discretization due to different resolution
dis = 0.02;
stepsize = 0.0025;
max_number_vertices = 8000;

% poor resolved
for Re = [7500, 10000]
    for p = [1,2,3]
        for n = [16+1]
        
            %no model - u
            ufile = ['NM_results/nm_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/u-xy.csv'];
            vfile = ['NM_results/nm_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/v-xy.csv'];
            uarray = csvread(ufile,1,0);
            varray = csvread(vfile,1,0);

            x = uarray(:,2);
            y = uarray(:,3);
            u = uarray(:,1);
            v = varray(:,1);

            % TRANSPOSE!!!!!
            X = reshape(x,[n,n]);
            X = X';
            Y = reshape(y,[n,n]);
            Y = Y';
            U = reshape(u,[n,n]);
            U = U';
            V = reshape(v,[n,n]);
            V = V';

            A =  figure('Renderer', 'painters', 'Position', [200 200 1000 1000]);
            hold on

            for loc = [0.005, 0.05]

                % line 1
                startx = 0:dis:1;
                starty = zeros(size(startx))+loc ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)

                % line 2
                starty = 0:dis:1;
                startx = zeros(size(starty))+loc ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)

                % line 3
                starty = 0:dis:1;
                startx = zeros(size(starty))+(1-loc) ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)
                set(gca,'XTick',[], 'YTick', [])
                set(gca,'XColor', 'none','YColor','none')
            end

            savename = ['fig/NM-',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'.png'];
            saveas(A,savename)
        end
        
    end
end

% stabilized
for Re = [7500, 10000]
    for p = [1,2,3]
        for n = [16+1]
        
            %no model - u
            ufile = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/u-xy.csv'];
            vfile = ['SS_results/ss_Re=',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'/v-xy.csv'];
            uarray = csvread(ufile,1,0);
            varray = csvread(vfile,1,0);

            x = uarray(:,2);
            y = uarray(:,3);
            u = uarray(:,1);
            v = varray(:,1);

            % TRANSPOSE!!!!!
            X = reshape(x,[n,n]);
            X = X';
            Y = reshape(y,[n,n]);
            Y = Y';
            U = reshape(u,[n,n]);
            U = U';
            V = reshape(v,[n,n]);
            V = V';

            A =  figure('Renderer', 'painters', 'Position', [200 200 1000 1000]);
            hold on

            for loc = [0.005, 0.05]

                % line 1
                startx = 0:dis:1;
                starty = zeros(size(startx))+loc ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)

                % line 2
                starty = 0:dis:1;
                startx = zeros(size(starty))+loc ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)

                % line 3
                starty = 0:dis:1;
                startx = zeros(size(starty))+(1-loc) ;

                h= streamline(X,Y,U,V,startx,starty,[stepsize, max_number_vertices]);
                set(h,'color','black')
                set(gca,'FontSize',40)
                set(gca,'XTick',[], 'YTick', [])
                set(gca,'XColor', 'none','YColor','none')
            end

            savename = ['fig/SS-',num2str(Re),'_p=',num2str(p),'_m=',num2str(n-1),'.png'];
            saveas(A,savename)
        end
        
    end
end

