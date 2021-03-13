clc
clear
close all

% integrate files
for p = [1,2,3]
    ss_p = zeros(6,2);
    count = 1;
    for m = [4,8,16,32,64,128]
        filename = ['error/ss_p=',num2str(p),'_Re=10_m=',num2str(m),'.csv'];
        error = csvread(filename, 0,0);
        ss_p(count,1) = error(1);
        ss_p(count,2) = error(2);
        count = count + 1;
    end 
    savename = ['ss_p',num2str(p),'.csv'];
    csvwrite(savename,ss_p)
end

for p = [1,2,3]
    file_p1 = ['ss_p',num2str(p),'.csv'];    
    phi_p1 = csvread(file_p1,0,0);
    for h = [1,2,3,4,5,6]
        if h == 1 
            res = 1/4;
            L2order = NaN;
            H1order = NaN;
        elseif h == 2
            res = 1/8;
            L2order = log( -phi_p1(h,1) + phi_p1(h-1,1) ) / log(-1/8+1/4);
            H1order = log( -phi_p1(h,2) + phi_p1(h-1,2) ) / log(-1/8+1/4);
        elseif h == 3
            res = 1/16;
            L2order = log( -phi_p1(h,1) + phi_p1(h-1,1) ) / log(-1/16+1/8);
            H1order = log( -phi_p1(h,2) + phi_p1(h-1,2) ) / log(-1/16+1/8);
        elseif h ==4
            res = 1/32;
            L2order = log( -phi_p1(h,1) + phi_p1(h-1,1) ) / log(-1/32+1/16);
            H1order = log( -phi_p1(h,2) + phi_p1(h-1,2) ) / log(-1/32+1/16);
        elseif h == 5
            res = 1/64;
            L2order = log( -phi_p1(h,1) + phi_p1(h-1,1) ) / log(-1/64+1/32);
            H1order = log( -phi_p1(h,2) + phi_p1(h-1,2) ) / log(-1/64+1/32);
        elseif h == 6
            res = 1/128;
            L2order = log( -phi_p1(h,1) + phi_p1(h-1,1) ) / log(-1/128+1/64);
            H1order = log( -phi_p1(h,2) + phi_p1(h-1,2) ) / log(-1/128+1/64);
        end
        
        fprintf('p=%d, h = %.3f, L2 = %.3e, H1 = %.3e, L2order= %.3f, H1order= %.3f\n', p , res, phi_p1(h,1), phi_p1(h,2), L2order, H1order)
    end
    fprintf('\n')
end