clc
clear 
close all

u_name = ['u_y.csv'];
u = csvread(u_name,1,0);

v_name = ['v_x.csv'];
v = csvread(v_name,1,0);

plot(u(:,end-1),u(:,1))
