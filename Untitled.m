y=[0.053 0.51 0.77 1.57 2.9 3.725 4.5 8.76 16.34 29.92 39.16 47.49 51.67 55.74 56.44 58.78 62.22 64.7 66.28 70.29 72.71 74.69 76.93 79.26 81.83 84.91 86.4 89.41]/100;
x=[0.004 0.04 0.05 0.12 0.23 0.31 0.39 0.79 1.61 4.16 7.41 12.64 17.41 25.75 27.3 33.24 42.09 48.92 52.68 61.02 65.64 68.92 72.36 75.99 79.82 83.87 85.97 89.41]/100;
y_idea=0:0.1:1;
x_idea=y_idea;
x_F=0.115385;
x_D=0.859756;
delta=0.0001;
x_full=0:delta:1;
y_full=interp1(x,y,x_full,'Pchip');

plot(x,y);hold on;
plot(x_idea,y_idea);hold on;
scatter(x_D,x_D);hold on;plot(x_F*ones(11,1),0:0.1:1);hold on;
plot(x_full,y_full);hold on;

x_q=x_F;
y_q=y_full(1154);
%plot([x_q x_D],[y_q x_D]);
error=[];
for i=1:1:0.8/delta

    slope=(y_full(i+1)-y_full(i))/delta;
    x_t=delta*i;y_t=y_full(i);
    e=(x_D-y_t)/(x_D-x_t)-slope;
    error=[error e];

end
% xlim([0 1]);
% plot(1:length(error),error)
x_linear=0:delta:1;
y_linear=(x_D-0.8100)/(x_D-0.7863)*(x_linear-0.7863)+0.81;
plot(x_linear,y_linear);
R_min=(x_D-0.3556)/(0.3556-0.1154);

%plot(0.0001:0.0001:0.8,error)














