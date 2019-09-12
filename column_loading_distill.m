x_min=1.6244;
x_max=12.5214;
y_max=3728.81;
y_min=1010.806;
A=1.19044e-8;
B=6.09907e-5;
C=0.2838-0.141;
x=0:0.1:(x_max+1);
y=sqrt((C-B*x.^2)/A);
plot(x_min*[1 1],[y_min y_max]);hold on;
plot(x_max*[1 1],[y_min y_max]);hold on;
plot([x_min x_max],y_min*[1 1]);hold on;
plot([x_min x_max],y_max*[1 1]);hold on;
plot(x,y);hold on;
scatter(3.58460351978,3454.2487412474966)