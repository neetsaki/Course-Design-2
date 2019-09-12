x_min=1.6244450966154191;
x_max=12.52143;
y_max=3491.522882224483;
y_min=1191.9051890434162;
A=7.753923008455667e-9;
B=6.099069584024051e-5;
C=0.2833-0.141;
x=0:0.1:(x_max+1);
y=sqrt((C-B*x.^2)/A);
plot(x_min*[1 1],[y_min y_max]);hold on;
plot(x_max*[1 1],[y_min y_max]);hold on;
plot([x_min x_max],y_min*[1 1]);hold on;
plot([x_min x_max],y_max*[1 1]);hold on;
plot(x,y);hold on;
scatter(6.033082537782988,3261.6382560342036)