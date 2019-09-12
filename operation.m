y=[0.053 0.51 0.77 1.57 2.9 3.725 4.5 8.76 16.34 29.92 39.16 47.49 51.67 55.74 56.44 58.78 62.22 64.7 66.28 70.29 72.71 74.69 76.93 79.26 81.83 84.91 86.4 89.41]/100;
x=[0.004 0.04 0.05 0.12 0.23 0.31 0.39 0.79 1.61 4.16 7.41 12.64 17.41 25.75 27.3 33.24 42.09 48.92 52.68 61.02 65.64 68.92 72.36 75.99 79.82 83.87 85.97 89.41]/100;
y_idea=0:0.1:1;
x_idea=y_idea;
x_F=0.115385;
x_D=0.859756;
delta=0.0001;
x_full=0:delta:1;
y_full=interp1(x,y,x_full,'Pchip');
x_B=0.0003915426781519186;
q=1;

plot(x,y);hold on;
plot(x_idea,y_idea);hold on;
scatter(x_D,x_D);hold on;plot(x_F*ones(11,1),0:0.1:1);hold on;
plot(x_full,y_full);hold on;

%operation line
x_o=0:0.01:1;
R=3.7161;
y_o_dist=R/(R+1)*x_o+x_D/(R+1);
R_prime=(x_F-x_B)/(x_D-x_F)*(R+q)+q;
y_o_ex=R_prime/(R_prime-1)*x_o-1/(R_prime-1)*x_B;
plot(x_o,y_o_dist);hold on;plot(x_o,y_o_ex)