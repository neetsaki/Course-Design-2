distill=of(1:27);
strip=of(28:37);
dist_av=mean(distill);
str_av=mean(strip);
plot(2:1:38,of);hold on;
plot(2:1:28,dist_av*ones(27,1));hold on;
plot(29:1:38,str_av*ones(10,1))