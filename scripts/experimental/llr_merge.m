pp1=.5;
pm1=.1;

pp2=.99;
pm2=.001;

L1 = log(pp1/pm1)
L2 = log(pp2/pm2)

Lres = log((pp1+pp2)/(pm1+pm2))

Lguess = 2*atanh((tanh(L1/2)+tanh(L2/2))/2)

Larith = mean([L1 L2])
Lgeo = geomean([L1 L2])

%%
pp1= linspace(0,1);
pm1= 

pp2=.99;
pm2=.001;