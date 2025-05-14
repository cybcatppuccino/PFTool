default(parisizemax,1000000000);
thelist = [10,20,30,40,50,60,70,80,90,100,110,120,140,150,160,180,200];
l = length(thelist);


counter = 0;
for (i=1, l, \
[N,k,chi] = [thelist[i], 4, Mod(1, thelist[i])]; \
mf = mfinit([N,k,chi],0);\
lf = mfeigenbasis(mf);\
for (j=1, length(lf), \
f = lf[j]; \
coef = mfcoefs(f, 89); \
if (type(coef[1]) <> "t_POLMOD", \
counter = counter + 1;\
print("L4[", counter, "]=", [N, j, coef] ); \
););); 