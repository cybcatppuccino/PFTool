default(parisizemax,1000000000);
thelist = [1024, 1080, 1152, 1200, 1215, 1280, 1296, 1350, 1440, 1458, 1536, 1600, 1620, 1728, 1800, 1920, 1944, 2025, 2048, 2160, 2187, 2304, 2400, 2430, 2560, 2592, 2700, 2880, 2916, 3072, 3200, 3240, 3456, 3600, 3645, 3840, 3888, 4050, 4096, 4320, 4374, 4608, 4800, 4860];
l = length(thelist);

counter = 716;
for (i=1, l, \
[N,k,chi] = [thelist[i], 2, Mod(1, thelist[i])]; \
mf = mfinit([N,k,chi],0);\
lf = mfeigenbasis(mf);\
for (j=1, length(lf), \
f = lf[j]; \
coef = mfcoefs(f, 12); \
if (type(coef[1]) <> "t_POLMOD", \
counter = counter + 1;\
print("L2[", counter, "]=", [N, j, coef, lfun(lfunmf(mf, f), 1)] ); \
););); 

counter = 1138;
for (i=1, l, \
[N,k,chi] = [thelist[i], 4, Mod(1, thelist[i])]; \
mf = mfinit([N,k,chi],0);\
lf = mfeigenbasis(mf);\
for (j=1, length(lf), \
f = lf[j]; \
coef = mfcoefs(f, 12); \
if (type(coef[1]) <> "t_POLMOD", \
counter = counter + 1;\
print("L4[", counter, "]=", [N, j, coef, lfun(lfunmf(mf, f), 1), lfun(lfunmf(mf, f), 2)] ); \
););); 