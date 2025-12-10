function [tspan, x0, f,maxconc]=get_batch_condition(batch)

feeds={' Br1'	0.0015	0.0044
        'Br2'	0.0015	0.0256
        'Br3'	0.0085	0.0044
        'Br4'	0.0085	0.0256
        'Br5'	0.0000	0.0150
        'Br6'	0.0100	0.0150
        'Br7'	0.0050	0.0000
        'Br8'	0.0050	0.0300
        'Br9'	0.0050	0.0150};

tspan=0:12:240;
x0= [4    0    1.11    6.06   11.17   4.47   2.5    72.55   6.21   0.97  0.23   0.65   1.50   8.17  5.22   10.84  5.75  2.23  2.12  3.42  8.98  10.70  5.20  2.32  7.62   4.4e-3    6.7e-5    1.9e-7     3.1e-7     7e-6    6e-4   4.2e-8   5.6e-8  4e-7    2.2e-4   1.5e-6    8.9e-6   7.5e-7   9.5e-7    9.3e-7   2e-7    8.8e-7    5e-7  2.4e-8    8e-8  3E-9     0.270];
maxconc=[21.38 1179.39 16.16 27.38 11.17 4.47 3.23 72.55 21.01 2.31 36.96 22.63 10.19 44.71 22.36 26.56 5.75 3.45 37.94 22.68 39.40 66.62 18.01 2.32 53.33]';
mat=ccdesign(2); %9 experiments
miupreDoE=rescale(mat(1:9,1),0,0.01); %Rescale matrix elements between 0 and 0.01 (desired miu) Feed=miupre*X0*V0*exp(miupre*t)
FeedpostDoE=rescale(mat(1:9,2),0,0.03); %Rescale matrix elements between 0 and 0.03 (Feed rate post induction)
for i=batch
    miupre=miupreDoE(i);
    Feedpost=FeedpostDoE(i);
    X0=0.2;
    V0=2;
    f=(miupre*X0*V0*exp(miupre*[1:8]));
end
f(9:21)=feeds{1,3};

end