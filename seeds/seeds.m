clear
clf
cla

N = 200; 


% MILLISECONDS + RANDI WITH PAUSE
seeds_milli_p = zeros(N,1);
firsts_milli_p = zeros(N,1);
for i=1:N
    seed = sscanf(datestr(now, 'FFF'),'%d') * 1000;
    s = RandStream('mcg16807','Seed', seed1);
    RandStream.setGlobalStream(s);
    firsts_milli_p(i) = randi(100);
    seeds_milli_p(i) = RandStream.getGlobalStream.Seed;
    pause(rand(1));
end

plot(1:N, seeds_milli_p,'rx');

autocorr(seeds_milli_p)

autocorr(firsts_milli_p)

plot(1:N, firsts_milli_p, 'rx');




% MILLISECONDS + RANDI WITHOUT PAUSE
seeds_milli = zeros(N,1);
firsts_milli = zeros(N,1);
for i=1:N
    seed = sscanf(datestr(now, 'FFF'),'%d') * 1000;
    s = RandStream('mcg16807','Seed', seed1);
    RandStream.setGlobalStream(s);
    firsts_milli(i) = randi(100);
    seeds_milli(i) = RandStream.getGlobalStream.Seed;
end

plot(1:N, seeds_milli,'rx');

autocorr(seeds_milli)

autocorr(firsts_milli)

plot(1:N, firsts_milli, 'rx');

% MIXED MILLISECONDS + RANDI WITH PAUSE
seedss = zeros(N,1);
firsts = zeros(N,1);
for i=1:N
    seed1 = sscanf(datestr(now, 'FFF'),'%d') * 1000;
    s = RandStream('mcg16807','Seed', seed1);
    RandStream.setGlobalStream(s);
    seed2 = randi(1000);
    seed = seed1 + seed2;
    s = RandStream('mcg16807','Seed', seed);
    RandStream.setGlobalStream(s);
    firsts(i) = randi(100);
    seedss(i) = RandStream.getGlobalStream.Seed;
    pause(rand(1));
end


plot(1:N, seedss,'rx');

autocorr(seedss)

autocorr(firsts)

plot(1:N, firsts, 'rx');

[h,p,Qstat,crit] = lbqtest(firsts,'Lags',[1,5,6,10,15])
mean(firsts)
std(firsts)

% MIXED MILLISECONDS + RANDI WITHOUT PAUSE
seedss_without_p = zeros(N,1);
firsts_without_p = zeros(N,1);
for i=1:N
    seed1 = sscanf(datestr(now, 'FFF'),'%d') * 1000;
    s = RandStream('mcg16807','Seed', seed1);
    RandStream.setGlobalStream(s);
    seed2 = randi(1000);
    seed = seed1 + seed2;
    s = RandStream('mcg16807','Seed', seed);
    RandStream.setGlobalStream(s);
    firsts_without_p(i) = randi(100);
    seedss_without_p(i) = RandStream.getGlobalStream.Seed;
end

plot(1:N, seedss_without_p,'rx');

plot(1:N, firsts_without_p, 'rx');

autocorr(seedss_without_p)

autocorr(firsts_without_p)


% RNG SHUFFLE WITH FIXED PAUSE 0.5
seeds_rngshuffle_p_fixed = zeros(N,1);
firsts_rngshuffle_p_fixed = zeros(N,1);
for i=1:N
    rng shuffle
    firsts_rngshuffle_p_fixed(i) = randi(100);
    seeds_rngshuffle_p_fixed(i) = RandStream.getGlobalStream.Seed;
    pause(0.5);
end
plot(1:N, seeds_rngshuffle_pfixed,'rx');

plot(1:N, firsts_rngshuffle_pfixed, 'rx');

autocorr(seeds_rngshuffle_pfixed)

autocorr(firsts_rngshuffle_pfixed)


% RNG SHUFFLE WITHOUT PAUSE
seeds_rngshuffle = zeros(N,1);
firsts_rngshuffle = zeros(N,1);
for i=1:N
    rng shuffle
    firsts_rngshuffle(i) = randi(100);
    seeds_rngshuffle(i) = RandStream.getGlobalStream.Seed;
end

plot(1:N, seeds_rngshuffle,'rx');

plot(1:N, firsts_rngshuffle, 'rx');

autocorr(seeds_rngshuffle)

autocorr(firsts_rngshuffle)

[h,p,Qstat,crit] = lbqtest(firsts_rngshuffle,'Lags',[1,5,10,15])

% RNG SHUFFLE WITH RANDOM PAUSE
seeds_rngshuffle_p = zeros(N,1);
firsts_rngshuffle_p = zeros(N,1);
for i=1:N
    rng shuffle
    firsts_rngshuffle_p(i) = randi(100);
    seeds_rngshuffle_p(i) = RandStream.getGlobalStream.Seed;
    pause(rand(1));
end

plot(1:N, seeds_rngshuffle_p,'rx');

plot(1:N, firsts_rngshuffle_p, 'rx');

autocorr(seeds_rngshuffle_p)

autocorr(firsts_rngshuffle_p)

[h,p,Qstat,crit] = lbqtest(firsts_rngshuffle_p,'Lags',[1,5,9,10,15])
mean(firsts_rngshuffle_p)
std(firsts_rngshuffle_p)

% NO SHUFFLE NO PAUSE
seeds_simple = zeros(N,1);
firsts_simple = zeros(N,1);
rng shuffle
for i=1:N
    firsts_simple(i) = randi(100);
    seeds_simple(i) = RandStream.getGlobalStream.Seed;
end

plot(1:N, firsts_simple, 'rx');

autocorr(firsts_simple)
[h,p,Qstat,crit] = lbqtest(firsts_simple,'Lags',[1,5,9,10,15])
mean(firsts_simple)
std(firsts_simple)

% NO SHUFFLE NO PAUSE - NEW SEED generated from the generator
seeds_regen = zeros(N,1);
firsts_regen = zeros(N,1);
rng shuffle
for i=1:N
    s = RandStream('mcg16807','Seed', randi(1));
    RandStream.setGlobalStream(s);
    firsts_regen(i) = randi(100);
    seeds_regen(i) = RandStream.getGlobalStream.Seed;
end

plot(1:N, firsts_regen, 'rx');

autocorr(firsts_regen)
[h,p,Qstat,crit] = lbqtest(firsts_regen,'Lags',[1,5,9,10,15])
mean(firsts_regen)
std(firsts_regen)

% It seems that generating shuffling the seed inside the loop is quite
% bad if we do not introduce a delay between iterations. Seeds are too
% close with each others. The same it is true for using milliseconds or
% clock.

% Using 1 random initial seed outside of the loop generates, and not
% re-initialized the generator in each loop brings better results. But 
% this is not possible in the cluster.





