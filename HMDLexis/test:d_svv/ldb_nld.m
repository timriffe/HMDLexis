pwd
addpath('/data/commons/boe/HMD_Lexis_mat2R/MATLAB.hg');
% NLD is test country case, so pick up data from there
addpath('/data/wilmoth0/HMD/HMDWORK/NLD/InputDB');

%% run
disp('Read...');

indb_input('NLD');

%% load Input DB into Matlab workspace 
load NLD;

%[deaths, population, births]=indb_input('NLD');

deaths=deaths(:,1:end-1);
population=population(:,1:end-1);
births=births(:,1:end-1);

deaths=d_ma0(deaths);

disp('split by triangles');
deaths=d_s1x1(deaths, births);

%% pre
% save deaths before and after the test point

save('NLDdeaths_pre.mat','deaths');

%% split some vv deaths counts into Lexis triangles
deaths=d_svv(deaths);

% save deaths after
save('NLDdeaths_post.mat','deaths');

%%--------------------------------

%open age interval
disp('open age interval');
deaths=d_soainew(deaths);
deaths=d_long(deaths);

%distribution of unknown 
deaths=d_unk(deaths);

disp('population');
population=p_s5x1(population,deaths);
population=p_unk(population);
population=p_ey2ny(population);

population=p_precensal(population,deaths,[1850,1850]);

disp('Extinct cohort method');
population=p_srecm(population, deaths);

ldb_output(deaths, population, 'mNLD.txt', 'fNLD.txt', births);

d_printRA('NLD','Netherlands');

ldb_psum('mNLD.txt', 'totmales');
ldb_psum('fNLD.txt', 'totfemales');