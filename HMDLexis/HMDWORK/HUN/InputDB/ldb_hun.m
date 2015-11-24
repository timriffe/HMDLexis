[deaths, population, births]=indb_input('HUN');
% load Input DB into Matlab workspace 

deaths=deaths(:,1:end-1);
population=population(:,1:end-1);
births=births(:,1:end-1);

population=selif(population,(population(:,9)==0 & population(:,5)>1959) |  (population(:,9)~=0 & population(:,5)<1960) | (population(:,9)~=0 & population(:,5)>2001));

deaths=d_long(deaths);
% adds lines up to 130

deaths=d_ma0(deaths);
% replaces missings with 0

deaths=d_s1x1(deaths, births);
% split RR deaths counts into Lexis triangles

deaths=d_soainew(deaths);
% split deaths counts in open age intervals into Lexis triangles

deaths=d_unk(deaths);

population=p_unk(population);
% distribute population number of unknown age

population=p_s5x1(population,deaths);
% Splitting 5x1 to 1x1distribute population number of unknown age

population=p_ic(population, deaths, births);
% calculate intercensal estimates for 1990s

population=p_srecm(population,deaths);
% use SR/extinct cohort method

ldb_output(deaths, population, 'mHUN.txt', 'fHUN.txt', births);
% create Lexis files

d_printRA('hun','Hungary');

ldb_psum('mHUN.txt', 'totmales');
ldb_psum('fHUN.txt', 'totfemales');
 