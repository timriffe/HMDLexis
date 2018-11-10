pwd
addpath('/data/commons/boe/HMD_Lexis_mat2R/MATLAB.hg');
% BEL is test country case, so pick up data from there
addpath('/data/wilmoth0/HMD/HMDWORK/BEL/InputDB');

%% run
disp('Read...');

indb_input('BEL');

% load Input DB into Matlab workspace 
load BEL;

%% pre-test calc
%disp('open age interval');
for y=1865:1875
  if ~isempty(find(deaths(:,3)==y & deaths(:,5)==6));
      ind1=find(deaths(:,3)==y & deaths(:,5)==6 & deaths(:,4)==0 & deaths(:,2)==1);
      ind2=find(deaths(:,3)==y & deaths(:,5)==3 & deaths(:,4)==0 & deaths(:,2)==1);
      deaths(ind1,9)=deaths(ind1,9)-deaths(ind2,9);
      deaths(ind1,4)=1;
      deaths(ind1,7)=deaths(ind1,7)-1;
      ind1=find(deaths(:,3)==y & deaths(:,5)==6 & deaths(:,4)==0 & deaths(:,2)==2);
      ind2=find(deaths(:,3)==y & deaths(:,5)==3 & deaths(:,4)==0 & deaths(:,2)==2);
      deaths(ind1,9)=deaths(ind1,9)-deaths(ind2,9);
      deaths(ind1,4)=1;
      deaths(ind1,7)=deaths(ind1,7)-1;
  end
end

d=selif(deaths, deaths(:, end)==1);
d=d(:,1:end-1);
p=selif(population, population(:, end)==1);
p=p(:,1:end-1);
births=selif(births, births(:, end)==1);
births=births(:,1:end-1);

d=selif(d, d(:, 9)>=0);
d=d_unk(d);

d=d_s1x1t(d,births,tadj);
save rsd
d=d_s5x1vv(d);
d=d_s5x1(d);

%% pre
% save deaths before and after the test point
deaths=d;
save('BELdeaths_pre.mat','deaths');

%% split some vv deaths counts into Lexis triangles
d=d_svv(d);
deaths = d;
% save deaths after
save('BELdeaths_post.mat','deaths');


%% done
d=d_s1x1t(d,births,tadj);
d=d_soainew(d,0,1914:1918);
d=d_fill(d);
d=d_long(d);
d=d_ma01(d);
save rsd1

%% population
disp('population');
%p=selif(population, population(:,5)>1923);
p=p_ey2ny(p);
p=p_unk(p);
p1=p_srecmt(p, d,tadj);
p1=selif(p1,p1(:,3)>=85 & (p1(:,5)==1921 | p1(:,5)==1867));
p=selif(p,p(:,3)<85 | (p(:,5)~=1921 & p(:,5)~=1867));
p=[p;p1];
p1=selif(p,p(:,5)>1916);
p=selif(p,p(:,5)<=1916);
save rsp
p=p_postcensal(p,d,births,[1912 1915]);
p=p_precensal(p,d,[1841 1846],1);
p1=p_precensal(p1,d,[1919 1920],1);
p=[p;p1];
p=p_ict(p, d, births,tadj);
p1=p_ict(p1, d, births,tadj);
save rsp1
disp('Extinct cohort method');
p=p_srecmt(p, d,tadj);
save belout
ind=find(p(:,5)>=1915 & p(:,5)<=1918);
p(ind,8)=-1;
ind=find(d(:,3)>=1914 & d(:,3)<=1918);
d(ind,9)=-1;


%% output
ldb_outputt(d, p, 'mBEL.txt', 'fBEL.txt', births, tadj); 
d_printRA('bel','Belgium');

%% move lexis databases
for sexl = {'m' 'f'}
    sexlc = char(sexl);
    fname1 = [hmdpath(Name, 'indb') sexlc Name '.txt'];
    fname2 = [hmdpath(Name, 'lxdb') sexlc Name '.txt'];
    % cmd  = sprintf('xcopy "%s" "%s" /Y', fname1, fname2);
    delete(fname2);  copyfile(fname1, fname2, 'f');  delete(fname1);
end


%% code for user =  Kirill only
if ~isempty(regexpi(getenv('username'), 'Kirill'))
    rmpath(hmdglobalsf('LexisSoftwareFolder'));
    rmpath(hmdpath(Name, 'matlab'));
%     cd('E:\!Home\HMD\HMDWork\DNK\Matlab\');
%     cd('..\indb');
end

% if flgError
%     error('Cannot run script');
% end