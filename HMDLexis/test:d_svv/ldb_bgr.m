addpath('/data/commons/boe/HMD_Lexis_mat2R/MATLAB.hg');
% BGR is test country case, so pick up data from there
addpath('/data/commons/boe/HMDLexis.git/HMDLexis/HMDWORK/BGR/InputDB');
indb_input('BGR');
load BGR;
% load Input DB into Matlab workspace 

deaths=d_long(deaths);
% adds lines up to 130

deaths=d_ma0(deaths);
% replaces missings with 0

deaths=d_s1x1(deaths, births);
% split RR deaths counts into Lexis triangles

% save deaths before and after the test point
save('deaths_pre.mat','deaths');


deaths=d_svv(deaths);
% split some vv deaths counts into Lexis triangles

% save deaths after
save('deaths_post.mat','deaths');

