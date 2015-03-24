%% Profile the code and save HTML output
clear all, close all
profilename = 'profile_case12_origcode';
diary([profilename '.txt'])
profile on
simulate_uncertainty
profile off
p = profile('info')

mkdir(profilename)
profsave(p,profilename)
diary off

% Open up the profile viewer
profview(0,p)