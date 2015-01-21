%getscat_lp_lss.m 
%[Sequoia, October 31, 2011]
% Function to get the scattering from LISST-Portable / LISST-StreamSide % data file. Assumes files are in
% LISST-Portable or LISST-StreamSide format with zscat data in columns % 41:80
% This script does NOT use the Zdddhhmm.DAT file, only the Ldddhhmm.DAT
% file produced by the LISST-Portable / LISST-StreamSide
%
% Usage: [scat,tau,zsc,data,cscat,r] = getscat_lp_lss('datafile',readfile,itype)
%
%Input:
% datafile is *EITHER* a Ldddhhmm.dat binary datafile from the Portable % or
% StreamSide *OR* a n x 80 matrix with raw data from a Portable or
% StreamSide
%
% readfile = 1 treats daatfile as a file to be read. All other values
% treats it as a matrix that is read from memory.
%
% iType: 2 for type B, 3 for type C
%
%Output:
% scat is the raw scattering signature in digital counts (i.e. NOT corrected
% with the ringarea file!!)
% tau is the optical transmission
% zsc is the zscat data
% data is the raw data file, converted from binary
% cscat is corrected scattering
% r is the zsc laser pow/laser ref ratio
%
% OAM, 5/27/2009
% OAM, 7/28/2010 
% OAM, 8/22/2011 - removed dcal as input - unnecessary, as it can be
% specified by i_type
% OAM 10/20/2011 - removed conversion to LISST-100X data format

function [scat,tau,zsc,data,cscat,r] = getscat_lp_lss(datafile,readfile,itype)

if nargin == 1
readfile = 1;
itype = 3;
end

if itype == 3
dcal = [1.0000000e+000 1.0038000e+000 9.9360000e-001 1.0027000e+000,...
9.9720000e-001 9.9570000e-001 9.9030000e-001 9.9430000e-001,...
9.9290000e-001 9.9000000e-001 9.9290000e-001 9.9300000e-001,...
9.9150000e-001 9.9300000e-001 9.9230000e-001 9.9090000e-001,...
1.1032000e+000 1.1123000e+000 1.2430000e+000 1.1562000e+000,...
1.3273000e+000 1.1999000e+000 1.0740000e+000 1.7489000e+000,...
1.5382000e+000 2.5109000e+000 2.5468000e+000 3.5504000e+000,...
3.9338000e+000 4.9731000e+000 5.7183000e+000 8.7255382e+000];
else

dcal = [1.0000000e+000 1.0038000e+000 9.9360000e-001 1.0027000e+000,...
9.9720000e-001 9.9570000e-001 9.9030000e-001 9.9430000e-001,...
9.9290000e-001 9.9000000e-001 9.9290000e-001 9.9300000e-001,...
9.9150000e-001 9.9300000e-001 9.9230000e-001 9.9090000e-001,...
1.1032000e+000 1.1123000e+000 1.2430000e+000 1.1562000e+000,...
1.3273000e+000 1.1999000e+000 1.0740000e+000 1.7489000e+000 ,...
1.5382000e+000 2.5109000e+000 2.5468000e+000 3.5504000e+000,...
3.9338000e+000 4.9731000e+000 5.7183000e+000 6.9546000e+000];
end

if readfile == 1
data = tt2mat(datafile,80);%read in the 80-column wide datafile
else
data = datafile;
end

rows = size(data,1);
zsc = data(:,41:80);
r = zsc(:,33)./zsc(:,36);
tau = data(:,33)./r./data(:,36);
scat = ones(rows,32);%pre-allocate scat matrix
cscat = ones(rows,32);

for i = 1:rows
scat(i,:) = (data(i,1:32)/tau(i))-(zsc(i,1:32)*data(i,36)/zsc(i,36));
cscat(i,:) = dcal.*scat(i,:);
cscat(cscat<0)=0;%negative cscats are not possible, so set them to 0.
end

