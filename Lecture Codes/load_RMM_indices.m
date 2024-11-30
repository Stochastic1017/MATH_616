% Loads RMM indices for SSA, SSA-CP

% First day of data to use
year_first=1999;
month_first=1;
day_first=1;

% Last day of data to use
year_last=2013;
month_last=12; %12;
day_last=31; %31;

% Open file
% fid=fopen('/Users/hrogrosky/Documents/Desktop_files/rmm_1974toRealtime.txt','r');
fid=fopen('rmm_1974toRealtime.txt','r');

% Ignore the first two lines of the RMM file. Those are 'header' lines.
fgetl(fid);
fgetl(fid);

% Scan each line until you find the first day to read.
% Note: Last column of text (string) will not be read properly.
formatSpec='%i %i %i %g %g %i %g %s';
sizeLine=[1 8];
linenow=fscanf(fid,formatSpec,sizeLine);
% Read lines until we arrive at year_first
while ~feof(fid) && linenow(1)<year_first
  linenow=fscanf(fid,formatSpec,sizeLine);
end
% Read lines until we arrive at month_first
while ~feof(fid) && linenow(2)<month_first
  linenow=fscanf(fid,formatSpec,sizeLine);
end
% Read lines until we arrive at day_first
while ~feof(fid) && linenow(3)<day_first
  linenow=fscanf(fid,formatSpec,sizeLine);
end
%
% linenow should now contain the data for the first day.
%
% Now start putting together the array of selected data values.
%
% Fill in the first RMM1 and RMM2 values.
% Then add the other values by concatenating (appending to end of array).
%
RMM1=linenow(4);
RMM2=linenow(5);
not_last_day=    linenow(1)<year_last ...
              || linenow(2)<month_last ...
              || linenow(3)<day_last;
while ~feof(fid) && not_last_day
  linenow=fscanf(fid,formatSpec,sizeLine);
  RMM1=[RMM1 linenow(4)];
  RMM2=[RMM2 linenow(5)];
  not_last_day=    linenow(1)<year_last ...
                || linenow(2)<month_last ...
                || linenow(3)<day_last;
end

% Close the data file
fclose(fid);
