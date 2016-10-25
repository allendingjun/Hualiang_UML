function [ F57, S57 ] = readSparamArray( filename, p ,freq)
%%% Obtain both phase and magnitude

% % % % % % example for using
% % % % % % 
% % % % % % filename1=['./data/SingleLayer_P120_1.25THz'];
% % % % % % [F, S]=readSparamArray([filename1 '.txt'],1,1.25);
% % % % % % Mag_S1=abs(S);
% % % % % % Deg_S=angle(S)*180/pi;
% % % % % % 
% % % % % % Deg_S1=wrapTo360(Deg_S);


fid = fopen(filename);
%     f = [];
%     S = [];
textscan(fid,'%[^-]');
textscan(fid,'%[-]');
No=1;
F57=[];
S57=[];
while ~feof(fid)
    c = textscan(fid,'%f');
    f = [];
    S = [];
    f = horzcat(f,c{1,1}(1:3*p:size(c{1,1})));
    S = horzcat(S,c{1,1}(2:3*p:size(c{1,1})) .* exp(1i*c{1,1}(3:3*p:size(c{1,1}))*pi/180));
    
    Idx=find(f==freq);
    F57=[F57 f(Idx)];
    S57=[S57 S(Idx)];
    
    No=No+1;
    
    textscan(fid,'%[^-]');
    textscan(fid,'%[-]');
end
fclose(fid);
end

