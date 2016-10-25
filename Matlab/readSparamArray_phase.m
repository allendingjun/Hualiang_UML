function [ f, S ] = readSparamArray_phase( filename, p)
fid = fopen(filename);
%     f = [];
%     S = [];
textscan(fid,'%[^-]');
textscan(fid,'%[-]');

No=1;

f = [];
S = [];

while ~feof(fid)
    c = textscan(fid,'%f');
%     f = [];
%     S = [];
%     
    if(No==1)
     f = horzcat(f,c{1,1}(1:2*p:size(c{1,1})));
    end
    
    % S = horzcat(S,c{1,1}(2:3*p:size(c{1,1})) .* exp(1i*c{1,1}(3:3*p:size(c{1,1}))*pi/180));
    
    S = horzcat(S,c{1,1}(2:2*p:size(c{1,1})));
    
    textscan(fid,'%[^-]');
    textscan(fid,'%[-]');
    No=No+1;
end

fclose(fid);

end