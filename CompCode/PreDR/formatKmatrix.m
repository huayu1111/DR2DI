function formatKmatrix(D,file)
%D=[Y K];Y is the label,K is kernel matrix
%file is the output filename(Note it must be a char!)  

[m,n]=size(D);
filename=file;
fid=fopen(filename,'w');
for i=1:m
    fprintf(fid,'%5d',D(i,1));
    fprintf(fid,'%5d: %5d',0,i);
    for j=2:n
        fprintf(fid,'%5d: %f',j-1,D(i,j));
    end
    fprintf(fid,'\n');    
end

fclose(fid);
clear fid i j; 