for n=2:1:46
    [x, w] = quad_points_legendre(n);
    fname=strcat('leg_', num2str(n));
    fname=strcat(fname, '.txt');
    fileID = fopen(fname,'w');
    fprintf(fileID,'%6.2f %12.8f\n',x, w);
end
fclose(fileID);