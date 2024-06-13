function doimage(handle,dirp,format,name,cl)

filename = [dirp name '.' format];
 if strmatch(format,'epsc','exact')
     filename = filename(1:end-1);
 end
 set(handle, 'PaperPositionMode', 'auto');
%   print(handle, '-r0', [dirp name '.' format], ['-d' format]);
  print(handle, filename, ['-d' format]);
 if cl
    close(handle)
 end