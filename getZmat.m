function res = getZmat(molDir)
warning('off')
setenv('GAUSS_EXEDIR', 'C:/G09W/');
if ~exist([molDir,'.out'],'file')
    status = system(['C:/G09W/g09.exe ',molDir,'.gjf '...
        molDir,'.out']);
    if status == 1
        disp('something went wrong but cool, since this happens sometimes, just rerunning it')
        status = system(['C:/G09W/g09.exe ',molDir,'.gjf '...
            ,molDir,'.out']);
        disp(status);
    elseif status == 0
        disp('gaussian calc done')
    else
        disp(['something is wrong, status: ',num2str(status)]);
    end
else
    disp('found existing output file')
end
res = getG09geoms([molDir,'.out'],'last');

end