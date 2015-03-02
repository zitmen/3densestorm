classdef IO
    
    methods (Static)
        
        function [info,frames,width,height] = getImageInfo(filepath,cal)
            info = imfinfo(filepath);
            frames = 1:numel(info);
            width  = info(1).Width;
            height = info(1).Height;
            if cal.is_biplane
                if cal.divide_dim == 1
                    height = height / 2;
                else
                    width = width / 2;
                end
            end
        end
        
        function I = readImage(file,cal,frame)
            im = imread(file.path,frame,'Info',file.info);
            if cal.is_biplane
                I(:,:,1) = im(1:file.height,1:file.width);
                if cal.divide_dim == 1
                    I(:,:,2) = im(file.height+1:end,1:file.width);
                else
                    I(:,:,2) = im(1:file.height,file.width+1:end);
                end
            else
                I = im;
            end
        end
        
        function writeHeader(fpath,header)
            hstr = header{1};
            for i = 2:length(header), hstr = [hstr,',',header{i}]; end;
            fid = fopen(fpath,'w');
            fprintf(fid,'%s\r\n',hstr);
            fclose(fid);
        end
        
        function appendResults(fpath,data)
            dlmwrite(fpath,data,'-append','delimiter',',');
        end
        
        function sortResults(fpath,header)
            try
                data = dlmread(fpath,',',1,0);
                data = sortrows(data,[1,2,3,4,5,6]);
                IO.writeHeader(fpath,header);
                IO.appendResults(fpath,data);
            catch
                % this usually happens in case when the file is empty
            end
        end
        
    end
    
end
