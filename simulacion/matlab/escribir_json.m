function escribir_json(path,value)
%ESCRIBIR_JSON Escribe datos JSON en UTF-8.
    folder=fileparts(path); if ~isempty(folder)&&~isfolder(folder),mkdir(folder);end
    fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'Cannot open %s',path);
    cleaner=onCleanup(@()fclose(fid));
    fwrite(fid,jsonencode(value,'PrettyPrint',true),'char'); fwrite(fid,newline,'char');
end
