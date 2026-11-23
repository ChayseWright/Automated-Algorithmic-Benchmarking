function appendCSV(filename,rows)
%APPENDCSV Stream large explanatory tables without holding a cohort in memory.
if isempty(rows), return; end
folder=fileparts(filename); if ~exist(folder,'dir'), mkdir(folder); end
if ~exist(filename,'file'), tbench.writeCSV(filename,rows); return; end
fid=fopen(filename,'r','n','UTF-8'); if fid<0, error('tbench:Write','Cannot inspect %s.',filename); end
header=fgetl(fid); fclose(fid);
if ~strcmp(header,strjoin(fieldnames(rows)',',')), error('tbench:CSVSchema','Refusing to append mismatched CSV columns.'); end
fid=fopen(filename,'a','n','UTF-8'); if fid<0, error('tbench:Write','Cannot append %s.',filename); end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fields=fieldnames(rows);
for k=1:numel(rows)
    cells=cell(1,numel(fields));
    for j=1:numel(fields)
        value=rows(k).(fields{j});
        if ischar(value), cells{j}=['"' strrep(value,'"','""') '"'];
        elseif (isnumeric(value)||islogical(value))&&isscalar(value), cells{j}=sprintf('%.15g',double(value));
        else, error('tbench:CSV','Streamed CSV needs scalar/text fields: %s.',fields{j});
        end
    end
    fprintf(fid,'%s\n',strjoin(cells,','));
end
end
