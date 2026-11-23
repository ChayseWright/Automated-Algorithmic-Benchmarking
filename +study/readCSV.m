function [headers,rows] = readCSV(filename)
%READCSV Quoted fields, doubled quotes, CRLF and UTF-8 BOM are supported.
fid=fopen(filename,'rb'); if fid<0, error('tbench:CSV','Cannot read clinical CSV.'); end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8')';
if numel(bytes)>=3&&isequal(bytes(1:3),uint8([239 187 191])), bytes=bytes(4:end); end
text=native2unicode(bytes,'UTF-8');
rows=cell(0,1); row=cell(1,0); field=''; quoted=false; k=1;
while k<=numel(text)
    c=text(k);
    if c=='"'
        if quoted&&k<numel(text)&&text(k+1)=='"', field=[field '"']; k=k+1; %#ok<AGROW>
        else, quoted=~quoted;
        end
    elseif c==','&&~quoted
        row{end+1}=field; field=''; %#ok<AGROW>
    elseif (c==char(10)||c==char(13))&&~quoted
        row{end+1}=field; field='';
        if c==char(13)&&k<numel(text)&&text(k+1)==char(10), k=k+1; end
        if any(~cellfun('isempty',row)), rows{end+1}=row; end %#ok<AGROW>
        row=cell(1,0);
    else, field=[field c]; %#ok<AGROW>
    end
    k=k+1;
end
if quoted, error('tbench:CSV','Unclosed quoted CSV field.'); end
if ~isempty(field)||~isempty(row), row{end+1}=field; rows{end+1}=row; end
if isempty(rows), error('tbench:CSV','Empty clinical manifest.'); end
headers=cellfun(@strtrim,rows{1},'UniformOutput',false); rows=rows(2:end);
if numel(unique(headers))~=numel(headers), error('tbench:CSV','Duplicate column names.'); end
for k=1:numel(rows)
    if numel(rows{k})~=numel(headers), error('tbench:CSV','Wrong column count at row %d.',k+1); end
end
end
