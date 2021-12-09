function printprogress(curr, total)
    prog=curr/total*100;
    n=round(prog/10);
    str1=repmat('#',1,n);
    str2=repmat('|',1,10-n);
    str=strcat('[',str1,'>',str2,']');
    if curr==1
        fprintf(1,'Progress %3.0f%% %s',prog,str)
    end
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress %3.0f%% %s',prog,str)
    if(curr==total)
        fprintf('\n')
    end
end