function proc(a,b,h)
    load 0/ccx;
    load 0/ccy;
    
    a=str2double(a);
    b=str2double(b);
    h=str2double(h);
    
    cx=a+h/2;
    cy=b+h/2;
    
    ncells=length(ccx);
    
    partition=zeros(ncells,1);
    
    for ic=1:ncells
        if(ccx(ic)<cx && ccy(ic)<cy)
            partition(ic)=0;
        elseif(ccx(ic)>=cx && ccy(ic)<cy)
            partition(ic)=1;
        elseif(ccx(ic)<cx && ccy(ic)>=cy)
            partition(ic)=2;
        else
            partition(ic)=3;
        end
    end
    
    fid=fopen('proc.dat','w');
    fprintf(fid,'%d\n',partition);
    fclose(fid);
end
