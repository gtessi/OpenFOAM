function fixFloat(a,b,h,lx,ly,nx,ny)
    a=str2double(a);
    b=str2double(b);
    h=str2double(h);
    lx=str2double(lx);
    ly=str2double(ly);
    nx=str2double(nx);
    ny=str2double(ny);
    
    dx=lx/nx;
    dy=ly/ny;
    
    if(((a/dx)-floor(a/dx))>=0.5)
        a=floor(a/dx)*dx+dx;
    else
        a=floor(a/dx)*dx;
    end
    
    if(((b/dy)-floor(b/dy))>=0.5)
        b=floor(b/dy)*dy+dy;
    else
        b=floor(b/dy)*dy;
    end
    
    if(dx<=dy)
%         if(((h/dx)-floor(h/dx))>=0.5)
%             h=floor(h/dx)*dx+dx;
%         else
            h=floor(h/dx)*dx;
%         end
    else
%         if(((h/dy)-floor(h/dy))>=0.5)
%             h=floor(h/dy)*dy+dy;
%         else
            h=floor(h/dy)*dy;
%         end
    end
    
    fid=fopen('fixedValues.dat','w');
    fprintf(fid,'%2.10f\n',a);
    fprintf(fid,'%2.10f\n',b);
    fprintf(fid,'%2.10f\n',h);
    fclose(fid);
end
