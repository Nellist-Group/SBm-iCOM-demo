detectordim=[682 682];
scandim=[13 13];
z=299;

m=zeros(detectordim(1),detectordim(2),scandim(1),scandim(2));

for i=1:scandim(1)
    for j=1:scandim(2)
        filename=strcat('4DSTEM_z=',num2str(z),'_A_pp_',num2str(i),'_',num2str(j),'_Diffraction_pattern_682x682.bin')
        fid =fopen(filename,'r','b');
        %m(:,:,j,i) = fread(fid,[682 682],'uint32');
        m(:,:,j,i) = fread(fid,[682 682],'real*4');

        fclose(fid);

    end       

end

save(strcat('STO_z=',num2str(z),'A_m.mat'),'m')
