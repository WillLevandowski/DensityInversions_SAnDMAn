%% Gravity kernels
i=1;
if have_geometry==0;
   kernel=zeros(length(thickness),length(Lon),length(Lon));
    d=[(0:10:200)'; (220:20:500)'; 1000;2000;3000; max(max(distances))+3000]*1000;
    dyl=ywidth/2;
    dyr=-ywidth/2;
%      topo_weight_8km=E/8000;
%     topo_weight_0km=1-topo_weight_8km;
    for h=1:length(thickness)
        gravity_template=0*d;
        gravity_template_tall=0*d;
        dz1=zs(h)*1000;
        dz2=zs(h+1)*1000;
%         dz1_tall=zs(h)*1000+8000;
%         dz2_tall=zs(h+1)*1000+8000;
        dxl=d+xwidth/2;
        dxr=d-xwidth/2;
        for i=1:length(d)            
            gravity_template(i) = gravprism(1,dxl(i),dxr(i),dyl,dyr,dz1,dz2);
%             gravity_template_tall(i) = gravprism(1,dxl(i),dxr(i),dyl,dyr,dz1_tall,dz2_tall);
        end
        gravity_template(end-1:end)=0;
        k1=interp1(d/1000,gravity_template,distances);
%         k2=interp1(d/1000,gravity_template_tall,distances);
%         k3=k1;
%         for i=1:length(Lon)
%             k3(i)=k1(i,:)*topo_weight_0km+k2(i,:)*topo_weight_8km;
%         end
        kernel(h,:,:)=k1;
        fprintf(['Gravity kernel ' num2str(round(100*h/length(thickness))) ' percent done \n']); 
    end
    % save kernel kernel 
    % dlmwrite('topo_grav', [Lon Lat topo_grav-mean(topo_grav)], '\t')
    %%%  Match to observations
    gridObs=0*Lon; numG=0*Lon;
    F=scatteredInterpolant(ObservedBouguer(:,1),ObservedBouguer(:,2),ObservedBouguer(:,3),'linear','nearest');
    distance=Rd/111e3;
    for j=1:length(Lon)
         if mod(j,500)==0; fprintf([num2str(j) ' stations out of ' num2str(length(Lon)) ' have gravity matched \n']); end

        
            nearby=find((abs(ObservedBouguer(:,2)-Lat(j))+abs(ObservedBouguer(:,1)-Lon(j)))<distance);
           
        gridObs(j)=sum(ObservedBouguer(nearby,3),1)/length(nearby);
        if isempty(nearby); gridObs(j)=F(Lon(j),Lat(j)); end
        numG(j)=length(nearby);
    end
    GravCorr=0.04193*E/1000.*2670; GravCorr=GravCorr-mean(GravCorr);
    gridObs=gridObs-mean(gridObs)+GravCorr;
dlmwrite('ObsGrav', [Lon Lat gridObs-mean(gridObs)], '\t')
dlmwrite('numG', [Lon Lat numG], '\t')
end
%% Things that should be repeated even if we already have the geometry
    initdens_dm=0*initdens;
    for i=1:length(zs)-1
        initdens_dm(:,i)=initdens(:,i)-mean(initdens(:,i));
    end

    grav=0*Lon; 
    for i=1:length(Lon)
        if mod(i,500)==0; fprintf([num2str(i) ' stations out of ' num2str(length(Lon)) ' have gravity predicted \n']); end
        g=zeros(length(Lon),length(thickness));
        for j=1:length(Lon)
                for k=1:length(thickness)
                        g(j,k)=initdens_dm(j,k)*kernel(k,i,j);
                end
        end
    grav(i)=sum(sum(g));  
    end
    PredGrav=grav-mean(grav);
    GravASL=0.04193*E/1000.*surfdens;
    GravASL=GravASL-mean(GravASL);
    PredGrav=PredGrav+GravASL;
    dlmwrite('PredGrav', [Lon Lat PredGrav], '\t')

    resGrav=PredGrav-gridObs;
    resGrav=resGrav-mean(resGrav);
    dlmwrite('resGrav', [Lon Lat resGrav-mean(resGrav)], '\t')
    %%
    surfdens_init=surfdens;
         sed_adj_needed=-resGrav./(0.04193.*(E/1000));
         sed_adj=sed_adj_needed;
        surfdens(ignored)=surfdens(ignored)+sed_adj(ignored);

        foo=find(surfdens>max(surface_density));
        surfdens(foo)=max(surface_density);

        foo=find(surfdens<1800+150/2*E/1000);
        surfdens(foo)=1800+150/2*E(foo)/1000;

        sed_adj=surfdens-surfdens_init;
        dlmwrite('sed_adj', [Lon Lat sed_adj], '\t') 
        dlmwrite('surfdens', [Lon Lat surfdens], '\t') %%%% TopoASL will be added in a second in Setup_MC_Bootstrap
        GravASL=GravASL+0.04193.*E/1000.*sed_adj;
        GravASL=GravASL-mean(GravASL);

        PredGrav=grav-mean(grav);
        PredGrav=PredGrav+GravASL;
        dlmwrite('PredGrav', [Lon Lat PredGrav], '\t')
        resGrav=PredGrav-gridObs;
        resGrav=resGrav-mean(resGrav);
        dlmwrite('resGrav2', [Lon Lat resGrav-mean(resGrav)], '\t')

    
    G=resGrav;
     G=G-mean(G);
