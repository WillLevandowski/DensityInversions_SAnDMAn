%%%%% *S*imulated *An*nealing *D*ensity *M*odeling *A*lgorithm, *N*eat!
%%%%% A simulated annealing algorithm to "invert" (really a forward model)
%%%%% gravity and topography for lithospheric density structure. Please
%%%%% cite Levandowski et al., G-cubed 2015. Direct any questions or
%%%%% concerns to wlevandowski@usgs.gov; wlevando@gmail.com. 
cc

%% Step 0: Initialize 
%%% Physical parameters
Default_ElasticThickness=30; %%% default value of elastic thickness, in km
rhoa=3200; %%% asthenospheric density
H0=-2.4; %%% Lachenbruch and Morgan, 1990
Youngs=70e9; %%%% Young's modulus
Poissons=0.25; %%%% Poisson's ratio

%%%~~~~~~~~~~~~Model Geometry and Spatial Extent~~~~~~~~~~~~%%%
%%% The lithosphere is divided into a number of layers that begin at
%%% user-defined (and flexible) depths, "zs" below sea-level. As written,
%%% the code isolates material above sea-level, assuming that its density
%%% is already known accurately. Feel free to change or contact Will for
%%% help therewith.
zs=[0 5 10 15 20 30 40 50 60 70 85 100 120 140 160 180 200]';
MaxDepth=max(zs);
thickness=diff(zs);
NumLayers=length(thickness);
%%% The model region, in Longitude/Latitude
left=75;
right=93;
top=45;
bottom=35.5;

to_overwrite=input('Do you want to overwrite existing files? \n    0 to append, 1 to overwrite: \n');
ObservedElevation=load('tarim_elevation.txt');
ObservedBouguer=load('tarim_gravity.txt');
moho=load('tarim_moho.txt');
load initial_positions; 
load initial_density; 
load initial_surfdens; 

Longitude_input=initial_positions(:,1);
Latitude_input=initial_positions(:,2);
Density_input=initial_density;
surface_density=initial_surfdens(:,3);
have_geometry=0;    %%% a logical that allows repeated simulations with the 
                    %%% same geometry. A BIG time-saver!!
new_starting_model=0;   %%% if == 0; use the same starting model. 
                        %%% if == 1; use a different starting model each
                        %%% simulation. Rare, but recommended (c.f.,
                        %%% Levandowski et al., 2015 G-cubed and
                        %%% 2017 Nature Communications)
variable_elastic_thickness=1;
if variable_elastic_thickness == 0
    Te=Default_ElasticThickness; 
end

%%%%%%%%% User-defined parameters for the simulated annealing
%%%%%%%%% Tolerated misfits
tol=2.5; %%% User-defined: half of the maximum allowable gravity residual (mGal)
        %%% it is assumed that the maximum allowable topography is ten
        %%% times as many meters. For example; tol=5 means max. grav.
        %%% residual = 10 mGal and max. topo. residual = 100 meters. 
maxGrav=tol*2;
maxTopo=tol/50;
outliers=0; %%% the number of nodes that can keep gravity or topography 
            %%% residuals above tolerance
num_reals=5; %%% how many realizations to run with a given input model 
             %%% and geometry            
num_sims=5 ; %%% the number of geometries to run through: 
            %%% The total size of posterior distribution of 3D density 
            %%% will therefore be num_reals*num_sims
%%%%%%%%% Allowable density perturbations each step and overall
% max_crust=250; %%% twice the maximum change to crustal density in one step
max_crust=200; %%% twice the maximum change to crustal density in one step
% max_mantle=40;
max_mantle=50;
max_drho_crust=200; %%% the overall maximum allowable change to the starting model
max_drho_mantle=50; 
%%% For example, max_crust=200 & max_drho_crust=150 mean that each time the
%%% algorithm hits on a crustal cell, it will try perturbations to the
%%% current model of up to +/- 100 kg/m^3, but the allowable range for a
%%% cell with initial density of 2700 kg/m^3 is limited to 2550-2850.


%% Step 1: Set up problem's geometry

spacing_average=40;
for simulations=1:num_sims
    tic
    if have_geometry == 0
        spacing=spacing_average-10*rand; %%% cell width, in km 
        xright=(right-left)*cosd(bottom)*111;
        xleft=-(rand)/2*spacing;
        nx=ceil((xright-xleft)/spacing);
        ytop=(top-bottom)*111+spacing*(rand)/2;
        ybottom=-spacing*rand/2;
        ny=ceil((ytop-ybottom)/spacing);

        X=linspace(xleft,xright,nx);
        Y=linspace(ybottom,ytop,ny);

        [xg, yg]=meshgrid(X,Y);
        X=reshape(xg,nx*ny,1);
        Y=reshape(yg,nx*ny,1);
        xwidth=max(diff(X))*1000;
        ywidth=mode(diff(Y))*1000;
        Rd=sqrt(xwidth*ywidth/pi); %%% the radius of a circle with same area

        Lat=Y/111+bottom;
        Lon=left+X./(111*cosd(Lat));
        okay=find(Lon<=right);  X=X(okay); Y=Y(okay); Lon=Lon(okay); Lat=Lat(okay);
        ignored=find(Lon<73 | Lon>95 | Lat<35.5 | Lat>44.25); %%% a chance to say that you don't care about reproducing
                    %%% grav & topo in some regions but still want to consider
                    %%% their density when reconciling the rest of the region.
                    ignored=[];
        ignore=0*Lon; ignore(ignored)=1; notignored=find(ignore==0);
        NumNodes=length(Lon);

         distances=zeros(NumNodes);
        for i=1:NumNodes
            for j=1:NumNodes
                latdiff=Y(i)-Y(j);
                longdiff=X(i)-X(j);
                distances(i,j)=sqrt(latdiff^2+longdiff^2);
            end
        end
        E=0*Lon;
        M=0*Lon;
        %%%% Find the average elevation in each box
        for i=1:NumNodes
            nearby=[]; distance=0.5*(xwidth+ywidth)/111e3; %% convert meters width to degrees
            while isempty(nearby)
                nearby=find(abs(ObservedElevation(:,1)-Lon(i))<distance & abs(ObservedElevation(:,2)-Lat(i))<distance);
                distance=distance+0.1; %%% if nothing has been found, 
                                       %%% it will search a larger radius next time
            end
            E(i)=mean(ObservedElevation(nearby,3),1);
        end
            %%%% Find the average moho in each box, added by Yangfan
            for i=1:NumNodes
            nearby=[]; distance=0.5*(xwidth+ywidth)/111e3; %% convert meters width to degrees
            while isempty(nearby)
                nearby=find(abs(moho(:,1)-Lon(i))<distance & abs(moho(:,2)-Lat(i))<distance);
                distance=distance+0.1; %%% if nothing has been found, 
                                       %%% it will search a larger radius next time
            end
                M(i)=mean(moho(nearby,3),1);
            end

    end
    %% Step 2: Interpolate the density model
    if have_geometry == 0 || new_starting_model == 1
        %%%%% and interpolate the density model to these points
        initdens=zeros(NumNodes,NumLayers);
        for i=1:NumLayers
            F=scatteredInterpolant(Longitude_input, Latitude_input, Density_input(:,i),'linear','nearest');
            initdens(:,i)=F(Lon,Lat);
        end
        F=scatteredInterpolant(Longitude_input, Latitude_input,surface_density,'linear', 'nearest');
        surfdens=F(Lon,Lat); dlmwrite('surfdens', [Lon Lat surfdens], '\t')
    end

    %% Step 3: Calculate topographic Green's functions and smooth 
    if have_geometry==0
        if variable_elastic_thickness == 1
            load Elastic_Thickness %%% Lon Lat Elastic_Thickness(km)
            okay=find(Elastic_Thickness(:,3)>=25);
            Elastic_Thickness=Elastic_Thickness(okay,:);
            F=scatteredInterpolant(Elastic_Thickness(:,1),Elastic_Thickness(:,2),Elastic_Thickness(:,3),'linear', 'nearest');
            Te=F(Lon,Lat);%%% 2D elastic thickness
            D=Youngs*((Te*1000).^3)/(12*(1-Poissons^2)); %%% 2D flexural rigidity
			flex_param=(D/((rhoa-2000*rand)*9.8)).^0.25/1000; %%% 2D flexural parameter
               %%% To first order, the flexural influence of a load beneath point i on the surface at point j
               %%% can be approximated by using the distance between i&j and the minimum strength 
               %%% between them. 
               %%% That is, even if a piece of lithosphere is strong, its loads are not transmitted efficiently
               %%% across very weak zones           
		   beta=zeros(NumNodes); 
            for i=1:NumNodes
                l=NumNodes;
                if mod(i,500)==0; fprintf(['Calculating distances ' num2str(round(100*i/l)) ' percent done \n']); end
                for j=1:NumNodes
                    beta(i,j)=min([flex_param(i) flex_param(j)]);
                end
            end
        else
            D=Youngs*((Te*1000).^3)/(12*(1-Poissons^2));
            beta=(D/((rhoa-2000*rand)*9.8)).^0.25/1000;
        end 

        coeff=distances./beta;


        load bessel_table
        beta=1000*beta; norm_bessel=0*coeff; 
        ber=interp1(berx(:,1),berx(:,2),coeff);
        bei=interp1(beix(:,1),beix(:,2),coeff);
        ker=interp1(kerx(:,1),kerx(:,2),coeff);
        kei=interp1(keix(:,1),keix(:,2),coeff);
        ker_prime=interp1(kerpB(:,1),kerpB(:,2),Rd./beta);
        ber_prime=interp1(berpB(:,1),berpB(:,2),Rd./beta);
        kei_prime=interp1(keipB(:,1),keipB(:,2),Rd./beta);
        bei_prime=interp1(beipB(:,1),beipB(:,2),Rd./beta);

        bw=((Rd./beta).*ber_prime.*ker - (Rd./beta).*bei_prime.*kei);
        %%%%% Watts, 2001 Cambridge, eqns. 3.54 & 3.55
        for i=1:length(Lon) 
            if variable_elastic_thickness == 1
                bw(i,i)=((Rd./beta(i,i)).*ker_prime(i,i).*ber(i,i)) - (Rd./beta(i,i)).*kei_prime(i,i).*bei(i,i) +1;
            else
                bw(i,i)=((Rd./beta).*ker_prime.*ber(i,i)) - (Rd./beta).*kei_prime.*bei(i,i) +1;
            end
        end

        for i=1:length(Lon); norm_bessel(i,:)=bw(i,:)/sum(bw(i,:)); end

        for i=1:length(Lon)
            if sum(bw(:,i))>1; bw(:,i)=bw(:,i)/sum(bw(:,i)); end
            if sum(bw(i,:))>1; bw(i,:)=bw(i,:)/sum(bw(i,:)); end
            if sum(bw(:,i))<0.25; bw(:,i)=0.25*bw(:,i)/sum(bw(:,i)); end
        end

        besselweights=bw;
        clear ber bei ker kei bei_prime ber_prime ker_prime kei_prime coeff beta bessel_table
    end

     %% Step 4: Calculate Airy elevation
     H_BSL=0*Lon;
    for i=1:NumNodes
        H_BSL(i)=(rhoa-initdens(i,:))/rhoa*thickness;
    end
    H_ASL=E/1000.*(rhoa-surfdens)/rhoa;
    H=H_BSL+H_ASL;
    PredTopo=H+H0;
    RoughRestopo=PredTopo-E/1000; 
    
    %%%% Adjust the density so that the average residual topography across
    %%%% the region is 0. 
    initdens=initdens+mean(RoughRestopo)*rhoa/max(zs);
    H_BSL=0*Lon;
    for i=1:NumNodes
        H_BSL(i)=(rhoa-initdens(i,:))/rhoa*thickness;
    end
    H=H_BSL+H_ASL;
    PredTopo=H+H0;
     
         RoughRestopo=PredTopo-E/1000; 
        dlmwrite('RoughRestopo1', [Lon Lat RoughRestopo], '\t')
    %%%% see variable *ignored* defined on line 110; here, we will set the
    %%%% density of this region such that the Airy residual topography is 0.
    %%%% This step is helpful in diminshing edge effects without
    %%%% introducing the computational cost normally associated with
    %%%% expanding the model to poorly constrained regions
     
        
       %% Step 5: Gravity in a separate script
    GravCalc_Final
    H_BSL=0*Lon;
    for i=1:NumNodes
        H_BSL(i)=(rhoa-initdens(i,:))/rhoa*thickness;
    end
        H_ASL=E/1000.*(rhoa-surfdens)/rhoa;
     H=H_BSL+H_ASL;
    PredTopo=H+H0;
    RoughRestopo=PredTopo-E/1000; 
   
    
     %%%%% Moral of a cautionary tale: 
    %%%%% If the average RoughRestopo is not near zero, something is
    %%%%% systematically wrong with the input density model. I STRONGLY suggest
    %%%%% trying to figure that issue out before feeding the density model to
    %%%%% SAnDMAN; the cludges above are meant as an insurance policy. 
            dlmwrite('RoughRestopo', [Lon Lat RoughRestopo], '\t')
         SmoothRestopo=0*Lon;
                fprintf('Calculating residual topography \n')
        for i=1:length(Lon)
            SmoothRestopo(i)=besselweights(i,:)*RoughRestopo;
        end
        Restopo=SmoothRestopo-mean(SmoothRestopo);
        
        RawGravResid=G;
        RawTopoResid=Restopo;
        dlmwrite('RawRestopoGrid', [Lon Lat Restopo], '\t')
        dlmwrite('RawResGravGrid', [Lon Lat G], '\t')

        
        
    %% Step 6: Run the darn thing!
    adjs=zeros(NumNodes,NumLayers+1,num_reals);
    for realizations=1:num_reals
        Current_Topo_Resid=RawTopoResid; 
        Current_Grav_Resid=RawGravResid; 

        adj=zeros(NumNodes,NumLayers+1); count=0; trials=0; mark=zeros(length(Lon),1);

        CTR_care=Current_Topo_Resid; 
         CTR_care(ignored)=0;
        CGR_care=Current_Grav_Resid; 
         CGR_care(ignored)=0;


        while  length(find(abs(CGR_care)>maxGrav))+length(find(abs(CTR_care)>maxTopo))>outliers && count<100000
            count=count+1; rand1=rand; 
            %%%%% The likelihood of a node's selection depends on its current
            %%%%% residuals...
            search1=rand1*max(abs(Current_Topo_Resid));
            search2=rand1*max(abs(Current_Grav_Resid));
            wrong=[find( abs(Current_Topo_Resid)>search1 );find(abs(Current_Grav_Resid)>search2)];
            %%%%% And the lucky winner is...
            beneath=ceil(rand*length(wrong)); beneath=wrong(beneath);

                %%% output the status to stdout
                if mod(count,5000)==0  %%% change frequency as desired
                    fprintf([num2str(simulations) ' Simulations:  \n' ]);
                     fprintf([num2str(realizations) ' Realizations, ' num2str(count) ' attempts:  \n']);
                    fprintf([num2str(round(spacing)) ' km spacing.  \n']);
                    fprintf([ '      ' num2str(  [  num2str(length(    find(   abs(CGR_care)>tol*2)))  ' Grav wrong. As bad as ' num2str(max(abs(CGR_care)))]) ' mGal \n'])
                    fprintf([ '      ' num2str(  [ num2str(length(    find(   abs(CTR_care)>maxTopo)))  ' Topo wrong. As bad as ' num2str(max(abs(CTR_care))*1000)]) ' meters \n '])
                    fprintf([ '        '   num2str(length(    find(   abs(CGR_care)>tol)))  ' Grav kinda wrong. \n'])
                    fprintf([ '      '   num2str(mean(abs(CGR_care)))  ' mGal on average  \n '])
                    fprintf([ '         '   num2str(length(    find(   abs(CTR_care)>tol/100)))  ' Topo kinda wrong. \n '])
                    fprintf([ '      '   num2str(mean(abs(CTR_care)*1000))  ' meters on average  \n '])
                    dlmwrite('CurrentGravResid', [Lon Lat Current_Grav_Resid], '\t')
                    dlmwrite('CurrentTopoResid', [Lon Lat Current_Topo_Resid], '\t')
                end


                %%%% The node has been selected, now pick which cells in that
                %%%% column to perturb and how much to perturb them. Then, decide
                %%%% which cell-perturbation pair is best (even if it is worse
                %%%% than the current density model).
                 gravweight=35*(length(find(abs(CTR_care)>maxTopo))+1)/(length(find(abs(CGR_care)>maxGrav))+1);
                 %%%% the weighting factor to reconcile the different orders of
                 %%%% magnitude of gravity and topography. Making the weight
                 %%%% adaptive is really a big help. Feel free to change that
                 %%%% coefficient out front. Ideally, the number of gravity
                 %%%% residuals > maxGrav should roughly equal the number of topo
                 %%%% residuals > maxTopo throughout the simulated annealing.
                 mark(beneath)=mark(beneath)+1;
                 attempts=min([2+floor(mark(beneath)/5) length(thickness)+1]);
                 %%%% the number of cell-density perturbation pairs to try
                 layer=zeros(attempts,1);
                 drho=layer;
                 var_Topo=layer;
                 var_Grav=layer;
                 Topo_Check=zeros(length(Lon),attempts);
                  Gravity_Check=zeros(length(Lon),attempts);

              for attempt=1:attempts
                  %%%% which cell beneath the node selected
                 layer(attempt)=ceil(rand*(NumLayers+1));
                   %%%% assign it a density perturbation relative to current model
                    drho(attempt)=max_crust*(rand-0.5);
                    if  drho(attempt)+adj(beneath,layer(attempt))>max_drho_crust 
                        drho(attempt)=-rand;
                    end
                    if  adj(beneath,layer(attempt))>max_drho_crust 
                        drho(attempt)=-rand;
                    end
                    if drho(attempt)+adj(beneath,layer(attempt))<-max_drho_crust 
                        drho(attempt)=rand;
                    end
                    if adj(beneath,layer(attempt))<-max_drho_crust 
                        drho(attempt)=rand;
                    end
                    if zs(layer(attempt))>=M(beneath) && layer(attempt)<length(thickness)+1 
                        drho(attempt)=max_mantle*(rand-0.5);
                        if  (drho(attempt)+adj(beneath,layer(attempt)))>max_drho_mantle;
                            drho(attempt)=-rand;
                        end
                        if  adj(beneath,layer(attempt))>max_drho_mantle;
                            drho(attempt)=-rand;
                        end
                        if drho(attempt)+adj(beneath,layer(attempt))<-max_drho_mantle;
                            drho(attempt)=rand;
                        end
                        if adj(beneath,layer(attempt))<-max_drho_mantle;
                            drho(attempt)=rand;
                        end
                    end
              end

            %%%% Calculate the attendant gravity and flexural topography changes
                for density=1:length(drho)
                    if layer(density)==length(thickness)+1
                        topo_changes=besselweights(:,beneath)*E(beneath)/1000*drho(density)/rhoa;
                        grav_changes=0*Lon;
                        grav_changes(beneath)=0.04193*E(beneath)/1000*drho(density);
                    else
                    topo_changes=besselweights(:,beneath)*thickness(layer(density))*drho(density)/rhoa;
                    grav_changes=kernel(layer(density),:,beneath)' *drho(density);
                    end
                    Topo_Check(:,density)=Current_Topo_Resid-topo_changes; 
                    TC1=Topo_Check(:,density); 
                    %%%% this is a little trick: I don't really care if a node has
                    %%%% a gravity residual of 0.5 mGal or 2 mGal. I only care
                    %%%% about those where the misfit is rather large. So if the
                    %%%% misfit is less than "tol" for gravity or "tol"/100 for
                    %%%% topography, the residual does not factor into the decision
                    %%%% about which cell-density perturbation pair is "best"
                    good=find(TC1<0 & TC1>-tol/100); TC1(good)=0;
                    good=find(TC1>0 & TC1<tol/100); TC1(good)=0;
                    rest=find(TC1<0 ); TC1(rest)=TC1(rest)+tol/100;
                    rest=find(TC1>0 ); TC1(rest)=TC1(rest)-tol/100;
                    var_Topo(density)=var(TC1);

                    Gravity_Check(:,density)=Current_Grav_Resid+grav_changes; 
                    GC1=Gravity_Check(:,density); 
                    good=find(GC1<0 & GC1>-tol); GC1(good)=0;
                    good=find(GC1>0 & GC1<tol); GC1(good)=0;
                    rest=find(GC1<0 ); GC1(rest)=GC1(rest)+tol;
                    rest=find(GC1>0 ); GC1(rest)=GC1(rest)-tol;
                    var_Grav(density)=var(GC1/gravweight);
                end
               %%%% Ladies and gentlemen: I present to you the highly disputed
               %%%% oddly weighted champion of the iteration!!!!
                 
               best=find((var_Grav+100).*(var_Topo+100)==min((var_Grav+100).*(var_Topo+100)) ,1,'first' );
                 %%%% those +100s up there regularize/stabilize things 

                adj(beneath,layer(best))=adj(beneath,layer(best))+drho(best);
                Current_Topo_Resid=Topo_Check(:,best);
                Current_Grav_Resid=Gravity_Check(:,best);
                CTR_care=Current_Topo_Resid; 
                CTR_care(ignored)=0;
                CGR_care=Current_Grav_Resid; 
                 CGR_care(ignored)=0;
        end
    adjs(:,:,realizations)=adj;
    end
    save adjs.mat adjs
    clear adj;
    adj=mean(adjs,3);
    have_geometry=0;
    adj_all=[]; positions_all=[]; initdens_all=[]; surfdens_all=[];
    if simulations>to_overwrite
        load adj_all
        load positions_all
        load initdens_all
        load surfdens_all
    end
    adj_all=[adj_all; adj];
    positions_all=[positions_all; [Lon Lat]];
    initdens_all=[initdens_all; initdens];
    surfdens_all=[surfdens_all; surfdens];
    dlmwrite('adj_all', adj_all, '\t')
    dlmwrite('positions_all', positions_all, '\t')
    dlmwrite('initdens_all', initdens_all, '\t')
    dlmwrite('surfdens_all', surfdens_all, '\t')
    toc


surfdens_all=surfdens_all+adj_all(:,end);
adj_all=adj_all(:,1:end-1);
finaldens=initdens_all+adj_all;
finaldens_all=finaldens;

finalspace=0.3;
Lons=Lon;
Lats=Lat; 
inits=[]; pos=[]; found=[]; finals=[];surfdensity=[];
for i=1:length(Lons)
        distance=finalspace/2+0.01;
        nearby=find(abs(positions_all(:,1)-Lons(i) )*cosd(Lats(i))<=distance & abs(positions_all(:,2)-Lats(i))<=distance);
        finals(end+1,:)=mean(finaldens(nearby,:),1);
        pos(end+1,:)=[Lons(i) Lats(i)];
        surfdensity(end+1,1)=mean(surfdens_all(nearby));
        found(end+1,1)=length(nearby);
        inits(end+1,:)=mean(initdens_all(nearby,:),1);
end
okay=find(found ); 
inits=inits(okay,:);
pos=pos(okay,:);
finaldens=finals(okay,:);
surfdensity=surfdensity(okay);
adjs=finaldens-inits;
% 

for i=1:NumLayers
    rhofinal=strcat('rho',num2str(i),'_final');
    dlmwrite(rhofinal, [pos finaldens(:,i)], '\t');
    rhoadj=strcat('rho',num2str(i),'_adj');
    dlmwrite(rhoadj, [pos adjs(:,i)], '\t');
    rhoinit=strcat('rho',num2str(i),'_init');
    dlmwrite(rhoinit, [pos inits(:,i)], '\t')
end

dlmwrite('rho0_final', [pos surfdensity], '\t')

dlmwrite('final_positions', pos, '\t')
dlmwrite('final_density', finaldens, '\t')
dlmwrite('final_surfacedensity', surfdensity, '\t')

dlmwrite('rho_0_15', [pos (finaldens(:,1)+finaldens(:,2)+finaldens(:,3))/3], '\t')
dlmwrite('rho_15_30', [pos finaldens(:,4)/3+finaldens(:,5)*2/3], '\t')
dlmwrite('rho_30_50', [pos finaldens(:,6)/2+finaldens(:,7)/2], '\t')
dlmwrite('rho_50_70', [pos finaldens(:,8)/2+finaldens(:,9)/2], '\t')
dlmwrite('rho_70_120', [pos finaldens(:,10)*15/50+finaldens(:,11)*15/50+finaldens(:,12)*20/50], '\t')
dlmwrite('rho_120_200', [pos finaldens(:,13)/4+finaldens(:,14)/4+finaldens(:,15)/4+finaldens(:,16)/4], '\t')

end