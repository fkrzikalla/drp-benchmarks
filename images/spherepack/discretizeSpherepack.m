% This is a pack of 621 grains with an initial box size of 102^3. The GD
% simulation compacts the pack to specific stress, using servo-control 
% boundary conditions. 

% Author: Fabian Krzikalla, Stanford University
% Date: July 1, 2012

simpath='./';
modpath='./';

% grain positions, the number 1 refers to the reference state
positions=load(strcat(simpath,'final_loc_cn1.dat'));

% lenght of the box reference state
boxdim=[55.185506 55.404375 55.540488]; % for quartz 10MPa

% Material properties for computing forces from displacements / overlaps
Gg=44;
PR=0.074194;
Kg=Gg*(2+2*PR)/(3-6*PR);
Eg=2*Gg*(1+PR);
Eeff=Eg/(2*(1-PR^2));

N=size(positions,1);      % 621 spheres
rad=3.5;                    % grain radius
V=N*4/3*pi*rad^3;           % volume not considering overlap

L=geomean(boxdim); Vbox=prod(boxdim);
porosity=(Vbox-V)/Vbox;   % undeformed initial porosity

%% Analyze distances between grains

% FIRST, duplicate grains close to boundaries
N2=N;

% identify grains to be duplicated
[iFrt]=find(positions(:,1)<2*rad);          % front
[iBck]=find(positions(:,1)>(boxdim(1)-2*rad));  % back
% compute index ranges for duplicate grains
nFrt=size(iFrt,1); rFrt=(N2+1):(N2+nFrt); N2=N2+nFrt;
nBck=size(iBck,1); rBck=(N2+1):(N2+nBck); N2=N2+nBck;
% copy grains to the corresponding opposite side
positions(rFrt,:)=positions(iFrt,:);
positions(rFrt,1)=positions(iFrt,1)+boxdim(1);
positions(rBck,:)=positions(iBck,:);
positions(rBck,1)=positions(iBck,1)-boxdim(1);

% now do the same for right and left boundaries in the 2-direction
[iRig]=find(positions(:,2)<2*rad);          % right
[iLef]=find(positions(:,2)>(boxdim(2)-2*rad));  % left
nRig=size(iRig,1); rRig=(N2+1):(N2+nRig); N2=N2+nRig;
nLef=size(iLef,1); rLef=(N2+1):(N2+nLef); N2=N2+nLef;
positions(rRig,:)=positions(iRig,:);
positions(rRig,2)=positions(iRig,2)+boxdim(2);
positions(rLef,:)=positions(iLef,:);
positions(rLef,2)=positions(iLef,2)-boxdim(2);

% and the top and bottom boundaries in the 3-direction
[iTop]=find(positions(:,3)<2*rad);          % top
[iBot]=find(positions(:,3)>(boxdim(3)-2*rad));  % bottom
nTop=size(iTop,1); rTop=(N2+1):(N2+nTop); N2=N2+nTop;
nBot=size(iBot,1); rBot=(N2+1):(N2+nBot); N2=N2+nBot;
positions(rTop,:)=positions(iTop,:);
positions(rTop,3)=positions(iTop,3)+boxdim(3);
positions(rBot,:)=positions(iBot,:);
positions(rBot,3)=positions(iBot,3)-boxdim(3);

% SECOND, compute 3-component distance matrix
distance=zeros(N2,N2,3);
for i=1:N2
for j=1:3
    distance(:,i,j)=positions(:,j)-positions(i,j);
end
end

distanceabs=sqrt(sum(distance.^2,3));       % absolute distances
conts=((distanceabs<rad*2)&(distanceabs~=0)); % contact matrix
conum=sum(conts,2);                         % coordination numbers

%% Analyse contact forces 
% Now that we know the contact points and distances, we can find the
% contact points and associated contact radii.

% find the list of contact grain pairs and compute the contact radius
[grain1,grain2,overlap]=find(triu(conts.*(2*rad-distanceabs)));
radContact = sqrt(rad*overlap/2);   % contact radii
posContact = (positions(grain1,1:3)+positions(grain2,1:3))/2; % positions
dirContact = (positions(grain1,1:3)-positions(grain2,1:3))/2; % directions
numContact = size(posContact,1);  % total number of contact points

kn     = 4/3*Eeff*sqrt(rad);  % contact stiffness
forcen = kn*overlap.^(3/2); % normal force arising from contact overlap
pressn = forcen./(pi.*radContact.^2); % average contact pressure

%% Discretize geometry 

% Generate grid
dx=rad/50;
xrange=dx/2:dx:(boxdim(1)-dx/2); 
yrange=dx/2:dx:(boxdim(2)-dx/2); 
zrange=dx/2:dx:(boxdim(3)-dx/2); 
nx=size(xrange,2);
ny=size(yrange,2);
nz=size(zrange,2);

% Start with empty model
model=zeros(nx,ny,nz);

% Populate model with digital spheres
for i=1:N2    
    pos=positions(i,1:3);   % coordinate of sphere center
    ipos=round(pos/dx+0.5); % index of sphere center location
    % Determine index range in which to operate
    xlo=round(ipos(1)-rad/dx-1); if(xlo<1  ); xlo=1 ; end;
    xhi=round(ipos(1)+rad/dx+1); if(nx <xhi); xhi=nx; end;
    ylo=round(ipos(2)-rad/dx-1); if(ylo<1  ); ylo=1 ; end;
    yhi=round(ipos(2)+rad/dx+1); if(ny <yhi); yhi=ny; end;
    zlo=round(ipos(3)-rad/dx-1); if(zlo<1  ); zlo=1 ; end;
    zhi=round(ipos(3)+rad/dx+1); if(nz <zhi); zhi=nz; end;
    for I=xlo:xhi; X=xrange(I);
    for J=ylo:yhi; Y=yrange(J);
    for K=zlo:zhi; Z=zrange(K);
        distancesquare=(X-pos(1))^2+(Y-pos(2))^2+(Z-pos(3))^2;
        if (distancesquare < rad^2); model(I,J,K)=1; end
    end
    end
    end
end

%% Save discretized geometry to model file

fp=fopen(strcat(modpath,'segmented-788x791x793.raw'),'w');
fwrite(fp,model,'uint8'); 
fclose(fp);

%% Matlab figures 

figure(1);clf; hist(conum(1:N),0:12); 
   xlim([-2 15]); title('Coordination number')
figure(2);clf; hist(radContact/rad,(0:4e-3:2e-1))
    xlim([-0.02 0.15]); title('Contact radius / Grain radius')
figure(3);clf; hist(overlap/(2*rad),0:4e-4:3e-2)
    xlim([-0.002 0.015]); title('Overlap / Grain diameter')
figure(4);clf;imagesc(model(:,:,1));axis square;

figsize=[300 200];
for i=1:4
    set(i,'position',[0 (i-1)*(figsize(2)+40)+50 figsize])
    colormap summer
end

print -f1 -r600 -depsc figure01.eps
print -f2 -r600 -depsc figure02.eps
print -f3 -r600 -depsc figure03.eps
print -f4 -r600 -depsc figure04.eps

%% Start the matVTK gui.

try 
    if (isdef_matVTKguiserver)
        % do nothing
    end
catch err
    ! /usr/local/matVTK/matVTK1.0/gui/linux64/matvtk_guiserver &
    pause(1);
    addpath('/usr/local/matVTK/matVTK1.0/matlab/');
    VTK=vtkinit();
    configS.backgroundColor = ones(1,3)*1.0;
    configS.cameraFocalPoint= boxdim/2;
    configS.cameraPosition=[2 -2 2].*boxdim;
    configS.cameraViewUp=[ 0 0 1];   
    vtkconfig(VTK,configS)
    clear hv hp hg ho;
    isdef_matVTKguiserver=true;
end

%% Plot sphere pack

% Remove objects from plot before replotting
try if(isnumeric(hv)); vtkdestroy(hv); end; end
try if(isnumeric(hp)); vtkdestroy(hp); end; end
try if(isnumeric(hg)); vtkdestroy(hg); end; end
try if(isnumeric(ho)); vtkdestroy(ho); end; end

% Plot volume data
hv=vtkplotvolume(VTK,model(1:4:end,1:4:end,1:4:end),...
    'opacity',1,'scale',4*dx*ones(1,3));

% Add boundary axes
hg=vtkgrid('gridFly','staticEdges', ...
    'gridXAxisLabelVisibility',0, ...
    'gridYAxisLabelVisibility',0, ...
    'gridZAxisLabelVisibility',0, ...
    'ambientColor', [0.9 0.1 0.1],...
    'volumeInterpolation','nearest');

% Overlay the volume plot with the sphere objects
points=((positions(1:N,1:3)));
pointLabels = linspace(1, 140, size(points, 1))';
pointColors = [0, 1, 0, 0;
                1, 1, 1, 1;
                0, 0, 0, 0];
hp=vtkplotpoints(VTK,points, pointLabels, pointColors, ...
    'pointSize',rad,...
    'thetaResolution', 16,...
    'phiResolution', 16,...
    'edgeColor',[0.5 0.5 0.5]);

% Include orientation arrows
ho=vtkorientation(VTK);    
