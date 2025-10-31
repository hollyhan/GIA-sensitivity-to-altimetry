function [md,md2ds]=world_mesher(varargin)
% This function creates a 3D mesh of the eath following the coastlines, based on a .25 degree resolution version of ETOPO1.
% Adapted from the code written by Lambert Caron
% Options (default), decription
%'model' (model()), model to plug mesh in
%'steps' ([1 2 3]), which parts of the routine are executed. 1:create coastline contours, 2:mesh 2d basins, 3:assemble 3d mesh
%'minpoint' (100), step 1 will discard contours with less than this amount of points prior to downscaling
%'smallconts' (200), step 1 will reduce contours that have less than this amount of points to triangles
%'proj' ([1 1]), projection for west and east hemisphere. 0:Sinusoidal, 1:Wagner VI, 2:Mollweide, 3:North Polar, 4:South Polar
%'coast_res' (250e3), coastline resolution in m
%'land_res' (700e3), continent resolution in m
%'ocean_res' (1000e3), ocean resolution in m
%'mergeres' (100), tolerance in m for points to be merged back in step 3
%'smoothness' (2), number of iteration to gradually accomodate coastline res to ocean/land
%'polecapS' (32), degrees of latitude from southpole where the antarctic basin will extend for meshing
%'arctic_divide', ([-180 -124.88 -114.08 -112.21  -108.8 -106.15 -102.61
%-99.082 -96.082 -95.682  -70.95 -60.622 -49.891 -2.2996  17.563  56.223
%74.153  108.88  135.07 168.67 180; 72  70.674  68.874  68.408 68.7 68.98
%68.68  69.474 72.073 74.206 74.206 67.275 52.814 64.211 74.14  76.541
%80.671 82.805   78.14  72.473 72]), [long;lat] of domain boundary for Arctic basin 
%'boundrefine' (2), number of iterations to make the basin boundaries seamless in terms of resoltion
%'extradom' (struct('lat',[],'long',[])), list of extra domains to be cut out from this mesh (so you can mesh it yourself an plug both back in a sea level model);
%'RefineMetric' (struct('H',[],'lat',[],'long',[])), boundary refinement metric
%'RefineError' (5), requested accuracy on RefineMetric
%'plotting' (false), triggers intermediary construction figures along the routine
%'coffset' (21.5), longitude of pacific divide relative to long=-180
%'coastline_contours' ([]), precomputed coastline contours, if step 1 is not run


	global plotting;

	options=pairoptions(varargin{:});
	
	load topo4.mat;

	[md,steps,minpoint,smallconts,coffset,proj,coastres,...
        land_res,ocean_res,mergeres, smoothness,polecapS,Nlim,...
        boundrefine,extradom,Refine,Refine_err,plotting,c2]=setdefaultoptions(options,...
        'model',model,...
        'steps',[1 2 3],...
        'minpoint',100,...
        'smallconts',200,...
        'coffset', 180-168.5,...
        'proj',[1 1],...
        'coast_res',250e3,...
        'land_res',700e3,...
        'ocean_res',1000e3,...
        'mergeres',100,...
        'smoothness',2,...
        'polecapS',32,...
        'arctic_divide', [-180 -124.88 -114.08 -112.21  -108.8 -106.15 -102.61...
        -99.082 -96.082 -95.682  -70.95 -60.622 -49.891 -2.2996  17.563  56.223...
        74.153  108.88  135.07 168.67 180; 72  70.674  68.874  68.408 68.7 68.98...
        68.68  69.474 72.073 74.206 74.206 67.275 52.814 64.211 74.14  76.541...
        80.671 82.805   78.14  72.473 72],...
        'boundrefine',2,...
        'extradom',struct('lat',[],'long',[]),...
        'RefineMetric',struct('H',[],'lat',[],'long',[]),...
        'RefineError', 5,...
        'plotting',false,...
        'coastline_contours',[]);

	if any(steps==1)
		%load('/Users/kyhan/Desktop/ISSM_Tutorial_Exercise/LIA-thickness/mesher_and_coastlines/coastline_coutours250km.mat')
		%c2=coastline_contours(coastres,minpoint,smallconts,coffset, extradom);
		load('coastline_coutours250km.mat')
        if ~any(steps>1)
            md=c2; %return contours instead of model
        end
    end

    
	if any(steps==2)

		disp('Meshing 2d');
		clear domain

		%atlantic divide
		lc=[0 0 0 -15 -25 -35 -35 -35 -35 -25 -20 -20 -20 -20 -20 -20 -20];
		lac=90-(0:length(lc)-1)*10;

		%Artic divide
        Nlim=Nlim';
% 	Nlim=[
% 	 -180.0000   72.0000
% 	 -124.8769   70.6738
% 	 -114.0789   68.8745
% 	 -112.2126   68.4081
%  	-108.8	    68.7
% 	 -106.15   68.98
% 	 -102.6141   68.68
% 	  -99.0816   69.4743
% 	  -96.0822   72.0733
% 	  -95.6823   74.2058
% 	  -76.9500+6   74.2058
% 	  -60.6221   67.2751
% 	  -49.8908   52.8138
% 	   -2.2996   61.2107+3
% 	   17.5634   72.1399+2
% 	   56.2229   69.5409+7
% 	   74.1529   77.6712+3
% 	  108.8798   77.8045+5
% 	  135.0750   73.1396+5
% 	  168.6688   72.4732
% 	  180.0000   72.0000];



		%[a,b]=interpsegment(Nlim(2:11,2)', Nlim(2:11,1)',0.25);
		%Nlim=[Nlim(1,:);[b' a'];Nlim(12:end,:)];

		r=md.solidearth.planetradius;
		
		if coffset~=0
			longi=-180+coffset;
			lati=interp1(Nlim(:,1),Nlim(:,2),longi);
			Nlim(end,:)=[];
			Nlim(end+1,:)=[longi lati];
			Nlim(Nlim(:,1)<longi,1)=Nlim(Nlim(:,1)<longi,1)+360;
			[~,ind]=sort(Nlim(:,1));
			Nlim=Nlim(ind,:);
			Nlim(end+1,:)=[longi+360 lati];
		end


		thres=ocean_res*180/pi/r;
		N=ceil((180-Nlim(1,2)-polecapS)/thres);

		exth=[];
		tol=0.1;
		for i=1:length(c2);
			exth=[exth, 90-c2(i).lat(c2(i).long+coffset<tol | c2(i).long>coffset+360-tol)];
		end
		exth(exth<90-Nlim(1,2))=[];
		exth(exth>180-polecapS)=[];

		exth=unique(exth);

		th=[Nlim(1,2) 90-sort(exth) -90+polecapS];

		[th,~]=interpsegment(th,th*0-180,thres);

		th=90-th;


		md2ds={};
		md2ds{1}=model;
		md2ds{2}=model;


	%%%%%%%% segments

		M2=ceil((lc(end)+180)*sind(th(end))*2/thres);

		[longi,lati]=polyxpoly(lc,lac,Nlim(:,1)',Nlim(:,2)');
		[longj,latj]=polyxpoly(lc,lac,[-180;180],-90+[polecapS polecapS]);

		seg1.long=th*0-180+coffset;
		seg1.lat=90-th;

		seg2.long=linspace(-180+coffset,lc(end),M2);
		seg2.lat=90-[seg2.long*0+th(end)];

		pos=find(lac<lati & lac>latj);
		pos=pos(end:-1:1);
		[seg3.lat,seg3.long]=interpsegment([latj lac(pos) lati], [longj,lc(pos),longi],thres);
		
		pos=find(Nlim(:,1)<longi);
		pos=pos(end:-1:1);
		[seg4.lat,seg4.long]=interpsegment([lati Nlim(pos,2)'], [longi,Nlim(pos,1)'],thres);

		M2=ceil((180-lc(end))*sind(th(end))*2/thres);

		seg5.long=linspace(lc(end),180+coffset,M2);
		seg5.lat=90-[seg5.long*0+th(end)];

		pos=find(Nlim(:,1)>longi);
		pos=pos(end:-1:1);

		[seg6.lat,seg6.long]=interpsegment([ Nlim(pos,2)' lati], [Nlim(pos,1)' longi],thres/2);

		domain1.long=[seg1.long seg2.long(2:end-1) seg3.long seg4.long(2:end-1) seg1.long(1)]';
		domain1.lat=[seg1.lat seg2.lat(2:end-1) seg3.lat seg4.lat(2:end-1) seg1.lat(1)]';

		domain2.long=[seg3.long(end:-1:1), seg5.long(2:end-1) seg1.long(end:-1:1)+360 seg6.long(2:end-1) seg3.long(end)]';
		domain2.lat=[seg3.lat(end:-1:1), seg5.lat(2:end-1) seg1.lat(end:-1:1) seg6.lat(2:end-1) seg3.lat(end)]';

		domainN.long=[seg4.long seg6.long(2:end)]';
		domainN.lat=[seg4.lat seg6.lat(2:end)]';

		domainS.long=[seg2.long seg5.long(2:end)]';
		domainS.lat=[seg2.lat seg5.lat(2:end)]';
		
		if plotting
			hold on
			plot(seg1.long,seg1.lat,'k',seg2.long,seg2.lat,'k',seg3.long,seg3.lat,'k',seg4.long,seg4.lat,'k',seg5.long,seg5.lat,'k',seg6.long,seg6.lat,'k')
		end

		md2ds{1}=mesh_partition(md2ds{1},domain1,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2, proj(1),90);
		md2ds{2}=mesh_partition(md2ds{2},domain2,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2, proj(2),-90);
		md2ds{3}=mesh_partition(model,domainN,extradom,Refine,Refine_err,coastres/2, land_res/2, ocean_res/2,smoothness, c2,3,0);
		md2ds{4}=mesh_partition(model,domainS,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2,4,0);

		for f=1:boundrefine
		%Refine boundaries and mesh again
		res1=boundtest(md2ds{1},smoothness)*180/pi;
		res2=boundtest(md2ds{2},smoothness)*180/pi;
		res3=boundtest(md2ds{3},smoothness)*180/pi;
		res4=boundtest(md2ds{4},smoothness)*180/pi;

		l1=length(seg1.long)-1;
		l2=length(seg2.long)-1;
		l3=length(seg3.long)-1;
		l4=length(seg4.long)-1;
		l5=length(seg5.long)-1;
		l6=length(seg6.long)-1;

		seg1res=[res1(1:l1) res2(l3+l5+ (l1:-1:1))];
		seg2res=[res1(l1+(1:l2)) res4(1:l2)];
		seg3res=[res1(l1+l2+(1:l3)) res2(l3:-1:1)];
		seg4res=[res1(l1+l2+l3+(1:l4)) res3(1:l4)];
		seg5res=[res2(l3+(1:l5)) res4(l2+(1:l5))];
		seg6res=[res2(l3+l5+l1+(1:l6)) res3(l4+(1:l6))];

		[seg1.lat,seg1.long]=interpsegment(seg1.lat,seg1.long,min(seg1res,[],2));
		[seg2.lat,seg2.long]=interpsegment(seg2.lat,seg2.long,min(seg2res,[],2));
		[seg3.lat,seg3.long]=interpsegment(seg3.lat,seg3.long,min(seg3res,[],2));
		[seg4.lat,seg4.long]=interpsegment(seg4.lat,seg4.long,min(seg4res,[],2));
		[seg5.lat,seg5.long]=interpsegment(seg5.lat,seg5.long,min(seg5res,[],2));
		[seg6.lat,seg6.long]=interpsegment(seg6.lat,seg6.long,min(seg6res,[],2));

		domain1.long=[seg1.long seg2.long(2:end-1) seg3.long seg4.long(2:end-1) seg1.long(1)]';
		domain1.lat=[seg1.lat seg2.lat(2:end-1) seg3.lat seg4.lat(2:end-1) seg1.lat(1)]';

		domain2.long=[seg3.long(end:-1:1), seg5.long(2:end-1) seg1.long(end:-1:1)+360 seg6.long(2:end-1) seg3.long(end)]';
		domain2.lat=[seg3.lat(end:-1:1), seg5.lat(2:end-1) seg1.lat(end:-1:1) seg6.lat(2:end-1) seg3.lat(end)]';

		domainN.long=[seg4.long seg6.long(2:end)]';
		domainN.lat=[seg4.lat seg6.lat(2:end)]';

		domainS.long=[seg2.long seg5.long(2:end)]';
		domainS.lat=[seg2.lat seg5.lat(2:end)]';


		md2ds{1}=mesh_partition(md2ds{1},domain1,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2, proj(1),90);
		md2ds{2}=mesh_partition(md2ds{2},domain2,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2, proj(2),-90);
		md2ds{3}=mesh_partition(model,domainN,extradom,Refine,Refine_err,coastres/2, land_res/2, ocean_res/2,smoothness, c2,3,0);
		md2ds{4}=mesh_partition(model,domainS,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2,4,0);
		end

		Npole=md2ds{3};
		Spole=md2ds{4};

	end

	if any(steps==3)

		disp('Projecting to 3d');
		md=model;md.mesh=mesh3dsurface();
		[lat1,long1]=revproje(md2ds{1}.mesh.x,md2ds{1}.mesh.y,proj(1));
		[lat2,long2]=revproje(md2ds{2}.mesh.x,md2ds{2}.mesh.y,proj(2));
		md.mesh.lat=[lat1;lat2;Npole.mesh.lat;Spole.mesh.lat];
		md.mesh.long=[long1-90;long2+90;Npole.mesh.long;Spole.mesh.long];
		r=md.solidearth.planetradius;
		md.mesh.x=r*cosd(md.mesh.lat).*cosd(md.mesh.long);
		md.mesh.y=r*cosd(md.mesh.lat).*sind(md.mesh.long);
		md.mesh.z=r*sind(md.mesh.lat);
		md.mesh.r=md.mesh.x*0+r;
		md.mesh.numberofelements=md2ds{1}.mesh.numberofelements+md2ds{2}.mesh.numberofelements+Npole.mesh.numberofelements+Spole.mesh.numberofelements;

		nv=[0 md2ds{1}.mesh.numberofvertices,md2ds{2}.mesh.numberofvertices,Npole.mesh.numberofvertices,Spole.mesh.numberofvertices];
		ne=[0 md2ds{1}.mesh.numberofelements,md2ds{2}.mesh.numberofelements,Npole.mesh.numberofelements,Spole.mesh.numberofelements];
		cnv=cumsum(nv);
		
		disp(['Vertices in basins: ' num2str(nv(2:end))]);
		disp(['Elements in basins: ' num2str(ne(2:end))]);
		disp(['Total mesh size: ', num2str(md.mesh.numberofelements), ' elements, ', num2str(cnv(end)), ' vertices']);

		md.mesh.numberofvertices=cnv(end);

		md.mesh.elements=[md2ds{1}.mesh.elements; md2ds{2}.mesh.elements+cnv(2);Npole.mesh.elements+cnv(3); Spole.mesh.elements+cnv(4)];

		md.mask.ocean_levelset=interp2([topop-180-360 topop-180 topop+180],90-topoth,double([topo topo topo]),md.mesh.long,md.mesh.lat);

		%plotmodel(md,'data', md.mask.ocean_levelset, 'edgecolor', 'k')

		[md,merged]=removeduplicatevertices(md,mergeres);
		
		%plotmodel(md,'data',test,'edgecolor', 'k');


		md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
		md.mesh.area=GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);

		md.mask.ocean_levelset=interp2([topop-180-360 topop-180 topop+180],90-topoth,double([topo topo topo]),md.mesh.long,md.mesh.lat);
	
		% Set topography same as the ocean levelset. But this is not always the case.
		md.geometry.bed = md.mask.ocean_levelset;
		md.geometry.surface = zeros(md.mesh.numberofvertices, 1);
		md.geometry.thickness = zeros(md.mesh.numberofvertices, 1);
		md.geometry.base = zeros(md.mesh.numberofvertices, 1);
	
		if plotting
			plotmodel(md,'data', md.mask.ocean_levelset, 'edgecolor', 'k', 'caxis', [-1 1])
		end

		res=zeros(md.mesh.numberofelements,1);
		d=[0 0 0];
		for e=1:md.mesh.numberofelements;
			for i=1:3;
				i1=md.mesh.elements(e,i);
				i2=md.mesh.elements(e,mod(i+1,3)+1);
				d(i)=great_circle(md.mesh.long(i1),md.mesh.lat(i1),md.mesh.long(i2),md.mesh.lat(i2));
			end
			res(e)=mean(d);
		end

		if plotting
			plotmodel(md,'data', res*r/1e3, 'edgecolor', 'k');
		end
	end
	clear plotting
end

function [res2,bounds]=boundtest(md,smoothness);

	bounds=find(md.mesh.vertexonboundary);
	bounds(end+1)=bounds(1);

	res2=zeros(length(bounds)-1,1);

	for i=1:length(bounds)-1;
		if smoothness>0
		bels{i}=md.mesh.vertexconnectivity(bounds(i:i+1),:);
		bels{i}=unique(bels{i}(:));bels{i}(bels{i}==0)=[];

		pels=[];
		for j=1:smoothness
			v=md.mesh.elements([pels; bels{i}],:);v=unique(v(:));v=setdiff(v,bounds(i:i+1));
			pels=md.mesh.vertexconnectivity(v,:);pels=unique(pels(:));pels(pels==0)=[];pels=setdiff(pels,bels{i});
		end

		res=zeros(length(pels),1);
		d=[0 0 0];
		for pe=1:length(pels);
			e=pels(pe);
			for k=1:3;
				i1=md.mesh.elements(e,k);
				i2=md.mesh.elements(e,mod(k,3)+1);
				d(k)=great_circle(md.mesh.long(i1),md.mesh.lat(i1),md.mesh.long(i2),md.mesh.lat(i2));
			end
			res(pe)=mean(d);
		end
		res2(i)=mean(res);
		else
		i1=bounds(i);i2=bounds(i+1);
		res2(i)=great_circle(md.mesh.long(i1),md.mesh.lat(i1),md.mesh.long(i2),md.mesh.lat(i2));
		end
	end
end

function d=great_circle(lon1,lat1,lon2,lat2)
	dl=lon2-lon1;
	dt=lat2-lat1;

	%d=sqrt(dl.^2+dt.^2);return;

	a=sind(dt/2).^2 + (1 - sind(dt/2).^2-sind(0.5*(lat1+lat2)).^2).*sind(dl/2).^2;
	d=2*asin(sqrt(a));
end

function [md,merged]=removeduplicatevertices(md,mergeres)

	dups={};
	merged=[];
	tol=mergeres;%m
	for i=1:md.mesh.numberofvertices;
		switch class(md.mesh)
		case 'mesh3dsurface'
			ind=find(abs(md.mesh.x-md.mesh.x(i))<tol & abs(md.mesh.y-md.mesh.y(i))<tol & abs(md.mesh.z-md.mesh.z(i))<tol);
		case 'mesh2d'
			ind=find(abs(md.mesh.x-md.mesh.x(i))<tol & abs(md.mesh.y-md.mesh.y(i))<tol);
		end
		if length(ind)>1
		dups{end+1}=ind;
		end
	end

	vertices=1:md.mesh.numberofvertices;
	del=[];
	for i=1:length(dups)
		vi=dups{i}(1);
		vd=dups{i}(2:end);
		for j=1:length(vd)
			ind=md.mesh.elements(:)==vd(j);
			md.mesh.elements(ind)=vi;
		end
		del=[del; vd];
		merged=[merged;vi];
	end

	merged=setdiff(merged,del);

	del=unique(del);
	vertices(del)=[];
	md.mesh.x(del)=[];
	md.mesh.y(del)=[];
	md.mesh.lat(del)=[];
	md.mesh.long(del)=[];
	md.mesh.r(del)=[];

	switch class(md.mesh)
		case 'mesh3dsurface'
			md.mesh.z(del)=[];
	end

	md.mesh.numberofvertices=length(vertices);

	elements2=md.mesh.elements;
	oldv=vertices;

	for i=1:length(merged)
		merged(i)=find(oldv==merged(i));
	end

	del=[];
	for i=1:md.mesh.numberofelements
		i1=find(oldv==md.mesh.elements(i,1));
		i2=find(oldv==md.mesh.elements(i,2));
		i3=find(oldv==md.mesh.elements(i,3));
		elements2(i,:)=[i1 i2 i3];
		if length(unique(elements2(i,:)))<3
			del=[del i];
		end
	end
	md.mesh.numberofelements=md.mesh.numberofelements-length(del);
	elements2(del,:)=[];

	md.mesh.elements=elements2;

end

function c2=coastline_contours(coastres,minpoint,smallconts,coffset,extradom)

	load topo4.mat;
	global plotting
	
	ind=find(topop-180<=-180+coffset);
	topo=topo(:,[ind(end)+1:end ind]);

	disp('Defining coastlines');
	cont=contourc(topop,90-topoth,double(topo),[0 0]);

	j=1;
	i=0;

	%Assembling contour structure
	c=struct();
	while j<length(cont)

		l=cont(2,j);
	
		i=i+1;
		c(i).long=cont(1,j+1:j+l);
		c(i).lat=cont(2,j+1:j+l);

		j=j+l+1;

	end


	%Removing tiny contours
	for i=1:length(c);
		cl(i)=length(c(i).lat);
	end
	x=find(cl<minpoint);
	c2=c;
	c2(x)=[];
	cl(x)=[];
	x=find(cl<smallconts);
    
    %% Removing landlocked contours
    for i=1:length(c2)
        for j=1:length(c2)
            if i~=j
            in=inpolygon(c2(i).long,c2(i).lat,c2(j).long,c2(j).lat);
            if sum(in)>=1
                del(i)=1;
            end
            end
        end
    end
    c2(del==1)=[];
    cl(del==1)=[];
    %%
    

		disp('Downscaling contours')

	for ind=1:length(c2);
		if ismember(ind,x)
			%Deal with islands and small contours: just get a triangle
			pos=[1 floor([cl(ind)/3 cl(ind)/3*2]) cl(ind)];
			pos=unique(pos);
			pos(pos==0)=[];
			c2(ind).lat=c2(ind).lat(pos);
			c2(ind).long=c2(ind).long(pos);
		else
			d=zeros(length(c2(ind).lat),1);
			for i=1:length(c2(ind).lat)-1;
				d(i)=great_circle(c2(ind).long(i),c2(ind).lat(i),c2(ind).long(i+1),c2(ind).lat(i+1));
			end

			d(end)=great_circle(c2(ind).long(end),c2(ind).lat(end),c2(ind).long(1),c2(ind).lat(1));

			d=d*6371e3;

			if d(end)==0;
				d(end)=[];
			end

			n=round(sum(d)/coastres);
			res=sum(d)/n;

			l=res:res:sum(d);

			pos=[1 round(interp1(cumsum(d),1:length(d),l))];
			pos=unique(pos);
			c2(ind).lat=c2(ind).lat(pos);
			c2(ind).long=c2(ind).long(pos);
		end
	end

	%checking whether this contour intersects another one
	del=[];
	for i=1:length(c2);
		ind=[];
		for j=1:i-1
			test=polyxpoly(c2(i).long,c2(i).lat,c2(j).long,c2(j).lat);
			if length(test>0)
				ind=[ind j];
			end
		end

		%delete subdomain i if it intersects several other subdomains
		% if not delete the subdomain with the least # of points
		if length(ind)>0
			k=i;
			if (length(ind)==1)
				if length(c2([ind]).lat)<length(c2(i).lat)
					k=ind;
				end
			end
			del=[del k];
		end

		%flipping segments that intersect
		fixed=false;
		nit=0;
		while ~fixed
			xi=c2(i).long;
			yi=c2(i).lat;
			[c2(i).long,c2(i).lat]=fixintersections(xi,yi);
			if sum(abs(c2(i).long-xi))==0 && sum(abs(c2(i).lat-yi))==0
				fixed=true;
			end
			nit=nit+1;
		end
		if nit>1
		disp(['contour ' num2str(i) ' : ' num2str(nit-1) ' iterations']);
		end
	end

	c2(del)=[];
	%removing interstection with domains;
	for i=1:length(c2);
		for j=1:length(extradom)
            if ~isempty(extradom(j).lat)
			[c2(i).long,c2(i).lat]=subtractpoly(c2(i).long,c2(i).lat,extradom(j).long+180-coffset,extradom(j).lat);
            end
		end
	end

	%Downgrading resolution
	%for i=1:length(c2);
	%	c2(i).lat=[c2(i).lat([1:skip:end-1 1])];
	%	c2(i).long=[c2(i).long([1:skip:end-1 1])];
	%end


	%Reassembling contour array and removing duplicates
	cols2=jet(length(c2));
	cont2=zeros(2,0);
	if plotting
	 figure(1)
	end
	for i=1:length(c2);
        c2(i).long=c2(i).long+coffset;
        if plotting
        	plot(c2(i).long-180,c2(i).lat,'color',cols2(i,:));
        end
        %cont2=[cont2 [c2(i).long;c2(i).lat]];
    end

    
    
    
% 	dup=[];
% 	for i=1:length(cont2)
% 		ind=find(cont2(1,:)==cont2(1,i) & cont2(2,:)==cont2(2,i));
% 		dup=[dup ind(2:end)];
% 	end
% 
% 	cont2(:,dup)=[];
% 
% 	del=find(cont2(1,:)==0);
% 	cont2(:,del)=[];

end

function md2d=mesh_partition(md2d,domain,extradom,Refine,Refine_err,coastres, land_res, ocean_res,smoothness, c2, proj,longoffset)
	global plotting
	r=md2d.solidearth.planetradius;
	inres=[coastres, ocean_res, land_res]*180/pi/r;
	inih=[min(inres) max(inres)];

	domains=struct;

	[domains(1).x,domains(1).y]=proje(domain.lat,domain.long+longoffset,proj);
	domains(1).nods=length(domains(1).x);

	j=1;
	for i=1:length(extradom);
		[x,y]=proje(extradom(i).lat,extradom(i).long+longoffset,proj);
		in=inpolygon(x, y, domains(1).x,domains(1).y);
		if sum(in)>1
			j=j+1;
			domains(j).x=x;
			domains(j).y=y;
			domains(j).nods=length(x);
		end
	end


	md2d=bamg(md2d, 'hmin', inih(1), 'hmax', inih(2), 'domain', domains,'NoBoundaryRefinement',1);

	%plotmodel(md2d,'data','mesh')
	conts=[];
	%add contours in main domain
	for i=1:length(c2)
		[x,y]=proje(c2(i).lat,c2(i).long-180+longoffset,proj);


		in=inpolygon(x, y, domains(1).x,domains(1).y);

		%remove them if they are inside extra domains, to be cut out from the mesh
		cut=in*0;
		for j=2:length(domains)
			cut=cut+inpolygon(x, y, domains(j).x,domains(j).y);
		end
		in(cut>0)=0;

		if sum(in)>1
			conts=[conts, i];
		end
	end

	cc=c2(conts);
	subdom=struct;
	del=[];

	for i=1:length(cc);
		[subdom(i).x,subdom(i).y]=proje(cc(i).lat(end:-1:1)',cc(i).long(end:-1:1)'-180+longoffset,proj);
		subdom(i).nods=length(cc(i).lat);

		%checking if points lie on the domain boundary, nudging them in if necessary
		test=subdom(i).x*0;
		[in,on]=inpolygon(subdom(i).x,subdom(i).y,domains(1).x,domains(1).y);
		test(in)=1;
		test(on)=0;
		out=find(test==0);
		subdom(i).x(out)=[];
		subdom(i).y(out)=[];
		subdom(i).nods=length(subdom(i).x);

		for j=2:length(domains)
			%remove points inside extra subdomains
			test=subdom(i).x*0;
			[in,on]=inpolygon(subdom(i).x,subdom(i).y,domains(j).x,domains(j).y);
			test(in)=1;
			out=find(test==1);
			subdom(i).x(out)=[];
			subdom(i).y(out)=[];
			subdom(i).nods=length(subdom(i).x);
		end

		%checking that last point of contour is also the first point
		if subdom(i).x(1)~=subdom(i).x(end) | subdom(i).y(1)~=subdom(i).y(end)
			subdom(i).x(end)=subdom(i).x(1);
			subdom(i).y(end)=subdom(i).y(1);
			subdom(i).nods=length(subdom(i).x);
		end

		subdom(i).closed=1;subdom(i).density=1;subdom(i).name='coast';subdom(i).Geometry='Polygon';
		subdom(i).BoundingBox=[min(subdom(i).x) min(subdom(i).y); max(subdom(i).x) max(subdom(i).y)];
		%hold on
		%plot(subdom(i).x,subdom(i).y,'linewidth',1.5);

		if subdom(i).nods<4;
			del=[del i];
		end

	end

	subdom(unique(del))=[];

	try
		md2d=bamg(md2d, 'hmin', inih(1), 'hmax', inih(2), 'domain', domains, 'subdomains',subdom,'KeepVertices',1,'NoBoundaryRefinementAllBoundaries',1);
	catch
		if plotting
			plotmodel(md2d,'data','mesh');
			hold on;
			for i=1:length(subdom);
				plot(subdom(i).x, subdom(i).y, 'linewidth',1.5);
			end
		end
		error('Problem meshing subdomains (usually that happens when coastlines intserect the mesh domain)');
	end


	land=md2d.mesh.x*0-1;
	for ind=1:length(cc);
		[x,y]=proje(cc(ind).lat,cc(ind).long-180+longoffset,proj);
		[in,on]=inpolygon(md2d.mesh.x,md2d.mesh.y,x,y);
		land(in & ~on)=1;
		%land(on)=1;
		land(on)=max([land(on) zeros(size(land(on)))],[],2);
	end



	r=md2d.solidearth.planetradius;
	Resolution=land*NaN;
	Resolution(land==1)=land_res *180/pi/r;
	Resolution(land==-1)=ocean_res *180/pi/r;

	md2d.mesh.vertexconnectivity=NodeConnectivity(md2d.mesh.elements,md2d.mesh.numberofvertices);

	for i=1:smoothness
		pos=find(Resolution~=land_res *180/pi/r & Resolution~=ocean_res *180/pi/r);
		v=md2d.mesh.vertexconnectivity(pos,:);
		v(v==0)=[];
		v=unique(v(:));
		v=md2d.mesh.elements(v,:);
		v=unique(v(:));
		Resolution(v)=0.5*(Resolution(v)+coastres*180/pi/r);
	end

	md2d.mesh.vertexonboundary=zeros(md2d.mesh.numberofvertices,1);
	for j=1:length(domains)
	for i=1:length(domains(j).x);
		v=find(md2d.mesh.x==domains(j).x(i) & md2d.mesh.y==domains(j).y(i));
		md2d.mesh.vertexonboundary(v)=1;
	end
	end

	skip =true;

	[md2d.mesh.lat,md2d.mesh.long]=revproje(md2d.mesh.x,md2d.mesh.y,proj);
	md2d.mesh.long=md2d.mesh.long-longoffset;
	bres=boundtest(md2d,0)*2*180/pi;
	v=find(md2d.mesh.vertexonboundary);
	Resolution(v)=bres;
	if ~skip
	for i=1:smoothness
		els=md2d.mesh.vertexconnectivity(v,:);
		els=unique(els);els(els==0)=[];
		v=md2d.mesh.elements(els,:);v=unique(v(:));
		Resolution(v)=NaN;
	end
	end

	%md2d=bamg(md2d, 'subdomains', subdom,'hVertices',Resolution,'KeepVertices',1,'MaxCornerAngle',1e-15,'NoBoundaryRefinement',1);	


	[Refine.x,Refine.y]=proje(Refine.lat,Refine.long+longoffset,proj);
	Hi=griddata(Refine.x,Refine.y,Refine.H,md2d.mesh.x,md2d.mesh.y);
	%Resolution(~isnan(Hi))=NaN;

	R=Resolution;R(isnan(R))=-1;
	if plotting
		plotmodel(md2d,'data',R, 'figure', get(gcf,'Number')+1,'edgecolor','k')
	end

	md2d=bamg(md2d,'subdomains', subdom, 'field', Hi, 'err', Refine_err,'gradation', 1.15,'hmaxVertices',Resolution,'KeepVertices',1,'MaxCornerAngle',1e-15,'NoBoundaryRefinement',1);

	md2d.mesh.vertexconnectivity=NodeConnectivity(md2d.mesh.elements,md2d.mesh.numberofvertices);

	md2d.mesh.vertexonboundary=zeros(md2d.mesh.numberofvertices,1);
	for i=1:length(domains(1).x);
		v=find(md2d.mesh.x==domains(1).x(i) & md2d.mesh.y==domains(1).y(i));
		md2d.mesh.vertexonboundary(v)=1;
	end

	[md2d.mesh.lat,md2d.mesh.long]=revproje(md2d.mesh.x,md2d.mesh.y,proj);
	md2d.mesh.long=md2d.mesh.long-longoffset;

end

	
function [lat2,long2]=interpsegment(lat,long,res);

	if length(res)==1
		res=repmat(res, [length(lat)-1,1]);
	end	

	lat2=[];
	long2=[];
	for i=1:length(lat)-1
		s1=i:i+1;
		d=great_circle(long(i),lat(i),long(i+1),lat(i+1))*180/pi;
		N=round(d/res(i));
		if N==0;N=1;end
		
		lati=interp1([0 1], lat(s1), linspace(0,1,N+1));
		longi=interp1([0 1], long(s1), linspace(0,1,N+1));
		lat2=[lat2 lati(1:end-1)];
		long2=[long2 longi(1:end-1)];
	end
	lat2(end+1)=lat(end);
	long2(end+1)=long(end);
end
	
function [xi,yi]=fixintersections(x,y)

	xi=x;yi=y;

	if length(x)<=4
		return;
	end

	ind=[];

	N=length(xi);
	for i=1:N-1;
		s1=[i i+1];
		if i==1
			M=N-2;
		else
			M=N-1;
		end
		for j=i+2:M;
			s2=[j j+1];

			%checkintersection			
			a1 = diff(y(s1)) / diff(x(s1)); 
			b1 = y(s1(1)) - a1*x(s1(1));
			a2 = diff(y(s2)) / diff(x(s2)); 
			b2 = y(s2(1)) - a2*x(s2(1));
			X = (b1-b2)/(a2-a1);
	
			test= X>min(x(s1)) && X<max(x(s1)) && X>min(x(s2)) && X<max(x(s2));
			if (test)
				ind=[ind;s1 s2];
			end
		end
	end

	%ind(ind==N)=1

	for i=1:size(ind,1)
		xi(ind(i,:))=xi(ind(i,[1 3 2 4]));
		yi(ind(i,:))=yi(ind(i,[1 3 2 4]));
	end

end

function [lat,long] = revproje(x,y,projid)
	switch projid
		case 0
			lat=y;
			long=x./sind(90-lat);

		case 1
			psi=asin(sqrt(3)/180*y);
			long=x./cos(psi);
			lat=y;
		case 2
			th=asin(y/sqrt(2)/180*pi);
			long=x./(2*sqrt(2)*cos(th));
			lat=180/pi*asin((2*th+sin(2*th))/pi);
		case 3
			lat=90-sqrt(x.^2+y.^2);
			long=180/pi*atan2(y,x);
		case 4
			lat=-90+sqrt(x.^2+y.^2);
			long=180/pi*atan2(y,x);
	end
end

function [x,y] = proje(lat,long,projid)
	switch projid
		case 0 %Sinusoidal
			x=sind(90-lat).*long;
			y=lat;
		case 1 %Wagner VI
			x=long.*sqrt(1-3*(lat/180).^2);
			y=lat;
		case 2 %Mollweide
			th=lat*pi/180;
			for i=1:100
				th=th-(2*th+sin(2*th)-pi*sind(lat))./(4*cos(th).^2);
			end
			x=2*sqrt(2).*long.*cos(th);
			y=sqrt(2)*sin(th)*180/pi;
		case 3 %North Polar
			x=(90-lat).*cosd(long);
			y=(90-lat).*sind(long);
		case 4 %South Polar
			x=(90+lat).*cosd(long);
			y=(90+lat).*sind(long);
	end
end
		
function	varargout=setdefaultoptions(options,varargin);
		N=length(varargin)/2;
		for j=0:N-1;
			if exist(options,varargin{2*j+1})
				varargout{j+1}=getfieldvalue(options,varargin{2*j+1});
			else
				varargout{j+1}=varargin{2*j+2};
			end
		end
end

function [x1,y1]=subtractpoly(x1,y1,x2,y2)

	global plotting
	if plotting
		hold off
		plot(x1,y1,x2,y2)
	end

	%transpose contour if not column vectors
	flip1=false;
	if size(x1,1)==1
		flip1=true;
		x1=x1';
		y1=y1';
	end

	if size(x2,1)==1
		x2=x2';
		y2=y2';
	end

	x1i=x1;y1i=y1; %backup


	%create a contour around contour to create a buffer zone around it: this is where new points will be place on contour 1
	dy=diff(y2);
	dx=diff(x2);
	
	dn=sqrt(dx.^2+dy.^2);
	if max(dn>0.5);
		dx=dx*0.5/max(dn);
		dy=dy*0.5/max(dn);
	end
	

	xm=0.5*(x2(2:end)+x2(1:end-1));
	ym=0.5*(y2(2:end)+y2(1:end-1));


		%figure where the outside normal is
	ox=1;
	oy=1;

	xn=xm+oy*dy;
	yn=ym+ox*dx;
    
	in1=inpolygon(xn,yn,x2,y2);

	if sum(in1)>0
		ox=-1;
	xn=xm+oy*dy;
	yn=ym+ox*dx;
	in1=inpolygon(xn,yn,x2,y2);
	end

	if sum(in1)>0
		oy=-1;
	xn=xm+oy*dy;
	yn=ym+ox*dx;
	in1=inpolygon(xn,yn,x2,y2);
	end

	if sum(in1)>0
		ox=1;oy=-1;
	xn=xm+oy*dy;
	yn=ym+ox*dx;
	in1=inpolygon(xn,yn,x2,y2);
	end

	in2=inpolygon(x1,y1,x2,y2);
	in1=inpolygon(x2,y2,x1,y1);
	if plotting
		hold on
		plot(x1(in2),y1(in2),'o',x2(in1),y2(in1),'o')
	end

	% find contour 1 and 2 intersections
	[xi,yi,ii,jj]=polyxpoly(x1,y1,x2,y2);
    ii=[floor(ii) floor(jj)];
%     [x1 y1]
%     [x2 y2]
%     ii
	if mod(size(ii,1),2)~=0
		error ('Error in subtractpoly: odd number of intersections, are the contours closed?')
	end

	%for each pair of intersection (we go in and out of contour 2), edit points of contour 1 to follow contour 2 so we stay outside of it
	scale=1.5;
    
    find(in1)
	while ~isempty(ii)
	if plotting
		plot(x1,y1,'.-')
	end

	if ii(2,1) < ii(1,1)
		ii(1:2,:)=ii([2 1],:);
		xi([1 2])=xi([2 1]);
		yi([1 2])=yi([2 1]);
	end

	%hold on; plot(x1(ii(:,1)),y1(ii(:,1)),'o')
    %ii
	s=1;o=1;
	if ii(1,2)>ii(2,2); s=-1; o=-1; end

    path1=ii(1,2):s:ii(2,2);

	%take the shortest # of sergments
	if s==1
		path2=[ii(2,2):s:size(x2,2) 1:s:ii(1,2)];
	else
		path2=[ii(2,2):s:1 size(x2,2):s:ii(1,2)]; 
    end

	if isempty(path1) && isempty(path2)
		%if we are here, we have no points of contour 1 inside 2 but segments ...
        % still intersect so we are toast;
		%increase the size of the buffer through "scale" and try again.
		scale=scale*1.1;
		x1=x1i;y1=y1i;
		[xi,yi,ii,jj]=polyxpoly(x1,y1,x2,y2);
        ii=[floor(ii) floor(jj)];
		continue
    else
        if isempty(path1)
            ind=path2;
        elseif isempty(path2)
            ind=path1;
        else
            if length(path2)<length(path1)
                ind=path2;
            else
                ind=path1;
            end
        end
    end

    if sum(in1(ind(2:end))==0) %if the segment we are about to add
                                %from polygon 2 is not in polygon 1 at all,
                                % we are starting from the wrong
                                %intersection, get to the next one and try
                                %again
        ii=ii([2:end 1],:);
        continue;
    end

    
%     if length(ind)>2
        imid=ind(2:end-1);
        xp=[xi(1)+oy*dy(ind(1))*scale;
        xm(imid)+oy*dy(imid)*scale;
        xi(2)+oy*dy(ind(end))*scale;
        ];
        yp=[yi(1)+ox*dx(ind(1))*scale;
        ym(imid)+ox*dx(imid)*scale;
        yi(2)+ox*dx(ind(end))*scale;
        ];
%     else
% 
%         xp=xm(ind)+oy*dy(ind)*scale;
%         yp=ym(ind)+ox*dx(ind)*scale;
%     end

	if plotting
		plot(xp,yp,'x-')
	end
    x1([ii(1,1) ii(2,1)+1])
    y1([ii(1,1) ii(2,1)+1])

	x1=[x1(1:ii(1,1)); xp ; x1(ii(2,1)+1:end)];
	y1=[y1(1:ii(1,1)); yp ; y1(ii(2,1)+1:end)];


	[xi,yi,ii,jj]=polyxpoly(x1,y1,x2,y2);
    ii=[floor(ii) floor(jj)];

	end

	if flip1
		x1=x1';
		y1=y1';
	end

end
