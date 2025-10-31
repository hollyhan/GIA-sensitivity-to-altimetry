function [md, transfer] = mergeMeshes(md_regional_3D, md_global_extruded_3D)
    % Adapted from a code written by Lambert Caron
    % This function merges a 3D regional mesh with a 3D global mesh.
    %
    % Inputs:
    %    - md_regional_3D: model object with 3D regional mesh.
    %    - md_global_extruded_3D: model object with 3D global mesh where Greenland mesh is extruded
    % Outputs:
    %    - md: merged global 3D mesh model
    %    - transfer: the list of indices to refer to the index of regional mesh on the merged global mesh, size of md_global_extruded_3D.mesh.numberofverticies+md_regional_3D.mesh.numberofvertices
	% col 1 is for unmerged mesh, col2 is merged mesh. output of the function 'removeduplicatevertices'

	md=model;
    md.mesh=mesh3dsurface();
	r=md.solidearth.planetradius;

	md.mesh.numberofelements=md_global_extruded_3D.mesh.numberofelements+md_regional_3D.mesh.numberofelements;
	md.mesh.x=[md_global_extruded_3D.mesh.x; md_regional_3D.mesh.x];
	md.mesh.y=[md_global_extruded_3D.mesh.y; md_regional_3D.mesh.y];
	md.mesh.z=[md_global_extruded_3D.mesh.z; md_regional_3D.mesh.z];
	md.mesh.lat=asind(md.mesh.z/r);
	md.mesh.long=atan2d(md.mesh.y,md.mesh.x);
	md.mesh.r=md.mesh.x*0+r;
	md.mesh.numberofvertices=length(md.mesh.x);

	md.mesh.elements=[md_global_extruded_3D.mesh.elements; md_regional_3D.mesh.elements+md_global_extruded_3D.mesh.numberofvertices];
	[md,merged,transfer]=removeduplicatevertices(md,100);
end

%%
function [md,merged,transfer]=removeduplicatevertices(md,mergeres)
    % 'transfer' is the list of indices to refer to the index of regional mesh on the global mesh., of size of md_global_extruded_3D.mesh.numberofverticies+md_regional_3D.mesh.numberofvertices
	% col 1 has values that go up to the size of 'md_global_extruded_3D.mesh.numberofverticies+md_regional_3D.mesh.numberofvertices(converted to a 3D surface)',
	% and the values is larger than in col2 because there are overlapping vertices still. 
	% col 2 has values that go up to the size of 'md_global.mesh.numberofverticies', after the overlapping vertices have been removed in the 3D regional mesh.
	% So, the values on col2 represent index values on md_regional and md_world mesh to refer to corresponding vertice indices on the global mesh.

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
	transfer=[vertices' vertices'];
	del=[];
	for i=1:length(dups)
		vi=dups{i}(1);
		vd=dups{i}(2:end);
		transfer(vd,2)=vi;
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
	
	for i=1:length(transfer)
		transfer(i,2)=find(oldv==transfer(i,2));
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