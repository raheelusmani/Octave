clear all;
close all;
clc;

abc4_sqr
mu=1;
f= 0;
nNode = size(msh.POS,1);
nElem = size(msh.QUADS,1);

soln  = zeros(nNode,1);
ndof=1
totaldof = nNode* ndof
bound_nodes = msh.LINES(:,1:2)
bound_nodes = bound_nodes(:)
bound_nodes = unique(bound_nodes)

dofs_fixed = bound_nodes;
dofs_free = setdiff(1:totaldof, dofs_fixed);
elem_node_conn = msh.QUADS(:,1:4);
elem_dof_conn = elem_node_conn;

Kglobal =zeros(totaldof,totaldof)
Fglobal = zeros(totaldof,1)
soln = zeros(totaldof,1)
soln_specified = zeros(totaldof,1);

for ii = 1:size(dofs_fixed,1)
   dof = dofs_fixed(ii);
   xc = msh.POS(dof,1);
   yc = msh.POS(dof,2);
   
   soln_specified(dof) = (cosh(pi*yc)-coth(pi)*sinh(pi*yc))*sin(pi*xc);
end
  
for elnum = 1:nElem
  
  [Klocal,Flocal]= fournoded(msh.POS, elem_node_conn(elnum,:), mu, f)
  
  rows = elem_dof_conn (elnum,:);
  cols = rows; 
  
  Kglobal(rows,cols)= Kglobal(rows,cols) + Klocal;
  Fglobal(rows,1)   =  Fglobal(rows,1)   + Flocal;
  
end

for ii= 1:size(dofs_fixed,1)
   dof = dofs_fixed(ii);
  
  Fglobal = Fglobal - Kglobal(:,dof)*soln_specified(dof);
  soln(dof) = soln_specified(dof);
end  
Kglobal
Fglobal
soln

soln(dofs_free) = Kglobal (dofs_free,dofs_free)\ Fglobal(dofs_free);

    trisurf(msh.QUADS(:,1:4), msh.POS(:,1) , msh.POS(:,2) , soln , 'edgecolor' , 'none' , 'facecolor' , 'interp');
  colormap jet;
  view(2); colorbar;

  