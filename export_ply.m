%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Description..... : 	Save the reconstruction PLY format, in 
%						mesh.obj, mesh.mtl and mesh.png
%						INPUT : XYZ, N, RHO -- nrows x ncols x 3
%						
%	Auteur ......... : 	Fotios Logothetis (adapted from Yvain Queau) 
%
%	Date de création : 	06/10/2014
%	Date de modif... :  06/01/2014 par Yvain (moins de points, aussi 
% 						precis)
% connectivity_info is either a N x 3 vertices list or a MxN segmentation
% mask. in the later case, full connectivity is assumed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export_ply(XYZ,connectivity_info,filename, rho_rgb)
	    
    export_rho=0;
    
     if(size(connectivity_info,2)==3 )
        disp('No meshing, just triangle input mode...')
        vertices=XYZ';
        [nverts,pts_per_vert]=size(vertices);
        assert(pts_per_vert==3);
        face_vertices=connectivity_info;
        [nfaces,~]=size(face_vertices);
        
         if (size(rho_rgb,1)>0) 
            export_rho=1;
            rho_r=rho_rgb(:,1);
            rho_g=rho_rgb(:,2);
            rho_b=rho_rgb(:,3);
            
            rho_r=uint8(255*rho_r);
            rho_g=uint8(255*rho_g);
            rho_b=uint8(255*rho_b);
            rho_rgb=double(cat(2,rho_r,rho_g,rho_b));
            %need a combined vector of everything to be able to write in
            %one command
            comb=[vertices,rho_rgb];
         end
     else 
        %need to perform meshing 
        disp('Meshing ...')
        mask=connectivity_info;
        
        if (size(rho_rgb,1)>0) 
            export_rho=1;
            rho_r=imresize(rho_rgb(:,:,1),size(mask));
            rho_g=imresize(rho_rgb(:,:,2),size(mask));
            rho_b=imresize(rho_rgb(:,:,3),size(mask));
        end 
        [nrows,ncols] = size(mask);
    %% SWITCH TO OPENGL AXIS- i.e. negative z towards incide the screen        
        X = XYZ(:,:,1);
        Y = -XYZ(:,:,2);
        Z = -XYZ(:,:,3);

        indices_mask = find(mask>0);
        [~,~]=ind2sub(size(mask),indices_mask);
        indices = zeros(size(mask));
        indices(indices_mask) = 1:length(indices_mask);		
        mask=[mask;zeros(1,size(mask,2))];
        mask=[mask,zeros(size(mask,1),1)];
       
        vertices = [X(indices_mask),Y(indices_mask),Z(indices_mask)];
        if export_rho
            rho_r=uint8(255*rho_r(indices_mask));
            rho_g=uint8(255*rho_g(indices_mask));
            rho_b=uint8(255*rho_b(indices_mask));
            rho_rgb=double(cat(2,rho_r,rho_g,rho_b));
            %need a combined vector of everything to be able to write in
            %one command
            comb=[vertices,rho_rgb];
        end
        
        indices_lower_triangle = find(mask(1:end-1,1:end-1)>0 & mask(2:end,2:end) & mask(2:end,1:end-1)); 
        [I_lt,J_lt] = ind2sub([nrows ncols],indices_lower_triangle);	
        indices_bas = sub2ind([nrows ncols],I_lt+1,J_lt);
        indices_bas_droite = sub2ind([nrows ncols],I_lt+1,J_lt+1);
        face_vertices = [indices(indices_lower_triangle),indices(indices_bas),indices(indices_bas_droite)];
	
        indices_upper_triangle = find(mask(1:end-1,1:end-1)>0 & mask(2:end,2:end) & mask(1:end-1,2:end));
        [I_ut,J_ut] = ind2sub([nrows ncols],indices_upper_triangle); 
        indices_droite = sub2ind([nrows ncols],I_ut,J_ut+1);
        indices_bas_droite = sub2ind([nrows ncols],I_ut+1,J_ut+1);
        face_vertices = [face_vertices;...
					 indices(indices_upper_triangle),indices(indices_bas_droite),indices(indices_droite)];	
    
        face_vertices=face_vertices-1;
                 
        [nverts,~]=size(vertices);
        [nfaces,~]=size(face_vertices);
     end 
    
    fprintf(1, 'writing %s ... \n',filename);
    %% NOTE the uppercase W for not flushing all the time  
    % TODO binary https://uk.mathworks.com/help/matlab/ref/fopen.html
	fileID = fopen(filename,'W'); 
    
    fprintf(fileID,'ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\n',nverts); %
    if export_rho
        fprintf(fileID,'property uchar red\nproperty uchar green\nproperty uchar blue\n');
    end
    fprintf(fileID,'element face %d\nproperty list uint8 int32 vertex_indices\nend_header\n', nfaces);   
    
    %single command is much faster but a bit tricky to syntax properly http://ocw.uci.edu/upload/files/mae10_w2011_lecture08.pdf
    if export_rho
        fprintf(fileID,'%f %f %f %d %d %d\n',comb');
    else
        fprintf(fileID,'%f %f %f\n',vertices');
    end
    fprintf(fileID,'3 %d %d %d\n',face_vertices');
    
%     for ii=1:nverts
%         fprintf(fileID,'%f %f %f',vertices(ii,1),vertices(ii,2),vertices(ii,3));
%         if export_rho
%             fprintf(fileID,' %d %d %d',rho_r(ii), rho_g(ii),rho_b(ii));
%         end
%         fprintf(fileID,'\n');
%     end
    
%      for ii=1:nfaces
%         fprintf(fileID,'3 %d %d %d\n',face_vertices(ii,1),face_vertices(ii,2),face_vertices(ii,3));
%     end   
    fclose(fileID); 
end


