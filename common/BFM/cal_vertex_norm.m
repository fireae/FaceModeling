function normal = cal_vertex_norm(shape, faces, index)
 faces1 = faces(:,1);
 triangles =  find(faces1 == index);
 tmp = cross(shape(:,faces(triangles,3))-shape(:,faces(triangles,1)),shape(:,faces(triangles,2))-shape(:,faces(triangles,1)));
 normal = sum(tmp,2);
 
 faces2 = faces(:,2);
 triangles =  find(faces2 == index);
 tmp = cross(shape(:,faces(triangles,1))-shape(:,faces(triangles,2)),shape(:,faces(triangles,3))-shape(:,faces(triangles,2)));
 normal = normal + sum(tmp,2);
 
 faces3 = faces(:,3);
 triangles =  find(faces3 == index);
 tmp = cross(shape(:,faces(triangles,2))-shape(:,faces(triangles,3)),shape(:,faces(triangles,1))-shape(:,faces(triangles,3)));
 normal = normal + sum(tmp,2);
 
 normal = normal/sqrt(sum(normal.^2));
end