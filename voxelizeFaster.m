function voxelizeFaster(filename, scale)
  if nargin < 2
    scale = 1;  
  end
  if nargin < 1
    filename = 'sampleData/arbatpuodis.obj';
  end
  [FC, VR] = fileReader(filename);

  mn = min(VR);
  VR = VR - growVertically(mn, size(VR, 1));
  VR = VR * scale; % relative voxel size is controlled by scaling 3D model

  maxCoordinate = max(max(VR));

  k = 1; % size of the edge of bounding box, power of two
  while k < maxCoordinate
    k = k * 2;
  end

  figure(1); clf; hold on, grid on, axis equal; h(1)=gca;
  patch('Faces',FC,'Vertices',VR, 'FaceColor',[0,0,1]);
  figure(2); clf; hold on, grid on, axis equal; h(2)=gca;
  linkprop(h, 'view');
  
  tic
  split(VR, FC, 1:size(FC,1), [0,0,0], [k k k]);
  toc
  
%   pause;
end

function split(VR, FC, FCi, p1, p2)
% Recursively splits AABB into 8 octants and finds intersecting faces

  size = p2(1) - p1(1);
  
  if size == 1 % 1 is a minimum size for voxel
    cube_plot(p1, 1, 'r');
  else 
    for i = 0:7    
      [q1, q2, FCi1] = getFacesIntersectingOctant(VR, FC, FCi, p1, p2, i);
      
      if ~isempty(FCi1) % if octant has intersecting faces
        split(VR, FC, FCi1, q1, q2);
      end
    end
  end
end

function [b] = isIntersecting(tr, p1, p2)
  eps = 1e-2;
  
  % move objects so box's center would match origin
  pm = (p1 + p2) / 2;
  p1 = p1 - pm + eps;
  p2 = p2 - pm - eps;
  tr(:,1) = tr(:,1) - pm(1);
  tr(:,2) = tr(:,2) - pm(2);
  tr(:,3) = tr(:,3) - pm(3);
  
  if hasPointInsideBox(tr, p1, p2)
    b = true;
    return;
  end

% https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tribox.pdf
% http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

  halfSize = (p2(1) - p1(1)) / 2;
  
  b = false;
  
  % test 1: bounding boxes of triangle and AABB
  if (min(tr(:,1)) > halfSize) || (max(tr(:,1)) < -halfSize) || ...
     (min(tr(:,2)) > halfSize) || (max(tr(:,2)) < -halfSize) || ...
     (min(tr(:,3)) > halfSize) || (max(tr(:,3)) < -halfSize)
    return;
  end
  
  % test 2: intersection of triangle's plane and AABB
  if ~isBoxIntersectingPlane(tr, p1, p2)
    return;  
  end
  
  % test 3: 9 other tests
  for i = 1:3
    for j = 1:3
      e = zeros(1, 3);
      e(i) = 1;
      f = tr(mod(j, 3)+1,:) - tr(j,:);
      a = cross3(e, f);
      
      p = zeros(1, 3);
      for k = 1:3
        p(k) = dot3(a, tr(k,:));
      end
      
      r = sum(halfSize * abs(a));
      
      if (min(p) > r) || (max(p) < -r)
        return;
      end
    end
  end
  
  b = true; % all tests passed
end

function [b] = isBoxIntersectingPlane(tr, p1, p2)
% Gets every vertices' distance from the plane of triangle. The distance 
% can be positive or negative, depending of the size of the halfspace 
% the vertex is in. If all the vertices are on the same side, box doesn't 
% intersect triangle's plane.

  q = tr(1,:);
  f1 = tr(2,:) - q;
  f2 = tr(3,:) - q;
  n = cross3(f1, f2);
  
  r = 0;
  for i=0:7
     bits = double(bitget(int8(i), 1:3));
     v = bits.*p1 + (1 - bits).*p2;  % vertex of the box
     v = v - q; % vector from point 'q' to the vertex 'v'
     if dot3(n, v) > 0 % if projection is positive 
       r = r + 1;  
     end
     if r > 0 && r < i + 1
       break;
     end
  end
  
  b = (r ~= 8);
end

function [n] = cross3(a, b)
  n = zeros(size(a));
  n(1) = a(2)*b(3)-a(3)*b(2);
  n(2) = a(3)*b(1)-a(1)*b(3);
  n(3) = a(1)*b(2)-a(2)*b(1);
end

function [s] = dot3(a, b)
  s = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

function [b] = hasPointInsideBox(tr, p1, p2)
% Checks if at least one point of the trangle is inside AABB

  p1 = growVertically(p1, 3);
  p2 = growVertically(p2, 3);
  
  b = sum(sum((tr > p1) + (tr < p2), 2) > 5) > 0;
end

function [q1, q2] = getOctant(p1, p2, i) 
% Transforms an integer from 0 to 7 into one of 8 octects of box (p1, p2)
% Binary representation of 'i' means which half of the threen coordinates
% x, y and z to take. For example, binary(6)=110 means to take first half
% for coordinate x, but second half for coordinates y and z.

% https://en.wikipedia.org/wiki/Octant_(solid_geometry)

  halfsize = (p2(1) - p1(1)) / 2;
  q1 = double(bitget(int8(i), 1:3)) * halfsize + p1;
  q2 = q1 + halfsize;
end

function [q1, q2, FCi1] = getFacesIntersectingOctant(VR, FC, FCi, p1, p2, i)
% Finds indexes of faces from 'FCi' intersecting octant 'i' of the
% AABB (p1, p2). 

  size = p2(1) - p1(1);

  [q1, q2] = getOctant(p1, p2, i);
  
  FCi1 = zeros(length(FCi),1); % indexes of faces intersecting current octant
  n = 0;
  for j = 1:length(FCi)
    tr = VR(FC(FCi(j),:), :);

    if isIntersecting(tr, q1, q2)
      n = n + 1;
      FCi1(n) = FCi(j);

      % its enough to find a single face intersecting an octant of size 1
      if size == 2 
        break;  
      end
    end

  end
  
  FCi1 = FCi1(1:n);
end