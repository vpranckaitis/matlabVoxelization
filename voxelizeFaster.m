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
    for i = 1:8    
      [q1, q2, FCi1] = getFacesIntersectingOctant(VR, FC, FCi, p1, p2, i);
      
      if ~isempty(FCi1) % if octant has intersecting faces
        split(VR, FC, FCi1, q1, q2);
      end
    end
  end
end

function [b] = isIntersecting(tr, p1, p2)
  if hasPointInsideBox(tr, p1, p2)
    b = true;
    return;
  end

% https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tribox.pdf
% http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

  halfSize = (p2(1) - p1(1)) / 2;
  
  b = false;
  
  % test 1: bounding boxes of triangle and AABB
  [mn, mx] = minmax3(tr(:,1));
  if (mn > halfSize) || (mx < -halfSize), return; end
  [mn, mx] = minmax3(tr(:,2));
  if (mn > halfSize) || (mx < -halfSize), return; end
  [mn, mx] = minmax3(tr(:,3));
  if (mn > halfSize) || (mx < -halfSize), return; end
  
  % test 2: intersection of triangle's plane and AABB
  if ~isBoxIntersectingPlane(tr, p1, p2)
    return;  
  end
  
  % test 3: 9 other tests
  p = zeros(3, 1);

  f = tr(2,:)-tr(1,:);
  ff = abs(f);
  p(1) = -f(3)*tr(1,2) + f(2)*tr(1,3);
%   p(2) = -f(3)*tr(2,2) + f(2)*tr(2,3);
  p(3) = -f(3)*tr(3,2) + f(2)*tr(3,3);
  [mn, mx] = minmax2(p(1), p(3));
  r = halfSize * (ff(3) + ff(2));
  if mn > r || mx < -r, return; end
  
  p(1) = f(3)*tr(1,1) - f(1)*tr(1,3);
%   p(2) = f(3)*tr(2,1) - f(1)*tr(2,3);
  p(3) = f(3)*tr(3,1) - f(1)*tr(3,3);
  [mn, mx] = minmax2(p(1), p(3));
  r = halfSize * (ff(3) + ff(1));
  if mn > r || mx < -r, return; end
  
  p(1) = -f(2)*tr(1,1) + f(1)*tr(1,2);
%   p(2) = -f(2)*tr(2,1) + f(1)*tr(2,2);
  p(3) = -f(2)*tr(3,1) + f(1)*tr(3,2);
  [mn, mx] = minmax2(p(1), p(3));
  r = halfSize * (ff(2) + ff(1));
  if mn > r || mx < -r, return; end
  
  f = tr(3,:)-tr(2,:);
  ff = abs(f);
  p(1) = -f(3)*tr(1,2) + f(2)*tr(1,3);
  p(2) = -f(3)*tr(2,2) + f(2)*tr(2,3);
%   p(3) = -f(3)*tr(3,2) + f(2)*tr(3,3);
  [mn, mx] = minmax2(p(1), p(2));
  r = halfSize * (ff(3) + ff(2));
  if mn > r || mx < -r, return; end

  p(1) = f(3)*tr(1,1) - f(1)*tr(1,3);
  p(2) = f(3)*tr(2,1) - f(1)*tr(2,3);
%   p(3) = f(3)*tr(3,1) - f(1)*tr(3,3);
  [mn, mx] = minmax2(p(1), p(2));
  r = halfSize * (ff(3) + ff(1));
  if mn > r || mx < -r, return; end
  
  p(1) = -f(2)*tr(1,1) + f(1)*tr(1,2);
  p(2) = -f(2)*tr(2,1) + f(1)*tr(2,2);
%   p(3) = -f(2)*tr(3,1) + f(1)*tr(3,2);
  [mn, mx] = minmax2(p(1), p(2));
  r = halfSize * (ff(2) + ff(1));
  if mn > r || mx < -r, return; end
  
  f = tr(1,:)-tr(3,:);
  ff = abs(f);
%   p(1) = -f(3)*tr(1,2) + f(2)*tr(1,3);
  p(2) = -f(3)*tr(2,2) + f(2)*tr(2,3);
  p(3) = -f(3)*tr(3,2) + f(2)*tr(3,3);
  [mn, mx] = minmax2(p(2), p(3));
  r = halfSize * (ff(3) + ff(2));
  if mn > r || mx < -r, return; end
  
%   p(1) = f(3)*tr(1,1) - f(1)*tr(1,3);
  p(2) = f(3)*tr(2,1) - f(1)*tr(2,3);
  p(3) = f(3)*tr(3,1) - f(1)*tr(3,3);
  [mn, mx] = minmax2(p(2), p(3));
  r = halfSize * (ff(3) + ff(1));
  if mn > r || mx < -r, return; end
  
%   p(1) = -f(2)*tr(1,1) + f(1)*tr(1,2);
  p(2) = -f(2)*tr(2,1) + f(1)*tr(2,2);
  p(3) = -f(2)*tr(3,1) + f(1)*tr(3,2);
  [mn, mx] = minmax2(p(2), p(3));
  r = halfSize * (ff(2) + ff(1));
  if mn > r || mx < -r, return; end
  
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
  
  z = [0 0 0;
       0 0 1;
       0 1 0;
       0 1 1;
       1 0 0;
       1 0 1;
       1 1 0;
       1 1 1];
  
  r = 0;
  for i=1:8
     bits = z(i,:);
     v = bits.*p1 + (1 - bits).*p2 - q;  % vertex of the box
     if dot3(n, v) > 0 % if projection is positive 
       r = r + 1;  
     end
     if r > 0 && r < i
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

function [mn, mx] = minmax2(a, b) 
  if a < b
    mn = a; mx = b;
  else
    mn = b; mx = a;
  end
end

function [mn, mx] = minmax3(a)
  mn = a(1); mx = a(1);
  if a(2) < mn, mn = a(2); end
  if a(2) > mx, mx = a(2); end
  if a(3) < mn, mn = a(3); end
  if a(3) > mx, mx = a(3); end
end

function [b] = hasPointInsideBox(tr, p1, p2)
% Checks if at least one point of the trangle is inside AABB

  p1 = [p1; p1; p1];
  p2 = [p2; p2; p2];
  
  b = sum(sum((tr > p1) + (tr < p2), 2) > 5) > 0;
end

function [q1, q2] = getOctant(p1, p2, i) 
% Transforms an integer from 0 to 7 into one of 8 octects of box (p1, p2)
% Binary representation of 'i' means which half of the threen coordinates
% x, y and z to take. For example, binary(6)=110 means to take first half
% for coordinate x, but second half for coordinates y and z.

% https://en.wikipedia.org/wiki/Octant_(solid_geometry)

  z = [0 0 0;
       0 0 1;
       0 1 0;
       0 1 1;
       1 0 0;
       1 0 1;
       1 1 0;
       1 1 1];
  halfsize = (p2(1) - p1(1)) / 2;
  q1 = z(i,:) * halfsize + p1;
  q2 = q1 + halfsize;
end

function [q1, q2, FCi1] = getFacesIntersectingOctant(VR, FC, FCi, p1, p2, i)
% Finds indexes of faces from 'FCi' intersecting octant 'i' of the
% AABB (p1, p2). 

  size = p2(1) - p1(1);

  [q1, q2] = getOctant(p1, p2, i);
  
  % move objects so box's center would match origin
  eps = 1e-2;
  qm = (q1 + q2) / 2;
  qt1 = q1 - qm + eps;
  qt2 = q2 - qm - eps;
  
  FCi1 = zeros(length(FCi),1); % indexes of faces intersecting current octant
  n = 0;
  for j = 1:length(FCi)
    tr = VR(FC(FCi(j),:), :);
    tr(:,1) = tr(:,1) - qm(1);
    tr(:,2) = tr(:,2) - qm(2);
    tr(:,3) = tr(:,3) - qm(3);
    
    if isIntersecting(tr, qt1, qt2)
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