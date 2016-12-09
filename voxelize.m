function voxelize(filename, scale)
  if nargin < 2
    scale = 1;  
  end
  if nargin < 1
    filename = 'sampleData/arbatpuodis.obj';
  end
  [FC, VR] = fileReader(filename);

  mn = min(VR);
  VR = VR - growVertically(mn, size(VR, 1));
  VR = VR * scale;

  maxCoordinate = max(max(VR));

  volSize = 1;
  while volSize < maxCoordinate
    volSize = volSize * 2;
  end

  figure(1); clf; hold on, grid on, axis equal; h(1)=gca;
  patch('Faces',FC,'Vertices',VR, 'FaceColor',[0,0,1]);
  figure(2); clf; hold on, grid on, axis equal; h(2)=gca;
  linkprop(h, 'view');
  
  tic
  split(VR, FC, 1:size(FC,1), [0,0,0], ones(1,3)*volSize);
  toc
end

function split(VR, FC, fci, p1, p2) 
  size = p2(1) - p1(1);
  if size == 1
    cube_plot(p1, 1, 'r');
  else 
    for i = 0:7
      pd = (p2 - p1) / 2;
      q1 = double(bitget(int8(i), 1:3)) .* pd + p1;
      q2 = q1 + ones(1, 3) * size / 2;
      
      fci1 = [];
      for j = 1:length(fci)
        tr = VR(FC(fci(j),:), :);
        if isIntersecting(tr, q1, q2)
          fci1 = [fci1 fci(j)];
          if size == 2
            break;  
          end
        end
      end
      
      if ~isempty(fci1)
        split(VR, FC, fci1, q1, q2);
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
  pm = (p1 + p2) / 2;
  halfSize = (p2(1) - p1(1)) / 2;
  tr = tr - growVertically(pm, 3);
  
  % test 1
  for i = 1:3
    if (min(tr(:,i)) > halfSize) || (max(tr(:,i)) < -halfSize)
      b = false;
      return;
    end
  end
  
  % test 2
  if ~isIntersectingPlane(tr, p1, p2)
    b = false; 
    return;  
  end
  
  % test 3
  for i = 1:3
    for j = 1:3
      e = zeros(1, 3);
      e(i) = 1;
      f = tr(mod(i, 3)+1,:) - tr(i,:);
      a = cross(e, f);
      
      p = zeros(1, 3);
      r = 0;
      for k = 1:3
        p(k) = dot(a, tr(k,:));
        r = r + halfSize * abs(a(k));
      end
      
      if (min(p) > r) || (max(p) < -r)
        b = false;
        return;
      end
    end
  end
  b = true;
end

function [b] = isIntersectingPlane(tr, p1, p2)
  f1 = tr(2,:) - tr(1,:);
  f2 = tr(3,:) - tr(2,:);
  n = cross(f1, f2);
  r = 0;
  for i=0:7
     bits = double(bitget(int8(i), 1:3));
     pp = bits.*p1 + (1 - bits).*p2;
     if dot(n, pp - tr(1)) > 0
       r = r + 1;  
     end
  end
  b = (abs(r) ~= 7);
end

function [b] = hasPointInsideBox(tr, p1, p2)
% Checks if at least one point of the trangle is inside AABB
  p1 = growVertically(p1, 3);
  p2 = growVertically(p2, 3);
  
  isCoordinateBetween = ((tr >= p1) + (tr < p2)) > 1;
  isPointBetween = sum(isCoordinateBetween, 2) > 2;
  hasPointInside = sum(isPointBetween) > 0;
  
  b = hasPointInside;
end