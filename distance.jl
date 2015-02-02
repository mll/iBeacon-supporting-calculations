# Marek Lipert 2015

type Point
  x::Float64
  y::Float64
  visited::Int64
  cluster::Int64
end


type Rectangle
  parent::Int64
  x0::Float64
  x1::Float64
  y0::Float64
  y1::Float64
  children::Array{Int64}
  kind::Int64 # 1 - leaf, 0 - non-leaf
end

type Cluster
  points::Array{Int64}
end



distance(x::(Float64,Float64),y::(Float64,Float64)) = sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
distance(x,y) = sqrt((x.x-y.x)^2 + (x.y-y.y)^2)
distance(x,y::(Float64,Float64)) = sqrt((x.x-y[1])^2 + (x.y-y[2])^2)

distance(x::(Float64,Float64),y) = sqrt((x[1]-y.x)^2 + (x[2]-y.y)^2)

mean(sequence) = reduce(+,sequence)/length(sequence)
sigma(sequence,mean) = sqrt(reduce(+,map((x) -> (x-mean)^2,sequence))/(length(sequence)-1))


# helpers 

# solves quadratic equation for two intersecting circles

function compute_possible_points(r1,r2,p1,p2)
  x1 = p1[1]
  x2 = p2[1]
  y1 = p1[2]
  y2 = p2[2]

  # safety (epsilon cannot be too small for numerical stability)
  epsilon = 0.001
  if abs(x1-x2) < epsilon 
    x1 = x1 + epsilon
  end
  if abs(y1-y2) < epsilon 
    y1 = y1 + epsilon 
  end
  @assert distance(p1,p2) < r1 + r2 "Circles do not cross"
 
  a = (x2-x1)/(y2-y1)
  A = 1 + a^2
  
  S = ((r1^2 - r2^2) - ((x1^2-x2^2) + (y1^2-y2^2)))/(2*(y2-y1))
  B = 2*a*(y1 - S) - 2*x1
  C = (S - y1)^2 + x1^2 - r1^2
  
  delta = B^2 - 4 * A * C
  @assert delta > 0 
  
  x_temp = (-B - sqrt(delta))/(2*A)
  sol_1 = (x_temp , S - a * x_temp)
  
  x_temp = (-B + sqrt(delta))/(2*A) 
  sol_2 = (x_temp , S - a * x_temp)
  (sol_1, sol_2)
end

# checks if a circle of radius eps and center at point intersects the rectangle

function intersects(rectangle,point,eps)
  x0 = point.x-eps
  x1 = point.x+eps
  y0 = point.y-eps
  y1 = point.y+eps
  if distance(point,Point(x0,y0,0,-1)) > eps && distance(point,Point(x0,y1,0,-1)) > eps && distance(point,Point(x1,y0,0,-1)) > eps && distance(point,Point(x1,y1,0,-1)) > eps
    if point.x > x0 && point.x < x1 && point.y > y0 && point.y < y1
      return 1
    else
      return 0
    end
  else
    return 1
  end
end

# searches for all points that are inside circle of 'radius' radius centered on point 

function find_neighbours_inside_circle(rectangles,points,point,radius,rectangle)
    ret_val = Int64[]
    if rectangle.kind == 1
      for i in rectangle.children
        if distance(point,points[i]) < radius
          push!(ret_val,i)
        end
      end
      return ret_val
    end
    
    for i in rectangle.children
      child = rectangles[i]
      if intersects(child,point,eps) == 1
        append!(ret_val,find_eps(rectangles,points,point,radius,child))
      end
    end
    ret_val
end

function DBSCAN(D,eps,MinPts,rectangles)
   clusters = Cluster[]
   for i = 1:length(D)
     point = D[i]
     if point.visited > 0
       continue
     end 
     point.visited = 1
     NeighborPoints = find_neighbours_inside_circle(rectangles,D,point,eps,rectangles[end])
     if length(NeighborPoints) < MinPts
       point.visited = 2
     else
       cl = Cluster(Int64[])
       push!(clusters,cl)
       push!(cl.points,i)
       point.cluster = length(clusters)
       while length(NeighborPoints) > 0
         point2 = D[NeighborPoints[1]]
         if point2.visited != 1 
           point2.visited = 1
           NeighborPoints2 = find_neighbours_inside_circle(rectangles,D,point2,eps,rectangles[end])
           if length(NeighborPoints2) >= MinPts
             NeighborPoints = append!(NeighborPoints,NeighborPoints2)
           end
         end
         if point2.cluster == -1
           push!(cl.points,NeighborPoints[1])
           point2.cluster = length(clusters)
         end
         shift!(NeighborPoints)
       end
    end
  end  
  clusters
end


function build_Rtree(input_points,leaf_rectangle_size)
  
  # Nearest - X
  dataset = sort(input_points, by = (x) -> x.x)
  rectangles = Rectangle[]
  
  # building R-tree
  for ci = 1:ceil(length(dataset)/leaf_rectangle_size)
    points = []
    e = (ci == ceil(length(dataset)/leaf_rectangle_size)) ? length(dataset) : leaf_rectangle_size*ci 
    points = dataset[leaf_rectangle_size*(ci-1)+1:e] 
    
    x0 = points[1].x
    x1 = points[end].x
    y0 = reduce( (x,y) -> x.y < y.y ? x : y,points).y
    y1 = reduce( (x,y) -> x.y > y.y ? x : y,points).y
    push!(rectangles,Rectangle( -1, x0,y0,x1,y1,[(ci-1)*leaf_rectangle_size + 1:e] , 1))
  end
  
  c_begin = 1
  c_end   = length(rectangles)
  while (c_end - c_begin + 1) > 1
    no = ceil((c_end - c_begin + 1) / 2)
    for ci = 1:no
      if (c_end - c_begin + 1) % 2 == 0 || ci!=no
        cl = rectangles[2*ci-1 + c_begin - 1]
        cr = rectangles[2*ci + c_begin - 1]
        x0 = min(cl.x0,cr.x0)
        x1 = max(cl.x1,cr.x1)
        y0 = min(cl.y0,cr.y0)
        y1 = max(cl.y1,cr.y1)
        cl.parent = c_end + ci
        cr.parent = c_end + ci
        push!(rectangles,Rectangle(-1,x0,y0,x1,y1,[int(2*ci-1++ c_begin - 1):int(2*ci + c_begin - 1)],0))
      else
        cl = rectangles[2*ci-1 + c_begin - 1]
        x0 = cl.x0
        x1 = cl.x1
        y0 = cl.y0
        y1 = cl.y1
        cl.parent = c_end + ci  
        push!(rectangles,Rectangle(-1,x0,y0,x1,y1,[int64((2*ci-1 + c_begin - 1))],0))
      end
    end
    c_begin = c_end + 1
    c_end = c_end + no 
  end
  dataset,rectangles
end


# O(N log(N))

function smart_clusterizer(candidates)
  @assert length(candidates) > 1
  #             x        y    visited
  unpacked = Point[]

  for i = 1:length(candidates)
    push!(unpacked,Point(candidates[i][1][1],candidates[i][1][2],0,-1))
    push!(unpacked,Point(candidates[i][2][1],candidates[i][2][2],0,-1))
  end

  dataset , rectangles = build_Rtree(unpacked,2)
  
  clusters = DBSCAN(dataset,0.3,2,rectangles)
  
  @assert length(clusters) > 0
  
  max_cluster = clusters[1]
  for cluster in clusters
    if length(cluster.points) > length(max_cluster.points)
      max_cluster = cluster
    end
  end
  
  all_x = map( (p) -> dataset[p].x ,max_cluster.points) 
  all_y = map( (p) -> dataset[p].y ,max_cluster.points) 
  x = mean(all_x)
  y = mean(all_y)
  dx = sigma(all_x,x)
  dy = sigma(all_y,y)
  (x,y, dx, dy)
end

# Naive but simple
#
# Returns tuple (x,y,dx,dy)
#
# We do not control the error in this computation, hence -1

function compute_position_naive(ranges, beacon_positions)
  @assert length(ranges) == length(beacon_positions)
  @assert length(ranges) == 3 "We support only 3 beacons for now"
  bp = beacon_positions
  A = [ 
        2.0 * (bp[2][1] - bp[1][1])   2.0 * (bp[2][2] - bp[1][2]) ;
        2.0 * (bp[3][1] - bp[1][1])   2.0 * (bp[3][2] - bp[1][2])
      ]
  B = [ 
        (ranges[1]^2 - ranges[2]^2) - ((bp[1][1]^2-bp[2][1]^2) + (bp[1][2]^2-bp[2][2]^2)) ; 
        (ranges[1]^2 - ranges[3]^2) - ((bp[1][1]^2-bp[3][1]^2) + (bp[1][2]^2-bp[3][2]^2)) 
      ]
  sol = inv(A) * B
  (sol[1],sol[2],-1,-1)
end


# More robust position computation. Many beacons (>= 3) can be used to improve accuracy
#
# Returns tuple (x,y,dx,dy)
#

function compute_position(ranges, beacon_positions)
  @assert length(ranges) == length(beacon_positions)
  @assert length(ranges) >= 3 "At least three points are needed for unambigous position computation"
  
  candidates = ((Float64,Float64),(Float64,Float64))[]
  bp = beacon_positions
  for i = 1:length(bp)
   for j = 1:length(bp)
     if i >=j  
       continue
     end
     push!(candidates,compute_possible_points(ranges[i],ranges[j],bp[i],bp[j]))
   end
  end
  smart_clusterizer(candidates)
end



r = [1.0,sqrt(2.0),1.0]
p = [(1.0,0.0),(0.0,0.0),(0.0,1.0)]
  
println("Naive position: ")
println(compute_position_naive(r,p))
println("Good position: ")
println(compute_position(r,p))

println("Now with distortion")
eps = 0.1

r = [1+eps,sqrt(2)-eps,1]
p = [(1.0,0.0),(0.0,0.0),(0.0,1.0)]
  
println("Naive position: ")
println(compute_position_naive(r,p))
println("Good position: ")
println(compute_position(r,p))










