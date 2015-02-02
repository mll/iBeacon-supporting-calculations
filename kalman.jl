#  Marek Lipert 2015

mean(sequence) = reduce(+,sequence)/length(sequence)
sigma(sequence,mean) = sqrt(reduce(+,map((x) -> (x-mean)^2,sequence))/(length(sequence)-1))

# filters

function kalman(sequence)
  x0 = mean(sequence)
  s = 3.0                # measurement error
  e = 1000    
  r = 0.0
  for i in sequence
    p = e/(e+s)
    r = r + p * (i - r)
    e = (1 - p)*e
  end
  r
end

function motion_kalman(sequence)
  dt = 1.0      
  sigma_process = 0.3 # process error
  r = 3.0    # measurement error
  L = 1000
  F = [1 dt; 0 1]
  Q = [dt^3 / 3.0 dt^2 / 2; dt^2 / 2 dt ] .* sigma_process
  H = [1 0]
  R = [r]
  
  x0 = [0.0 ; 0.0]
  P0 = [L 0; 0 L]
  x1 = x0
  P1 = P0
  U = [1 0 ; 0 1]
  for z in sequence[2:]
    # predict
    x1 = F * x1             
    P1 = F * P1 * F' + Q    
    # update
    y = [z] - H * x1        
    S = H * P1 * H' + R     
    K = P1 * H'  ./ S[1,1]   
    x1 = x1 + K .* y[1]      
    P1 = (U - K*H) * P1     
  end
  (H * x1)[1] 
end
  
s0 = [25 7 14 8 13 20 ]
s1 = [3.2 3.1 3.2 3.3 3.5]
s2 = [3.2 4.2 2.8 6.1 2.7]

measures = (s0,s1,s2)

for l in measures
  print(l)
  @printf("mean %.2f sd %.2f kalman %.2f motion_kalman %.2f\n",mean(l),sigma(l,mean(l)),kalman(l),motion_kalman(l))
  println("------------------")
end



