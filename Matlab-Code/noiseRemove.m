function yClean=noiseRemove(y,n)
mtlb_noisy = y;
noise = n;
% Define Adaptive Filter Parameters
filterLength = 32;
weights = zeros(1,filterLength);
step_size = 0.004;
% Initialize Filter's Operational inputs
output = zeros(1,length(mtlb_noisy));
err = zeros(1,length(mtlb_noisy));
input = zeros(1,filterLength);
% For Loop to run through the data and filter out noise
for n = 1: length(mtlb_noisy),
      %Get input vector to filter
      for k= 1:filterLength
          if ((n-k)>0)
              input(k) = noise(n-k+1);
          end
      end
      output(n) = weights * input';  %Output of Adaptive Filter
      err(n)  = mtlb_noisy(n) - output(n); %Error Computation
      weights = weights + step_size * err(n) * input; %Weights Updating 
  end
yClean = err;
end