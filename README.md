# Computer-Vision-Edge-detection
This program is for finding edges in a picture.

First I create many morlet wavelets using different thetas and sigmas. Each wavelets has a real part
and an imaginary part.
Then I convolve the picture I am trying to find the edges in, which in this case is a picture of a circle, 
with the wavelets.
Then I make one matrix for real and imaginary each by finding the max value from all the convolved
matrixes in each corresponding pixel.
Then I make a clear-edge picture by finding the difference, either by ratio or subtracting.

You can find the documentation in more detail in Homework1.pdf
