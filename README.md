# Computer-Vision-Edge-detection
This program is for finding edges in a picture.

First creates many morlet wavelets using different thetas and sigmas. Each wavelets has a real part
and an imaginary part.
Then the picture we are trying to find the edges in, which in this case is the picture of the circle,
is convolved with a Gaussian Blur. 
Then the blurred images are convolved once again with 
