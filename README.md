# Lenia
Flow lenia with reintegration tracking and advection. But in 3D. And on Windows XP.

![hippo](https://media2.giphy.com/media/v1.Y2lkPTc5MGI3NjExbHc1bzl0cXFjbXl2NXBsODM3dm9rYWhydG5vZmR0ajFjOXJ5cWkwOSZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/4q1vuFa4fggxx7Nk7Z/giphy.gif)

Ok... But why?

About a year and a half ago I implemented Conway's game of life in C to get familiar with SDL2 and OpenGL, which would later lead to my experimentations in game engine development. It took a few hours and I had the classic cell grid at a massive scale. It was really fun to watch. Then I got a wonderful ![video](https://www.youtube.com/watch?v=6kiBYjvyojQ) in my recommended detailing lenia, a recent generalization of conway. This led to a project that I would spend the next month or so on. I implemented lenia in 2D OpenGL, with a neat menu where you could tweak all of the sim parameters. When I learned about flow lenia, I implemented that as well, since I thought it looked so organic and beautiful. I always knew I could push it further though, I just didn't have the tools yet. Fast forward a... while and I ended up with Windows XP on this old computer from the garage. I brushed the dust off and started programming on it. It is really delightful to retreat from society, the internet, whatever, and this machine further satisfies that. If I want to get code onto it, I need a flash drive. No internet, nothing. With all this stuff going on nowadays, I prefer it that way. Anyway, as someone who loves challenges, it came across my mind to implement lenia in 3D on this thing. Why not? So there you go, that's why. I just think it's fucking cool, man.

So, what's inside of it? A marching cubes isosurface representation of four panel flow lenia. Due to the heavy computations involved in flow lenia, translating this to 3D required a lot of optimizations. Way too many to list here, mainly vectorization and padding techniques. I will discuss this in more detail further down. I also wrote some comments at the top of lenia.cpp that might enlighten you. Also, goes without saying, yes it's "c++" but as you know, this is really just written in C.

Here's how we calculate a simulation step:
- T = Perform forward FFT on GRID per channel
- Temp = Per kernel, perform complex multiplication with the kernel's input channel T and the kernel forward FFT, KG
- U (Weighted Neighbor Sums) = Reverse FFT Temp
- H = for each U, compute growth using sigmoid bell function based on kernel m and s parameters
- Here is our first divergence. Optionally apply growth to GRID and get BACK.
- Second divergence. We have two methods to calculate gradients. Sobel is a bit biased in 3D, but it works for certain params. Central difference is a bit more uniform and boring. To help with smoothing it, there is an optional blurring pass that occurs here. This means another forward fft, complex mul, and reverse fft to convolve. To calculate the sobel gradient, we convolve. To calculate the central difference gradient, we simply subract over axes. Either way, we get nabla a, the gradient of the grid sum, and nabla u, the gradient of the growth.
- Calculate F (force). This is nabla_u * 1-alpha - nabla_a * alpha, where alpha is traditionally GRID **2
- Calculate MU (displacement). For each cell, use F to find out the position for the current timestep.
- Calculate V (velocity) = V + F - yV
- Reintegration. Accumulate mass into TT. Optionally use BACK instead of GRID. We have a lot of diverging here.
  - SRT Mode: Single reintegration. No actual transport, just destination into source. Makes noisy, bloblike patterns.
  - RT Mode: Classical reintegration. For each cell surrounding dst, accumulate into dst based on the difference between dst and src's MU. Apply normalization.
  - Advection Mode: Lagrangian advection. Use pos - V or F to find the position to pull mass from. Accumulate based on neighbors using anisotropy.
- Transposition. After reintegration we have to transpose the data and transform it to end the sim step. This includes slotting TT into GRID and transposing it into a format we can FFT. We also calculate some moving averages for some other sim options.

So with all of that out of the way, what else is done to make this faster? There's a ton of handmade vectorization. It basically makes the code unreadable, but it also makes it incredibly fast. Instead of doing modulus math to wrap around the grid, I pad most grids by a large margin. This margin is 4 so we can stay 16 byte aligned for speed. This means that we need to copy the opposite side of the grid into the borders, or zero the borders depending on the border type that is selected in simulation options. We also have pretty much completely branchless loops. This means many duplicate functions instead of if statements (see the cursed mesh.c). No, the compiler will not optimize this because these varaibles can be changed via ui during runtime. Yes, this was painful. There are 4 worker threads, meaning that a 4 core processor is basically necessary for this. 2 work on mesh generation and 2 work on simulation. When a new sim step is finished, we increase a counter so that the mesh threads never create the same mesh twice. There's probably a bunch I'm missing, feel free to read the code or ask me any questions you have personally.

Why 32 x 32 x 32 and not more? Well, couple of reasons. Literally out of DX9 memory. Thanks DX9. Also, 3D convolution math is heavy. Especially on this hardware and OS, hence the challenge. 36^3 is actually not bad. 40^3 would theoretically be 30fps with certain settings (if we weren't out of memory). You can get away with 44 x 24 x 44. Anyway, you can compile it for yourself and play with the dimensions, you just really need DX9, FFTW, and Windows XP. And of course only dimensions divisible by 4.

That's really it, I plan to have a release up, but it's really a beta binary. Some features are disabled / unresolved, I'm just ready to take a break from it. There's also rare thread conditions that case the display to stop updating on ui interactions. It's probably very fixable, but as I said I need a break. I will also be publishing a video detailing all of this as well. If you want to run it yourself and have any trouble, feel free to reach out to me on discord, I'm always open to talk.

**References**
Marching cubes tables coutesy of https://github.com/nihaljn/marching-cubes/
You can find out more about lenia here https://chakazul.github.io/lenia.html
You can find out more about flow lenia here https://sites.google.com/view/flowlenia/
I want to give a thank you to everyone that helped me out while streaming my development of this live, especially @SpeechH4CK for the help with alignment and simd. And thanks to everyone that gave me a good laugh along the way as well.
