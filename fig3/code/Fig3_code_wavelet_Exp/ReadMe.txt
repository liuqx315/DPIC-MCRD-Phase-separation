Cite this code package as:
Penny, G., Daniels, K.E., Thompson, S.E. (2013). Local pattern properties. Local properties of patterned vegetation: quantifying endogenous and exogenous effects. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences. doi:10.xxxx/rspa.xxxx.xxxx
*******************************************************


***INCLUDED FILES***
LocalPatternProps.m
LPPwindow.m
LPPmerge.m
LocalPattProp_Example.m
ExamplePattern1.tiff
ExamplePattern2.tiff

To run the algorithm to calculate the local pattern properties from a patterned image, the MATLAB function files LocalPatternProps.m, LPPwindow.m, and LPPmerge.m are needed. The MATLAB script LocalPattProp_Example.m is useful for setting the parameters and calling the functions. As prepared, it is ready analyze the example images, ExamplePattern1.tiff and ExamplePattern2.tiff.


***REQUIRED MATLAB TOOLBOXES***
Image processing toolbox
Statistics toolbox


***RUNNING THE EXAMPLE FILE***
LocalPattProp_Example.m is prepared to call and analyze the image ExamplePattern.tiff. Run each block of code sequentially to perform the analysis and view the results.


***FILE DESCRIPTIONS***

- LocalPatternProps.m
This function calculates the local wavelength and direction of a patterned image by successively applying a 2D Fourier algorithm to square sections of the image. The image is divided into a set of overlapping windows, and the wavelength and direction of the pattern in each window are calculated. The results from each window are combined into output arrays.

- LPPwindow.m
This function takes a patterned image as input, crops it to a single window (specified by i,j) and computes the 2D FFT and power spectrum. It then extracts the dominant wavelength and orientation from the power spectrum as outlined in Penny et al (2013). The uniqueness metric is also computed for both the wavelength and orientation. The mean power and max power of the window are also returned.

- LPPmerge.m
This function merges the results from LocalPatternProps.m, which contain overlapping windows. Additionally, this function includes filters to remove windows that do not meet criteria specified by MinPower and MinPatchSize. Once the filter is applied, the windows are merged such that overlapping windows are averaged for each cell. This process is applied to the pattern wavelength, orientation, and uniqueness metrics. The edges of the results are then trimmed by half the width of a window and the results are returned.

- LocalPattProp_Example.m
This script file calls the functions LocalPattProps.m and LPPmerge.m to calculate the local pattern wavelength and orientation (direction) for the example images. The code is prepared so that it can be run as-written to calculate the local properties for the example image.

- ExamplePattern1.tiff
This is an synthetic circular pattern, with wavelength steadily decreasing moving away from the center.

- ExamplePattern2.tiff
This pattern is part of a large patterned area in Fort Stockton, TX, which was analyzed in Penny et al (2013).


***ANALYZING OTHER PATTERNS:***
To analyze a different pattern, create a binarized version of the pattern image and import it into MATLAB. Using LocalPattProp_Example.m as a starting point, adjust the parameters (w, dL, minL, maxL) accordingly. Next, follow the steps in the example code (adjusting MinPower and MinPatchSize in step 4) to run the analysis and view results. See Penny et al (2013) and the comments in LocalPattProp_Example.m for more information.