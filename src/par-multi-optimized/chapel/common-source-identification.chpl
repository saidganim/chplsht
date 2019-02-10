/* Use assertions. */
use Assert;

use CyclicDist;

/* Use sorting. */
use Sort;

/* Allow us to read the contents of a directory. */
use FileSystem;

/* Time the execution. */
use Time;

/* Read in JPG images. */
use read_jpg;

/* Compute PRNU noise patterns. */
use prnu;

/* For FFT */
use FFTW;

/* For Math */
use Math;

/* Configuration parameters */
config const imagedir : string = "images";
config const writeOutput : bool = false;

/* Add a directory to a file name. */
proc addDirectory(fileName : string, dir : string) : string {
  return dir + "/" + fileName;
}

/* Get an array of the file names of the images in directory dir. */
proc getImageFileNames(dir : string) {
    var imageFiles = listdir(dir);
    sort(imageFiles);
    return addDirectory(imageFiles, dir);
}

/* Write a real array to a file. */
proc write2DRealArray(array : [] real, fileName :string) {
  assert(array.rank == 2);
  var file = open(fileName, iomode.cw);
  var channel = file.writer();

  for i in array.domain.dim(1) {
    for j in array.domain.dim(2) {
      channel.writef("%{#####.#####} ", array[i, j]);
    }
    channel.writeln();
  }
}

// proc dotProduct(ref C: [?DC] complex, ref A: [?DA] complex, ref B: [?DB] complex){

  
// }

proc main() {
  /* Obtain the images. */
  var imageFileNames = getImageFileNames(imagedir);

  /* n represents the number of images that have to be correlated. */
  var n = imageFileNames.size;

  /* h, w will represent the height and width of an image or PRNU noise pattern 
   * throughout the code.
   */
  var h, w : int;
  (h, w) = getDimensionsJPG(imageFileNames.front());

  /* Create a domain for the correlation matrix. */
  const corrDomain : domain(2)  dmapped Cyclic(startIdx=(1,1)) = {1..n, 1..n};
  var corrMatrix : [corrDomain] real;

  var overallTimer : Timer;

  const imageDomain: domain(2) = {0..#h,0..#w};
  var images : [imageFileNames.size][imageDomain] RGB; 
  for i in 1..imageFileNames.size do{
    readJPG(images[i], imageFileNames[i]);
  }

  writeln("Running Common Source Identification...");
  writeln("  ", n, " images");
  writeln("  ", numLocales, " locale(s)");

  overallTimer.start();

  /* Below is example code that contains part of the pieces that you need to 
   * compute the correlation matrix.
   * 
   * It shows how to read in a JPG image, how to compute a noise pattern.
   * Modify the code to compute the correlation matrix.
   */
  
  /* Create a domain for an image and allocate the image itself */
  for i in 1..imageFileNames.size do {
    writeln("Outer file  " , i);

        var image : [imageDomain] RGB; 
        image = images[i];
        /* Read in the first image. */
        var data : prnu_data;
        var prnuc : [imageDomain] real;
        var prnucomp : [imageDomain] complex;

        prnuInit(h, w, data);
        prnuExecute(prnuc, image, data);

        // for ii in 0..#w do {
        //   for jj in 0..#h do {
        //     prnucomp[jj, ii].re = prnuc[jj, ii];
        //     prnucomp[jj, ii].im = 0.0;
        //   }
        // }
        prnucomp = prnuc;
        
      coforall j in i + 1..imageFileNames.size do {
        var image2 : [imageDomain] RGB;

        /* Read in the first image. */
        image2 = images[j];
        var data2 : prnu_data;
        var prnu2 : [imageDomain] real;
        var prnu2rot : [imageDomain] complex;
        var product : [imageDomain] complex;

        prnuInit(h, w, data2);
        prnuExecute(prnu2, image2, data2);

  
        // Rotating the second prnu image and representing it as matrix of complex numbers
        coforall ii in 0..h - 1 do {
          for jj in 0..w - 1 do {
            var newy = h - ii - 1;
            var newx = w - jj - 1;
            prnu2rot[newy, newx].re = prnu2[ii, jj];
            prnu2rot[newy, newx].im = 0.0;
          }
        }

      // writeln("Before FTT element 100 100 == " , prnucomp[100,100]);
      // writeln("Before FTT rotated element 100 100 == " , prnu2rot[100,100]);

        
        var planForward = plan_dft(prnucomp, prnucomp, FFTW_FORWARD, FFTW_ESTIMATE);
        var planForward2 = plan_dft(prnu2rot, prnu2rot, FFTW_FORWARD, FFTW_ESTIMATE);
     
        execute(planForward);
        execute(planForward2);
        // writeln("After FTT element 100 100 == " , prnucomp[100,100]);
        // writeln("After FTT rotated element 100 100 == " , prnu2rot[100,100]);


        /* allocate a prnu_data record */
        // dotProduct(product, prnucomp, prnu2rot);

        forall (row, col) in imageDomain {
            var tmp:complex = 0 + 0i;
            tmp.re = (prnucomp[row, col].re * prnu2rot[row, col].re - prnucomp[row, col].im * prnu2rot[row, col].im);
            tmp.im = prnucomp[row, col].im *  prnu2rot[row, col].re + prnucomp[row, col].re * prnu2rot[row, col].im;
            product[row, col] = tmp;
        }
        // writeln("After Cross corellation element 100 100 == " , product[100,100]);


        var planBackward = plan_dft(product, product, FFTW_BACKWARD, FFTW_ESTIMATE);

        execute(planBackward);

        product = product / (w * h);
        // writeln("After IFFT and scaling element 100 100 == " , product[100,100]);


        var (maxVal, maxLoc) = maxloc reduce zip(abs(product), product.domain);
        
        // writeln("peak coord " , maxLoc);
        // writeln("peak val " , maxVal);


        var sum:real = 0;
        var totalNum = 0;
        for ii in 0..#w do {
          for jj in 0..#h do {
            if (abs(jj - maxLoc[1]) > 5 || abs(ii - maxLoc[2]) > 5) && !isnan(product(jj,ii).re){
              sum += product(jj,ii).re**2;
              totalNum += 1;
            }
          }
        }
        writeln("Energy " ,sum);

        sum /= totalNum;
        writeln("Energy weighted " ,sum);

        corrMatrix[i,j] = maxVal * maxVal / sum;

        corrMatrix[j,i] = maxVal * maxVal / sum;

        prnuDestroy(data2);
       
    }
     prnuDestroy(data);
  }

  overallTimer.stop();

  writeln("The first value of the corrMatrix is: ", corrMatrix[2,1]);
  writeln("Time: ", overallTimer.elapsed(), "s");
  var nrCorrelations = (n * (n - 1)) / 2;
  writeln("Throughput: ", nrCorrelations / overallTimer.elapsed(), " corrs/s");

  if (writeOutput) {
    writeln("Writing output files...");
    write2DRealArray(corrMatrix, "corrMatrix");
    // for now, also write the prnu noise pattern, can be removed
    // write2DRealArray(prnucomp, "prnu");
  }
}