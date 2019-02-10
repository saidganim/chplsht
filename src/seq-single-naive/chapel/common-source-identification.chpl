/* Use assertions. */
use Assert;

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

proc dotProduct(ref C: [?DC] complex, ref A: [?DA] complex, ref B: [?DB] complex){

  forall (row, col) in DC {
    // Zero out the value, in case C is reused.
    C[row, col].re = 0;
    C[row, col].im = 0;
    for i in DA.dim(2) do{
      var tmp:complex = 0 + 0i;
      tmp.re = (A[row, i].re * B[i, col].re - A[row, i].im * B[i, col].im);
      tmp.im = A[row, i].im *  B[i, col].re + A[row, i].re * B[i, col].im;
      C[row, col] += tmp;
    }
  }
}

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
  const corrDomain : domain(2) = {1..n, 1..n};
  var corrMatrix : [corrDomain] real;

  var overallTimer : Timer;

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
  for i in 0..imageFileNames.size do {
     const imageDomain: domain(2) = {0..#h,0..#w};
        var image : [imageDomain] RGB;

        /* Read in the first image. */
        readJPG(image, imageFileNames.pop_front());
        var data : prnu_data;
        var prnuc : [imageDomain] real;
        var prnucomp : [imageDomain] complex;

        prnuInit(h, w, data);
        prnuExecute(prnuc, image, data);

        for ii in 0..#w do {
          for jj in 0..#h do {
            prnucomp[ii, jj].re = prnuc[ii, jj];
            prnucomp[ii, jj].im = 0.0;
          }
        }
        
      for j in i + 1..imageFileNames.size do {
        var image2 : [imageDomain] RGB;

        /* Read in the first image. */
        readJPG(image2, imageFileNames.pop_front());
        var data2 : prnu_data;
        var prnu2 : [imageDomain] real;
        var prnu2rot : [imageDomain] complex;
        var product : [imageDomain] complex;

        prnuInit(h, w, data2);
        prnuExecute(prnu2, image2, data2);

        // Rotating the second prnu image and representing it as matrix of complex numbers
        for ii in 0..#w do {
          for jj in 0..#h do {
            var newx = w - ii;
            var newy = h - jj;
            prnu2rot[newx, newy].re = prnu2[ii, jj];
            prnu2rot[newx, newy].im = 0.0;
          }
        }
        
        var planForward = plan_dft(prnucomp, prnucomp, FFTW_FORWARD, FFTW_ESTIMATE);
        var planForward2 = plan_dft(prnu2rot, prnu2rot, FFTW_FORWARD, FFTW_ESTIMATE);
        execute(planForward);
        execute(planForward2);
        /* allocate a prnu_data record */
        dotProduct(product, prnucomp, prnu2rot);
        var planBackward = plan_dft(product, product, FFTW_BACKWARD, FFTW_ESTIMATE);
        execute(planBackward);
        product = product / (w * h);
        var (maxVal, maxLoc) = maxloc reduce zip(abs(product), product.domain);
        var sum = 0;
        var totalNum = 0;
        for ii in 0..#w do {
          for jj in 0..#h do {
            if abs(ii - mod(maxLoc,w)) > 11 || abs(jj - (maxLoc / w):int) > 11{
              sum += product(ii,jj);
              totalNum += 1;
            }
          }
        }
        sum /= totalNum;
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