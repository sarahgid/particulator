float[] firstRow(float[][] A) {
  return A[0];
}

float[] firstCol(float[][] A) {
  float[] x = new float[A.length];
  for (int i=0; i<x.length; i++) x[i] = A[i][0];
  return x;
}



float[] zeros(int I) {
  float[] a = new float[I];
  for (int i=0; i<I; i++) a[i] = 0;
  return a;
}
float[][] zeros(int J, int I) {
  float[][] a = new float[J][I];
  for (int j=0; j<J; j++) for (int i=0; i<I; i++) a[j][i] = 0;
  return a;
}
float[][][] zeros(int K, int J, int I) {
  float[][][] a = new float[K][J][I];
  for (int k=0; k<K; k++) for (int j=0; j<J; j++) for (int i=0; i<I; i++) a[k][j][i] = 0;
  return a;
}

float[] arrayFill(int I, float val) {
  float[] a = new float[I];
  for (int i=0; i<I; i++) a[i] = val;
  return a;
}
float[][] arrayFill(int J, int I, float val) {
  float[][] a = new float[J][I];
  for (int j=0; j<J; j++) for (int i=0; i<I; i++) a[j][i] = val;
  return a;
}
float[][][] arrayFill(int K, int J, int I, float val) {
  float[][][] a = new float[K][J][I];
  for (int k=0; k<K; k++) for (int j=0; j<J; j++) for (int i=0; i<I; i++) a[k][j][i] = val;
  return a;
}



float Inf = 1./0.;
float almostInf = 1e20;
boolean isfinite(float a) {
  return (a<almostInf) && (a>-almostInf);
}
boolean isnan(float a) {
  return !((a>0) || (a<=0));
}


void finitize(float[] a) {
  for (int i=0; i<a.length; i++) {
    if (!isfinite(a[i])) a[i] = 0; 
  }
}
void finitize(float[][] a) {
  for (int j=0; j<a.length; j++) {
    for (int i=0; i<a[0].length; i++) {
      if (!isfinite(a[j][i])) a[j][i] = 0; 
    }
  }
}
void finitize(float[][][] a) {
  for (int k=0; k<a.length; k++) {
    for (int j=0; j<a[0].length; j++) {
      for (int i=0; i<a[0][0].length; i++) {
        if (!isfinite(a[k][j][i])) a[k][j][i] = 0; 
      }
    }
  }
}



int findIndexBefore(float[] x, float xi) {
  // assumes x is monotonic and increasing
  // return -1 if xi is out of range
  if (xi < x[0]) return -1;
  if (xi > x[x.length-1]) return -1;
  int nbefore = 0;
  int nafter = x.length-1;
  while (nafter-nbefore > 1) {
    int nmid = (nbefore+nafter)/2;
    if (xi > x[nmid]) {
      nbefore = nmid;
    } else {
      nafter = nmid;
    }
  }
  return nbefore;
}
