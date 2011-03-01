import ucar.*;
import java.io.IOException;
import org.slf4j.LoggerFactory;


// opening and closing -----------------------------------------------


NetcdfFile nc_open(String filename) {
  NetcdfFile nc = null;
  try {
    nc = NetcdfFile.open(filename);
  } catch (IOException ioe) {
    println("error opening " + filename + ": " + ioe.toString());
  }
  return nc;
}


void nc_close(NetcdfFile nc) {
  if (nc==null) return;
  try {
    nc.close();
  } catch (IOException ioe) {
    println(ioe.toString());
  }
}


Variable nc_getVar(NetcdfFile nc, String varname) {
  if (nc==null) return null;
  Variable v = nc.findVariable(varname);
  if (v==null) {
//    println(varname + " not found");
    return null;
  } else {
    return v;
  }
}



// reading entire variables -----------------------------------------



float nc_readOne(NetcdfFile nc, String varname) {
  float result = 0./0;
  Variable v = nc_getVar(nc, varname);
  if (v != null) {
    try {
      result = v.readScalarFloat();      
    } catch (Exception e) {
      println(e.toString());
    } 
  }
  return result;
}



float[] nc_read1D(NetcdfFile nc, String varname) {
  float[] result = null;
  Variable v = nc_getVar(nc, varname);
  if (v != null) {
    try {
      Array a = v.read().reduce(); 
      int[] s = a.getShape();
      result = new float[s[0]];
      Index ind = a.getIndex();
      for (int i=0; i<s[0]; i++) {
        result[i] = a.getFloat(ind.set(i));
      }
    } catch (Exception e) {
      println(e.toString());
    } 
  }
  return result;
}

  
float[][] nc_read2D(NetcdfFile nc, String varname) {
  float[][] result = null;
  Variable v = nc_getVar(nc, varname);
  if (v != null) {
    try {
      Array a = v.read().reduce(); 
      int[] s = a.getShape();
      result = new float[s[0]][s[1]];
      Index ind = a.getIndex();
      for (int j=0; j<s[0]; j++) {
        for (int i=0; i<s[1]; i++) {
          result[j][i] = a.getFloat(ind.set(j,i));
        }
      }
    } catch (Exception e) {
      println(e.toString());
    } 
  }
  return result;
}


float[][] nc_read2Dfrom3D(NetcdfFile nc, String varname, int k0) {
  float[][] result = null;
  Variable v = nc_getVar(nc, varname);
  if (v != null) {
    try {
      int[] s = v.getShape();
      int[] origin = {k0, 0, 0};
      int[] shape = {1, s[1], s[2]};
      Array a = v.read(origin, shape).reduce(); 
      s = a.getShape();
      result = new float[s[0]][s[1]];
      Index ind = a.getIndex();
      for (int j=0; j<s[0]; j++) {
        for (int i=0; i<s[1]; i++) {
          result[j][i] = a.getFloat(ind.set(j,i));
        }
      }
    } catch (Exception e) {
      println(e.toString());
    } 
  }
  return result;
}



float[][][] nc_read3D(NetcdfFile nc, String varname) {
  float[][][] result = null;
  Variable v = nc_getVar(nc, varname);
  if (v != null) {
    try {
      Array a = v.read().reduce(); 
      int[] s = a.getShape();
      result = new float[s[0]][s[1]][s[2]];
      Index ind = a.getIndex();
      for (int k=0; k<s[0]; k++) {
        for (int j=0; j<s[1]; j++) {
          for (int i=0; i<s[2]; i++) {
            result[k][j][i] = a.getFloat(ind.set(k,j,i));
          }
        }
      }
    } catch (Exception e) {
      println(e.toString());
    } 
  }
  return result;
}
