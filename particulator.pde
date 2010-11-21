// particulator 2.6
// neil banas, nov 2010

/* new
  GETMRun
  removing backwards-in-time tracking (that may or may not be working anyway)
*/

/* to do
  make it easy to load tracer data from the run, and to interpolate tracers along trajectories (or from the current x,y to a particular z or cs)
  test higher-order stepping schemes (start with AB3)
  do a proper convergence test in the coastal domain; Puget Sound is so chaotically dispersive that it's a hard place to look at this
  divide seedParticles into cs and z versions
  less hokey way to determine dz in dKs/dz
  allow particles to be depth-trapped, not just surface-trapped (figure out what to do in shallow water: cs = -0.9?)
*/


boolean debug = true;
boolean useConfigFile = false;


void setup() {
  float tic = millis();
  
  if (useConfigFile) {
    String[] config = loadStrings("particulator_config.txt");
    String exptname = trimComment(config[0]);
    String runDir = trimComment(config[1]);
    int fileStart = int(trimComment(config[2]));
    int fileEnd = int(trimComment(config[3]));
    String outputDir = trimComment(config[4]);
  
    if (debug) println("running expt " + exptname + " on " + runDir + " " + fileStart + ".." + fileEnd + ", output to " + outputDir);
    
    if (exptname.equals("jdf2005")) jdf2005(runDir, fileStart, fileEnd, outputDir);
    if (exptname.equals("returnmaps2005")) returnmaps2005(runDir, fileStart, fileEnd, outputDir);
    if (exptname.equals("riverYear")) riverYear(runDir, fileStart, fileEnd, outputDir);
    
  } else {
    getmMap2D_set("/Users/neil/Desktop/", 1, 1);
  }
  
  if (debug) println("elapsed time " + (millis()-tic)/1000. + " sec");
  exit();
}



String trimComment(String S) {
  String[] ss = split(S,'!');
  if (ss==null || ss.length==0) {
    return null;
  } else {
    return trim(ss[0]);
  }
}
