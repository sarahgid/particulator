// particulator 3.0a
// neil banas, jun 2012

/* new
  multithreaded!
  note that since the old autoSave method is apparently incompatible with multithreading, when multithreading is on,
  autoSave simply saves once per _frame of input data_, not particle timestep * saveInterval. The options for this
  need to be rethought.
*/

/* to do
  test higher-order stepping schemes (start with AB3)
  do a proper convergence test in the coastal domain; Puget Sound is so chaotically dispersive that it's a hard place to look at this
  divide seedParticles into cs and z versions
  less hokey way to determine dz in dKs/dz
*/


boolean debug = true;
boolean useConfigFile = true;


void setup() {
  float tic = millis();
  if (useConfigFile) {
    Configuration config = new Configuration("particulator_config.txt");
    String exptname = config.getString("exptname");
    if (debug) config.echo();
    
    if (exptname.equals("surfaceBox")) surfaceBox(config);
    if (exptname.equals("jdf2005")) jdf2005(config);
    if (exptname.equals("returnmaps2005")) returnmaps2005(config);
    if (exptname.equals("riverYear")) riverYear(config);
    if (exptname.equals("getmMap2D_set")) getmMap2D_set(config);
    
  } else {
    pugetSoundTest();
  }
  
  if (debug) println("elapsed time " + (millis()-tic)/1000. + " sec");
  exit();
}
