// particulator 2.8.1
// neil banas, feb 2011

/* new
  getm2d loads one time slice at a time
*/

/* to do
  
  test higher-order stepping schemes (start with AB3)
  do a proper convergence test in the coastal domain; Puget Sound is so chaotically dispersive that it's a hard place to look at this
  divide seedParticles into cs and z versions
  less hokey way to determine dz in dKs/dz
  allow particles to be depth-trapped, not just surface-trapped
*/


boolean debug = true;
boolean useConfigFile = true;


void setup() {
  float tic = millis();
  if (useConfigFile) {
    Configuration config = new Configuration("particulator_config.txt");
    String exptname = config.getString("exptname");
    if (debug) config.echo();
    
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
