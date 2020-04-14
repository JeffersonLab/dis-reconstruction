R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(/sphenix/user/baschmoo/myfork/dis-reconstruction/detectors/detector_perfect_C.so); //For sphenix account
//R__LOAD_LIBRARY(/eic/data/baraks/dis-reconstruction/detectors/detector_perfect_C.so); //For eic account
//R__LOAD_LIBRARY(/work/halla/gmp12/baraks/dis-reconstruction/detectors/detector_perfect_C.so); //For JLAB account

void make_smeared_perfect(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false; //Need to comment this for eic account 
 
  std::string dirstr = "outfiles";
  std::string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_perfect_smeared.root";

  SmearTree(BuildDetector(),inputstr,outputstr);

}
