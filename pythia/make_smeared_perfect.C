R__LOAD_LIBRARY(libeicsmear.so);
R__LOAD_LIBRARY(/sphenix/user/baschmoo/myfork/dis-reconstruction/detectors/detector_perfect_C.so);

void make_smeared_perfect(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false; 
 
  std::string dirstr = "outfiles";
  std::string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_perfect_smeared.root";

  SmearTree(BuildDetector(),inputstr,outputstr);

}
