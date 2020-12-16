R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(/eic/data/baraks/dis-reconstruction/detectors/detector_perfect_C.so); //For eic account

void make_smeared_perfect(std::string dirstr,std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_perfect_smeared.root";

  SmearTree(BuildDetector(),inputstr,outputstr);

}
