R__LOAD_LIBRARY(libeicsmear);
//R__LOAD_LIBRARY(/sphenix/user/baschmoo/myfork/dis-reconstruction/detectors/smearHandBook_cxx.so); //For sphenix account
R__LOAD_LIBRARY(/eic/data/baraks/dis-reconstruction/detectors/smearHandBook_cxx.so); //For eic account
//R__LOAD_LIBRARY(/work/halla/gmp12/baraks/dis-reconstruction/detectors/smearHandBook_cxx.so); //For JLAB account

void make_smeared_handbook(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string dirstr = "outfiles";
  std:string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_handbook_smeared.root";

  SmearTree(BuildHandBookDetector(),inputstr,outputstr);

}
