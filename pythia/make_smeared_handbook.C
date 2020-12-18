R__LOAD_LIBRARY(libeicsmear);
//need to update for JLab and SPhenix accounts
R__LOAD_LIBRARY(/eic/data/baraks/dis-reconstruction/detectors/SmearHandBook_1_2_cxx.so); //For eic account

void make_smeared_handbook(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string dirstr = "outfiles";
  std:string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_handbook_smeared.root";

  SmearTree(BuildHandBook_1_2(),inputstr,outputstr);

}
