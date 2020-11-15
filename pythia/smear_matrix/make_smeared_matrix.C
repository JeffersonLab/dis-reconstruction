R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(/eic/data/baraks/dis-reconstruction/detectors/SmearMatrixDetector_0_1_cxx.so); //For eic account

void make_smeared_matrix(std::string dirstr,std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
 
  std::string inputstr = dirstr + "/" + filstr + ".root";
  std::string outputstr = dirstr + "/" + filstr + "_matrixsmear.root";

  SmearTree(BuildMatrixDetector_0_1(),inputstr,outputstr);

}
