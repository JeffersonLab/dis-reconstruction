//This script assumes perfect PID and perfect angular resolution


int energy_resolution_optimization(double E, double P, double eta, int PID, bool IsESmeared, bool IsPSmeared){

  //Setting the resolutions large makes it so that if only E or P is smeared, it will automatically say that the resolution of the smeared one is better.
  double sigma_P(100.); 
  double sigma_E(100.);

  //Check Tracking Resolution (0)
  if(IsPSmeared){
    if( -3.5 <= eta && eta < -2.5 ){
      sigma_P = sqrt( pow ( 0.001*P*P, 2) + pow ( 0.005*P, 2) );
    } else if( -2.5 <= eta && eta < -1.0 ){
      sigma_P = sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) );
    } else if( -1.0 <= eta && eta < 1.0 ){
      sigma_P = sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) );
    } else if( 1.0 <= eta && eta < 2.5 ){
      sigma_P = sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) );
    } else if( 2.5 <= eta && eta < 3.5 ){
      sigma_P = sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) );
    }
  }

  //Check ECal Resolution (1)
  if( (PID == 11 || PID == 22 || PID == -11) && IsESmeared){
    if( -3.5 <= eta && eta < -2.0 ){
      sigma_E = sqrt( pow ( 0.0*E,2 ) + pow( 0.02,2)*E); 
    } else if( -2.0 <= eta && eta < -1.0){
      sigma_E = sqrt( pow ( 0.0*E,2 ) + pow( 0.07,2)*E);
    } else if( -1.0 <= eta && eta < 3.5 ){
      sigma_E = sqrt( pow ( 0.0*E,2 ) + pow( 0.12,2)*E);
    }
  }

  //Check HCal Resolution (1)
  if( (PID != 11 || PID != 22 || PID != -11) && IsESmeared){
    if( -3.5 <= eta && eta < -1.0 ){
      sigma_E = sqrt(pow( 0.0*E, 2) + pow ( 0.5,2) *E);
    } else if( -1.0 <= eta && eta < 1.0){
      sigma_E = sqrt( pow( 0.07*E, 2) + pow( 0.85, 2)*E);
    } else if( 1.0 <= eta && eta < 3.5 ){
      sigma_E = sqrt(pow( 0.0*E, 2) + pow ( 0.5,2) *E);
    }
  }

  if(sigma_P > sigma_E){ //If cal resolution better, return 0 
    return 0; 
  } else if(sigma_P < sigma_E){ //If track resolution better, return 1
    return 1;
  } else{
    return std::nan("");
  }

}
