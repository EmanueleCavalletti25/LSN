/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

//EQUILIBRATION for ISING: find the right configuration for the spin; for the evolution
//Evolution mode is the second part of programm where it is observed the evolution of the system before and after equilibration
Mode getModeFromArg(const string& arg) {
    if (arg == "equilibration") return MODE1;
    if (arg == "evolution") return MODE2;
    return UNKNOWN;
}

//PROGRAMM BEGIN
int main (int argc, char *argv[]){

   if (argc < 3){
      cout << "Usage, select: " << argv [0] << " <equilibration> or <evolution> ; <PHASE> or <ISING>" << endl;
      return -1;
   }
  //Setting program usage in based of flag: equilibration and evolution 
  Mode mode = getModeFromArg(argv[1]);
  //Setting function of controll
  if (mode == UNKNOWN){
    cout << "Usage, select: " << argv [0] << " <equilibration> or <evolution> " << endl;
      return -1;
  }
  
  System SYS;
  
  //Set problem with phase or Ising problem
  string phase = argv[2];
  while (phase != "GAS" && phase != "LIQUID" && phase != "SOLID" && phase != "ISING"){
    
    cerr << "Error, select the right phase GAS LIQUID or SOLID; or select the model ISING" << endl;
    return -1;
  }
  
  SYS.set_phase_of_system(phase);
  double counter = 1.;
  bool k = false;

  //do..while cycle: necessary expecially for the equilibration.
  do{
    
    int nconf = 1;
    double T_attesa = SYS._set_T_attended(phase);
    SYS.initialize( counter, mode);
    SYS.initialize_properties(counter,mode);
    SYS.block_reset(0,mode,counter);
    cout << "temperatura ciclo: " << SYS.get_temp() << endl;
    

    switch (mode){
    case MODE1 :
    break;

    case MODE2 : 
    int threshold;
    if (phase == "GAS"){threshold = 400;}
    if (phase == "SOLID" || phase == "LIQUID"){threshold = 200;}
    if (phase == "ISING"){threshold = 100;}
     
      for (int i=0; i<threshold; i++){
         for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            
          }
        cout << i << " block ok!" << endl;
        }
   
    default :
    break;
  }
    
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
          //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
      cout << i << " block ok!" << endl;
   
      SYS.averages(i+1 , mode);
      
      SYS.block_reset(i+1, mode , counter);
    }
    
    switch(mode) {
    //The first case I have to see that temperature equilibrate to the desired temperature of the system
     case MODE1 :
    
    //Chose between rough or fine equilibration
    // k = SYS.equilibration_rough(T_attesa, counter);
    if (phase != "ISING"){
    k = SYS.equilibration_fine(T_attesa, counter);
    } else {k = true;}
     
     if (k == false){
      counter += 1.;  
      }
     
     break;
    
    case MODE2:
      k = true;
    if (phase == "ISING"){
      while (SYS.get_temp() <= 1.9){
        SYS.modify_temp_input("../INPUT/"+ SYS.get_phase() +"/input.dat", SYS.get_temp()+0.1 );
        k=false;
        if (SYS.get_temp()>2.0){k=true;}
        SYS.write_configuration();
        break;
        }
      }
      
      break;
    
    default :
      break;

  }

  } while( k == false );

  SYS.finalize(counter,mode);
  counter = 1.;
  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
