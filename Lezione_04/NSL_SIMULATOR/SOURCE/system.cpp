/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    // TO BE FIXED IN EXERCISE 6
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize( double counter, Mode mode){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/"+_phase+"/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/"+_phase+"/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  string filename_input;
  string filename_outputa;
  //Switch in order to divide equilibration case and evolution case: divide input file and output_acceptance
  switch (mode)
    {
    case MODE1:
      filename_input = "../INPUT/"+_phase+"/input_equilibration.dat";
      filename_outputa = "../OUTPUT/"+_phase+"/EQUILIBRATION/acceptance.dat";
    break;
    case MODE2:
      filename_input = "../INPUT/"+_phase+"/input.dat";
      filename_outputa = "../OUTPUT/"+_phase+"/acceptance.dat";
    break;
    default:
      break;
    }
  ofstream couta(filename_outputa); // Set the heading line in file ../OUTPUT/acceptance.dat or ../OUTPUT/EQUILIBRATION_acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();
    std::ifstream inputFile(filename_input);
    if (!inputFile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename_input << std::endl;
        exit(-1);  // Esce dal programma con codice di errore
    }
  ifstream input(filename_input); // Start reading ../INPUT/input.dat
  ofstream coutf;
  string output_f;
  double t_counter= counter; 
        ostringstream oss;
        oss << fixed << setprecision(2) << t_counter;
        string formatted_counter = oss.str();
  //Switch in order to divide equilibration case and evolution case: divide output_file in a direct for equilibration
    switch (mode){
      case MODE1:
      output_f = "../OUTPUT/"+_phase +"/EQUILIBRATION/output_" + formatted_counter + ".dat";
      break;
      
      case MODE2:
      output_f = "../OUTPUT/"+_phase +"/output_" + formatted_counter + ".dat";
      default:
      break;
  }
  coutf.open(output_f);
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
    //switch: equilibration, increase or decrease temperature of a fix factor, for searching adapting temperature
    //        evolution: chosen a temperature or a set of temperature evolved them.
      switch (mode){
        case MODE1:
        input >> _temp;
        if (counter != 1.){
          if (_phase == "SOLID" ||_phase == "LIQUID" ){
            if (_T_average < _T_attesa + _T_precision){
              _temp += 0.01 * (counter - 1.);
              } 
          }
        if (_phase == "GAS"){
            if(_T_average  > _T_attesa - _T_precision){
        _temp = _temp - 0.01 * (counter - 1.);   
            }
        }
      }
        break;
      case MODE2:
        input >> _temp;
        break;
      
      default:
        break;
      }
     _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);  //resize è una funzione di armadillo che rigenera la grandezza del vettoere
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_velocities(){
  if(_restart and _sim_type==0){
    ifstream cinf;
    cinf.open("../INPUT/"+_phase +"/CONFIG/velocities.in");
    if(cinf.is_open()){
      double vx, vy, vz;
      for(int i=0; i<_npart; i++){
        cinf >> vx >> vy >> vz;
        _particle(i).setvelocity(0,vx);
        _particle(i).setvelocity(1,vy);
        _particle(i).setvelocity(2,vz);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file velocities.in"<< endl;
    cinf.close();
  } else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0); //this means to put in a sdr where emerged only the fluctuation
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef); //reducing velocity, scale factor like norma --> reduced velocity
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
  }
  if(_sim_type == 0){
  double xold, yold, zold;
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System :: initialize_properties(double counter, Mode mode){ // Initialize data members used for measurement of properties
  string property;
  int index_property = 0;
  _nprop = 0;
  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  
  ifstream input;
  switch (mode){
    case MODE1:{
    input.open("../INPUT/"+ _phase +"/properties_equilibration.dat");
    break;
    }
    case MODE2:{
    input.open("../INPUT/"+ _phase +"/properties.dat");
    break;
    }
    default :{
      cout << "Error in reading System:: initialise_properties in System.cpp ";
      break;
    }
  }
  
  if (input.is_open()){
  
  while ( !input.eof() ){
      input >> property;
          if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/"+ _phase +"/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/"+ _phase +"/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/"+ _phase +"/total_energy.dat");
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        double _t_temporanea = _temp; 
        ostringstream oss;
        oss << fixed << setprecision(2) << _t_temporanea;
        string formatted_temp = oss.str();
        ofstream coutte;
        
        switch (mode){
        case MODE1:{
          coutte.open("../OUTPUT/"+ _phase +"/EQUILIBRATION/temperature_" + formatted_temp + ".dat");
          break;
        }
        case MODE2:{
          coutte.open("../OUTPUT/"+ _phase +"/temperature_" + formatted_temp + ".dat");
        break;
        }
        default:{
          cout << "Error in open file temperature.dat ";
          break;
        }
        }
        
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/"+ _phase +"/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/"+ _phase +"/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/"+ _phase +"/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/"+ _phase +"/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/"+ _phase +"/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        double t_counter = counter; 
        ostringstream oss;
        oss << fixed << setprecision(2) << t_counter;
        string formatted_counter = oss.str();

          switch (mode){
        case MODE1:{
          coutf.open("../OUTPUT/"+ _phase +"/EQUILIBRATION/output_" + formatted_counter + ".dat",ios::app);
          break;
        }
        case MODE2:{
          coutf.open("../OUTPUT/"+ _phase +"/output_" + formatted_counter + ".dat",ios::app);
        break;
        }
        default:{
          cout << "Error in open ENDPROPERTIES ";
          break;
        }
        }
        coutf.open("../OUTPUT/"+ _phase +"/output_" + formatted_counter + ".dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    } 
     input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;
  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System :: finalize(double counter, Mode mode){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
    double t_counter = counter; 
    ostringstream oss;
    oss << fixed << setprecision(2) << t_counter;
    string formatted_counter = oss.str();
    string output_f;
      switch (mode){
      case MODE1:
      output_f = "../OUTPUT/"+ _phase +"/EQUILIBRATION/output_" + formatted_counter + ".dat";
      break;
      case MODE2:
      output_f = "../OUTPUT/"+ _phase +"/output_" + formatted_counter + ".dat";
      break;
      default:
      break;
  }
    coutf.open(output_f,ios::app);
  
  coutf << "Simulation completed!" << endl;
  coutf.close();
  set_best_diff();
  return;
}

// Write current configuration as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/"+ _phase +"/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setw(16) << _particle(i).getposition(0,true)/_side(0)          // x
              << setw(16) << _particle(i).getposition(1,true)/_side(1)          // y
              << setw(16) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    this->write_velocities();
  } else {
    coutf.open("../OUTPUT/"+ _phase +"/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/"+ _phase +"/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

void System :: write_velocities(){
  ofstream coutf;
  coutf.open("../OUTPUT/"+ _phase +"/CONFIG/velocities.out");
  if(coutf.is_open()){
    for(int i=0; i<_npart; i++){
      coutf << setw(16) << _particle(i).getvelocity(0)          // vx
            << setw(16) << _particle(i).getvelocity(1)          // vy
            << setw(16) << _particle(i).getvelocity(2) << endl; // vz
    }
  } else cerr << "PROBLEM: Unable to open velocities.dat" << endl;
  coutf.close();
  return;
}

// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/"+ _phase +"/CONFIG/config.xyz");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/"+ _phase +"/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk, Mode mode, double counter){ // Reset block accumulators to zero
  ofstream coutf;
  string output_f;
  double t_counter = counter; 
  ostringstream oss;
  oss << fixed << setprecision(2) << t_counter;
  string formatted_counter = oss.str();
  if(blk>0){
      switch (mode){
      case MODE1:
      output_f = "../OUTPUT/"+ _phase +"/EQUILIBRATION/output_" + formatted_counter + ".dat";
      break;
      
      case MODE2:
      output_f = "../OUTPUT/"+ _phase +"/output_" + formatted_counter + ".dat";
      default:
      break;
  }
    coutf.open(output_f,ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double magnetization=0.0;
  double virial=0.0;

  //ESERCIZIO 4 -- PRESSION
  double pression_temp = 0.0;
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(dr < _r_cut){
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          // PRESSURE ... TO BE FIXED IN EXERCISE 4 

          if (_measure_pressure) pression_temp += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        
        }
      }
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_pressure) _measurement(_index_pressure) = _rho*_measurement(_index_temp)+ (48.0/(3.0 * _volume))*pression_temp/(double)(_npart);

// TO BE FIXED IN EXERCISE 4
  // MAGNETIZATION /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6

  _block_av += _measurement; //Update block accumulators
  return;
}

void System :: averages(int blk, Mode mode){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);
  _global_av  += _average;
  _global_av2 += _average % _average; // % -> element-wise multiplication


  switch(mode) {
  ///EQUILIBRATION CASE //////////////////////////////////////////////////////////////////// 
    case MODE1 :{
  // TEMPERATURE /////////////////////////////////////////////////////////////////////
    if (_measure_temp){
      double _t_temporanea = _temp;
      bool k = true;
      ostringstream oss;
      oss << fixed << std::setprecision(2) << _t_temporanea;
      string formatted_temp = oss.str();
      coutf.open("../OUTPUT/"+ _phase +"/EQUILIBRATION/temperature_" + formatted_temp + ".dat",ios::app);

      average  = _average(_index_temp);
      sum_average = _global_av(_index_temp);
      sum_ave2 = _global_av2(_index_temp);
      coutf << setw(12) << blk
            << setw(12) << average
            << setw(12) << sum_average/double(blk)
            << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
      coutf.close();
    }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
      double fraction;
      coutf.open("../OUTPUT/"+ _phase +"/EQUILIBRATION/acceptance.dat",ios::app);
      if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
      else fraction = 0.0; 
      coutf << setw(12) << blk << setw(12) << fraction << endl;
      coutf.close();
    break;
  }
    
///EVOLUTION CASE ////////////////////////////////////////////////////////////////////
    case MODE2:{
      // POTENTIAL ENERGY //////////////////////////////////////////////////////////
      if (_measure_penergy){
        coutf.open("../OUTPUT/"+ _phase +"/potential_energy.dat",ios::app);
        average  = _average(_index_penergy);
        sum_average = _global_av(_index_penergy);
        sum_ave2 = _global_av2(_index_penergy);
        coutf << setw(12) << blk 
              << setw(12) << average
              << setw(12) << sum_average/double(blk)
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
      }
      // KINETIC ENERGY ////////////////////////////////////////////////////////////
      if (_measure_kenergy){
        coutf.open("../OUTPUT/"+ _phase +"/kinetic_energy.dat",ios::app);
        average  = _average(_index_kenergy);
        sum_average = _global_av(_index_kenergy);
        sum_ave2 = _global_av2(_index_kenergy);
        coutf << setw(12) << blk
              << setw(12) << average
              << setw(12) << sum_average/double(blk)
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
      }
      // TOTAL ENERGY //////////////////////////////////////////////////////////////
      if (_measure_tenergy){
        coutf.open("../OUTPUT/"+ _phase +"/total_energy.dat",ios::app);
        average  = _average(_index_tenergy);
        sum_average = _global_av(_index_tenergy);
        sum_ave2 = _global_av2(_index_tenergy);
        coutf << setw(12) << blk
              << setw(12) << average
              << setw(12) << sum_average/double(blk)
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
      }
      // TEMPERATURE ///////////////////////////////////////////////////////////////
      if (_measure_temp){
        double _t_temporanea = _temp; 
        ostringstream oss;
        oss << fixed << std::setprecision(2) << _t_temporanea;
        string formatted_temp = oss.str();
        coutf.open("../OUTPUT/"+ _phase +"/temperature_" + formatted_temp +".dat",ios::app);
        average  = _average(_index_temp);
        sum_average = _global_av(_index_temp);
        sum_ave2 = _global_av2(_index_temp);
        coutf << setw(12) << blk
              << setw(12) << average
              << setw(12) << sum_average/double(blk)
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
      }
      // PRESSURE //////////////////////////////////////////////////////////////////
      // TO BE FIXED IN EXERCISE 4
      if (_measure_pressure){
          coutf.open("../OUTPUT/"+ _phase +"/pressure.dat",ios::app);
        average  = _average(_index_pressure);
        sum_average = _global_av(_index_pressure);
        sum_ave2 = _global_av2(_index_pressure);
        coutf << setw(12) << blk
              << setw(12) << average
              << setw(12) << sum_average/double(blk)
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
      }
      // GOFR //////////////////////////////////////////////////////////////////////
      // TO BE FIXED IN EXERCISE 7
      // MAGNETIZATION /////////////////////////////////////////////////////////////
      // TO BE FIXED IN EXERCISE 6
      // SPECIFIC HEAT /////////////////////////////////////////////////////////////
      // TO BE FIXED IN EXERCISE 6
      // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
      // TO BE FIXED IN EXERCISE 6
      // ACCEPTANCE ////////////////////////////////////////////////////////////////
      double fraction;
      coutf.open("../OUTPUT/"+ _phase +"/acceptance.dat",ios::app);
      if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
      else fraction = 0.0; 
      coutf << setw(12) << blk << setw(12) << fraction << endl;
      coutf.close();

    break;
    } //END MODE 2
    
    default :{
    cout << "Error in reading System:: Averages in System.cpp ";
      break;
    }  
  }

  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

//function to compare attended temperature with effective temperature of the system after the equilibration
bool System :: equilibration_fine(double att_Temp,double  counter ){
      double ave_temp,err_temp;
      double diff;
      ofstream cout_t;
      ave_temp = _global_av(_index_temp)/double(get_nbl());
      err_temp = error(_global_av(_index_temp), _global_av2(_index_temp), (get_nbl()));
      bool k = true;
      _T_average = ave_temp;
      cout << _temp << " : " << ave_temp << " "<< err_temp << endl;
      diff = fabs(ave_temp + err_temp - att_Temp);
      if ( diff  <= _T_precision ){ 
        double T_star = 0;
        k = false;
        cout_t << _temp << endl;
        cout_t.close();
        if (counter == 1. ){
          _best_diff = diff ;
          modify_temp_input("../INPUT/"+ _phase +"/input.dat", _temp );
          }
        if ( counter != 1. && diff < _best_diff ){
          _best_diff = diff;
          modify_temp_input("../INPUT/"+ _phase +"/input.dat", _temp );
        }
        return k; //if return true, I can stop the measure
      }
//finchè la temperatura è compresa in questo intervallo stampa
      cout_t.close();
      return k; 
  }

bool System :: equilibration_rough(double att_Temp, double counter ){
      double ave_temp,err_temp;
      ofstream cout_t;
      ave_temp = _global_av(_index_temp)/double(get_nbl());
      err_temp = error(_global_av(_index_temp), _global_av2(_index_temp), (get_nbl()));
      bool k = true;
      _T_average = ave_temp;
      cout << _temp << " : " << ave_temp << " "<< err_temp<< endl;
      double diff = fabs(ave_temp + err_temp - att_Temp);
      cout << "diff " << diff <<endl;

      if ( diff > _T_precision ){ 
        k = false;
        cout_t << _temp << endl;
        cout_t.close();
        return k; //if return true, I can stop the measure
      }
//finchè la temperatura è compresa in questo intervallo stampa
      cout_t.close();
      return k; 
  }

void System :: modify_temp_input(const string& filename, double T_STAR) {
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return;
    }
    ostringstream tempStream;
    string line;

    while (getline(inputFile, line)) {
        // Searching line named "TEMP"
        if (line.find("TEMP") != string::npos) {
            // Costruisci la nuova riga con T_STAR
            ostringstream newLine;
            newLine << "TEMP                   " << T_STAR;
            tempStream << newLine.str() << endl;
        } else {
            // copy line without modifying
            tempStream << line << endl;
        }
    }

    inputFile.close();

    // Scrivi il contenuto modificato in un nuovo file (o sovrascrivi quello esistente)
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the file for writing " << filename << std::endl;
        return;
    }

    outputFile << tempStream.str();
    outputFile.close();
}

void System :: set_best_diff(){
  _best_diff = 0. ;
}

string System :: get_phase(){
  return _phase;
}

void System :: set_phase_of_system(string phase){
  _phase = phase;
}

double System :: _set_T_attended(string phase){
  if (_phase == "GAS"){
    _T_attesa = 1.2;
    return _T_attesa;
  } else if( _phase == "LIQUID"){
    _T_attesa = 1.1;
    return _T_attesa;
  } else if (_phase == "SOLID"){
    _T_attesa = 0.8;
    return _T_attesa;
  }
  cout << "Error in setting _T_attended"<< endl;
  exit (-1);
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
