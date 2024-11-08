//187-198
#include "system.h"

// Implementation of the uniform Metropolis algorithm
////////////////////////////////UNIFORME//////////////////////////////////////////

System::System() {
   
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   
   ifstream input("seed.in");
   string property;
   
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            RND.SetRandom(seed,p1,p2);
         }
      }
      input.close();
      
   }
   else cerr << "PROBLEM: Unable to open seed.in" << endl;
    }; // Constructor for uniform Metropolis

    // Get and set methods for the delta value in uniform distribution
    double System::get_delta() {
    return unif_delta; // Return the current sigma for uniform distribution
    };

    void System::set_delta(double x) {
    unif_delta = x; // Set a new sigma value for uniform distribution
    };

    // Get and set methods for old coordinates in uniform distribution
    double System::get_pos_old() {
    return pos_old; // Return the previously stored coordinates
    };

    void System::set_pos_new(double x) {
    pos_old = x;
     // Update the old coordinates with new ones
    };

    // Get and set methods for temperature in uniform distribution
    double System::get_temp() {
    return _temp; // Return the previously stored coordinates
    };

    void System::set_temp(double x) {
    _temp = x;
     // Update the old coordinates with new ones
    };

        // Get and set methods for temperature in uniform distribution
    double System::get_energy_old() {
    return _energy_old; // Return the previously stored coordinates
    };

    void System::set_energy_new(double x) {
    _energy_old = x;
     // Update the old coordinates with new ones
    };

    



    // Implementation of the Metropolis algorithm for uniform distribution
    bool System::metropolis_method(Waves_function* WF) {
    
    if (WF == nullptr) {
    cerr << "Error: null WF pointer passed to initialize_acceptance!" << endl;
    return -1;
    }
    
  
    // Propose new coordinates using uniform random distribution
    double coordinate_proposed = RND.Rannyu(-1,1)*unif_delta;

    //ofstream out_p;
    //out_p.open("./OUTPUT/position.dat",ios::app);
    

    
    //coordinate_proposed.print("coordinate del punto proposto: ");
    
    // Calculate wave function values at current and proposed coordinates
    double x = WF->measure(get_pos_old());
    //cout << "x: " << x << endl;
    double y = WF->measure(coordinate_proposed);
    //cout << "y: " << y << endl;
    double z = RND.Rannyu(0,1); // Generate a uniform random number in [0,1]
    //out_p << x <<setw(12) <<y << endl;

    // Metropolis acceptance criteria
    if (pow(y / x, 2) >= 1.) {
        set_pos_new(coordinate_proposed); // Accept new coordinates
        //cout << "coordinate accetate" << endl;
        //out_p.close();
        return true;
    }
    else if (pow(y / x, 2) < 1. && z > pow(y / x, 2)) {
        set_pos_new(get_pos_old()); // Keep old coordinates
        //cout << "coordinate non accetate" << endl;
        //out_p.close();
        return false;
    }
    if (pow(y / x, 2) < 1. && z <= pow(y / x, 2)){
        set_pos_new(coordinate_proposed); // Accept new coordinates
        //out_p.close();
        return true;
        //cout << "coordinate accetate" << endl;
        
    }

    return false;
    
    };

    void System::initialize_acceptance(Waves_function* WF, double coordinates) {
    double counter = 0.;
    ofstream out_acc;
    
    out_acc.open("./OUTPUT/accettanza_unif_.dat");
    bool repeat = false;
    double coordinates_true = coordinates;

    int step = 100000;
    double s = get_delta(); // Ottieni il sigma iniziale
    double target_acceptance = 0.5; // Target di accettazione del 50%
    double tolerance = 0.01; // Tolleranza per l'accettazione
    s = WF->get_mu()+3*WF->get_sigma();
    //cout << "prima sigma" << s << endl;
    out_acc << "Accettanza for uniform method is start:"<< endl << " delta "<< " " << " accettanza " << endl;
    
    do {
        
        counter = 0; //reinitialize the cycle
        //double x = get_pos_old();
          //  cout << "pos_iniziale" << x << endl;
    
        for (int i = 0; i < step; i++) {
            double x = get_pos_old();
            //x.print("x: ");
            metropolis_method(WF);
            double y = get_pos_old();
            //y.print("y: ");
            // Confronto tra vettori, if they are the same means that the vector proposed is not accepted
            if (x != y) {
                counter += 1;
            }
        }

        double acceptance_rate = counter / (double) step;
        //cout << "delta: "<< get_delta() << " " <<" accettanza del ciclo: " << acceptance_rate*100. << endl;
        out_acc << get_delta() << " " << acceptance_rate*100 << endl;
         if (acceptance_rate > target_acceptance + tolerance) {
        s *= 1.0 + 0.1 * (acceptance_rate - target_acceptance); // Aumenta delta proporzionalmente all'eccesso di accettazione
        }else if (acceptance_rate < target_acceptance - tolerance) {
        s *= 1.0 - 0.1 * (target_acceptance - acceptance_rate); // Diminuisci delta proporzionalmente al deficit di accettazione
        } 
        if (acceptance_rate > target_acceptance - tolerance && acceptance_rate < target_acceptance + tolerance){
            repeat = true; // Fermati se il tasso è dentro il range accettabile
        }

        set_delta(s); // Aggiorna sigma
        set_pos_new(coordinates_true);
        //cout << "sigma aggiornata: " << get_delta() << endl; 

        } while (repeat == false);

    
    cout << endl << "delta accettata per uniforme: " << get_delta() << endl; //finished the acceptance fase I reset the initial coordinates
    out_acc << endl << "ACCETANZA FINISHED!"<< endl;
    out_acc.close();
};


void System:: print_coordinate_accepted(double pos){
    
    ostringstream oss_x;
    oss_x << setprecision(1) << fixed << pos;
    ofstream out_cor;
    
    // Conversione da coordinate sferiche a cartesiane
    
    out_cor.open("./OUTPUT/Coordinate_accettate_uniforme_"+oss_x.str()+"_.dat", ios:: app);
    out_cor << get_pos_old() << endl;
    
    out_cor.close();
    
};

//SIMULATED ANNEALING

bool System:: Boltzmann(double temp, double Energy_old , double Energy_new){
    double _beta = 1./_temp;
    double z = RND.Rannyu(0,1);
    double boltz = exp(-1.*(Energy_new-Energy_old)*_beta);

    if(boltz >= 1.){return true;
    }else{
        if(z < boltz){return true;
        }else{return false;}
    }
    return false;
};

void System:: simulated_annealing(double temp, Waves_function* WF , Waves_function* Energy){

    double temp_energy = 0.;
    double accp = 0;
    bool acc_annealing = false;
    vector <double> mean_err(2, 0.0);
    string mu  = "mu";
    string sigma = "sigma";
    int max_cicle = 1000000;
    int counter = 0;
    
    initialize_acceptance(WF,pos_old);
    do{
        temp_energy = 0.;
        vector <double> average_energy;
        vector <double> average_energy_2;
    
    for (int i = 0; i < 100 ; i++ ){
        for (int j = 0; j< 20000; j++){

            bool measure_acc = metropolis_method(WF);
            if(measure_acc == true){
                temp_energy += Energy->measure(pos_old);
                    accp += 1;
                }
        }
    
    temp_energy /= accp;
    //cout << "temp_energy " << temp_energy << endl;
    average_energy.push_back(temp_energy);
    average_energy_2.push_back(pow(temp_energy,2));
    }
    cout << "lunghezza average " << average_energy.size() << endl;

    //mean_err=Bloch_mean_method(average_energy,average_energy_2,100, mean_err);

    
    
    acc_annealing = Boltzmann(temp, get_energy_old(), mean_err[0]);
    if (acc_annealing == true){
        set_energy_new(mean_err[0]);
        cout << "passo accettato: "<< endl;
        WF->set_mu(WF->get_mu()+RND.Rannyu(-1,1));
        Energy->set_mu(WF->get_mu()+RND.Rannyu(-1,1));
        WF->set_sigma(WF->get_sigma()+RND.Rannyu(-1,1));
        Energy->set_sigma(WF->get_sigma()+RND.Rannyu(-1,1));
        }
        counter+=1;
        if (counter == max_cicle){
            cerr << "Too much cicle to find a new energy state: " << endl;
            exit (-1);
        }
    } while(acc_annealing == false);


    print_on_file(mu, WF->get_mu());
    print_on_file(sigma, WF->get_sigma());
    print_on_file("Energy_ave", mean_err[0], mean_err[1]);


}

void System:: cooling(double temp){
    double temperature = get_temp();
    set_temp(temperature*0.99);
}

void System::print_on_file(const string& filename, double content) {
    ofstream output;
    
    // Verifica se il file già esiste e se è vuoto
    ifstream check_file("./OUTPUT/" + filename + ".dat");
    bool is_empty = check_file.peek() == ifstream::traits_type::eof();
    check_file.close();
    
    // Apre il file in modalità append
    output.open("./OUTPUT/" + filename + ".dat", ios::app);
    
    if (output.is_open()) {
        // Se il file è vuoto, scrivi il nome del file come prima riga
        if (is_empty) {
            output << "TEMPERATURE: " << filename << endl;
        }
        // Scrivi il contenuto successivo
        output << _temp << " " << content << endl;
        output.close();
    } else {
        cerr << "Errore nell'aprire il file: " << filename << ".dat" << endl;
    }
}

void System::print_on_file(const string& filename, double mean, double err) {
    ofstream output;
    
    // Verifica se il file già esiste e se è vuoto
    ifstream check_file("./OUTPUT/" + filename + ".dat");
    bool is_empty = check_file.peek() == ifstream::traits_type::eof();
    check_file.close();

    // Apre il file in modalità append
    output.open("./OUTPUT/" + filename + ".dat", ios::app);
    
    if (output.is_open()) {
        // Se il file è vuoto, scrivi il nome del file come prima riga
        if (is_empty) {
            output << "TEMPERATURE: " << "ENERGY_AVE: " << "ERR:AVE: " << endl;
        }
        // Scrivi il contenuto successivo
        output << _temp << " " << mean << " " << err << endl;
        cout << "temperatura: " << _temp << " media: " << mean << " errore: " << err << endl;
        output.close();
    } else {
        cerr << "Errore nell'aprire il file: " << filename << ".dat" << endl;
    }
};

void System:: initialize_print(){
    ofstream out_m, out_s, out_e, out_pos;

    unique_ptr<Waves_function> WF_test = unique_ptr<Wave_Variational_test>(new Wave_Variational_test());

    ostringstream oss_mu,oss_s;
    oss_mu << fixed << std::setprecision(1) << WF_test -> get_mu();
    oss_s << fixed << std::setprecision(1) << WF_test -> get_sigma();


    //out_m.open("./OUTPUT/mu.dat");
    //out_m << "TEMPERATURE: " << "MU: " << "ERR: " << endl;
    //out_s.open("./OUTPUT/sigma.dat");
    //out_s << "TEMPERATURE: " << "SIGMA: " << "ERR: " << endl;
    //out_e.open("./OUTPUT/Energy_ave.dat");
    //out_e << "TEMPERATURE: " << "ENERGY_AVE: " << "ERR: " << endl;
    out_pos.open("./OUTPUT/WaveFunction_"+oss_mu.str()+"_"+oss_s.str()+"_.dat");
    out_pos << "POSITION: " << endl;

    //out_m.close();
    //out_s.close();
    //out_e.close();
    out_pos.close();

};





