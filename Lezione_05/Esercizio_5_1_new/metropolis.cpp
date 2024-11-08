//187-198
#include "metropolis.h"

// Implementation of the uniform Metropolis algorithm
////////////////////////////////UNIFORME//////////////////////////////////////////

metropolis_unif::metropolis_unif() {
   
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

    // Get and set methods for the sigma value in uniform distribution
    double metropolis_unif::get_sigma() {
    return unif_sigma; // Return the current sigma for uniform distribution
    };

    void metropolis_unif::set_sigma(double x) {
    unif_sigma = x; // Set a new sigma value for uniform distribution
    };

    // Get and set methods for old coordinates in uniform distribution
    vec metropolis_unif::get_coord_old() {
    return coord_unif_old; // Return the previously stored coordinates
    };

    void metropolis_unif::set_coord_new(vec x) {
    coord_unif_old = x;
     // Update the old coordinates with new ones
    };

    // Implementation of the Metropolis algorithm for uniform distribution
    void metropolis_unif::metropolis_method(Waves_function* WF, vec coordinates) {

    if (WF == nullptr) {
    cerr << "Error: null WF pointer passed to initialize_acceptance!" << endl;
    return;
    }
    if (coordinates.n_elem != 3) {
    cerr << "Error: coordinates vector must have exactly 3 elements!" << endl;
    return;
    }

    // Propose new coordinates using uniform random distribution
    vec coordinate_proposed = RND.random_newcoord_unif(get_coord_old(), unif_sigma);
    

    
    //coordinate_proposed.print("coordinate del punto proposto: ");
    
    // Calculate wave function values at current and proposed coordinates
    double x = WF->measure(get_coord_old());
    //cout << "x: " << x << endl;
    double y = WF->measure(coordinate_proposed);
    //cout << "y: " << y << endl;
    double z = RND.Rannyu(0,1); // Generate a uniform random number in [0,1]

    // Metropolis acceptance criteria
    if (pow(y / x, 2) >= 1.) {
        set_coord_new(coordinate_proposed); // Accept new coordinates
        //cout << "coordinate accetate" << endl;
    }
    if (pow(y / x, 2) < 1. && z > pow(y / x, 2)) {
        set_coord_new(get_coord_old()); // Keep old coordinates
        //cout << "coordinate non accetate" << endl;
    }
    if (pow(y / x, 2) < 1. && z <= pow(y / x, 2)){
        set_coord_new(coordinate_proposed); // Accept new coordinates
        //cout << "coordinate accetate" << endl;
    }
    
    };

    void metropolis_unif::initialize_acceptance(Waves_function* WF, vec coordinates, int wf) {
    double counter = 0.;
    ofstream out_acc;
    
    out_acc.open("./OUTPUT/accettanza_unif_"+to_string(wf)+"_.dat");
    bool repeat = false;
    vec coordinates_true = coordinates;

    int step = 100000;
    double s = get_sigma(); // Ottieni il sigma iniziale
    double target_acceptance = 0.5; // Target di accettazione del 50%
    double tolerance = 0.0001; // Tolleranza per l'accettazione

    //cout << "prima sigma" << s << endl;
    out_acc << "Accettanza for uniform method is start:"<< endl << " sigma "<< " " << " accettanza " << endl;
    
    do {
        
        counter = 0; //reinitialize the cycle
        vec x = get_coord_old();
            //x.print("x primario ");
    
        for (int i = 0; i < step; i++) {
            vec x = get_coord_old();
            //x.print("x: ");
            metropolis_method(WF, x);
            vec y = get_coord_old();
            //y.print("y: ");
            // Confronto tra vettori, if they are the same means that the vector proposed is not accepted
            if (!arma::all(x == y)) {
                counter += 1;
            }
        }

        double acceptance_rate = counter / (double) step;
        //cout << "sigma: "<< get_sigma() << " " <<" accettanza del ciclo: " << acceptance_rate*100. << endl;
        out_acc << get_sigma() << " " << acceptance_rate*100 << endl;
         if (acceptance_rate > target_acceptance + tolerance) {
        s *= 1.0 + 0.1 * (acceptance_rate - target_acceptance); // Aumenta sigma proporzionalmente all'eccesso di accettazione
        }else if (acceptance_rate < target_acceptance - tolerance) {
        s *= 1.0 - 0.1 * (target_acceptance - acceptance_rate); // Diminuisci sigma proporzionalmente al deficit di accettazione
        } 
        if (acceptance_rate > target_acceptance - tolerance && acceptance_rate < target_acceptance + tolerance){
            repeat = true; // Fermati se il tasso è dentro il range accettabile
        }

        set_sigma(s); // Aggiorna sigma
        set_coord_new(coordinates_true);
        //cout << "sigma aggiornata: " << get_sigma() << endl; 

        } while (repeat == false);

    
    cout << endl << "sigma accettata per uniforme: " << get_sigma() << endl; //finished the acceptance fase I reset the initial coordinates
    out_acc << endl << "ACCETANZA FINISHED!"<< endl;
    out_acc.close();
};

void metropolis_unif:: print_coordinate_accepted(int wf,double pos){
    vec coordinate = get_coord_old();
    ostringstream oss_x;
    oss_x << setprecision(1) << fixed << pos;
    ofstream out_cor;
    double r = coordinate[0];      // raggio
    double theta = coordinate[1];  // angolo polare(in radianti)
    double phi = coordinate[2];    // angolo azimutale (in radianti)

    // Conversione da coordinate sferiche a cartesiane
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);
    
    out_cor.open("./OUTPUT/Coordinate_accettate_uniforme_"+to_string(wf)+"_"+oss_x.str()+"_.dat", ios:: app);
    out_cor << x << " "<< y << " "<< z << " " << endl;
    
    out_cor.close();
    
};


///////////////////////////////GAUSSSSSS////////////////////////////////////////////////

// Implementation of the Gaussian Metropolis algorithm
metropolis_gauss :: metropolis_gauss() {
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
}; // Constructor for Gaussian Metropolis

// Get and set methods for the sigma value in Gaussian distribution
double metropolis_gauss::get_sigma() {
    return gauss_sigma; // Return the current sigma for Gaussian distribution
};

void metropolis_gauss::set_sigma(double x) {
    gauss_sigma = x; // Set a new sigma value for Gaussian distribution
};

// Get and set methods for old coordinates in Gaussian distribution
vec metropolis_gauss::get_coord_old() {
    return coord_gauss_old; // Return the previously stored coordinates
};

void metropolis_gauss::set_coord_new(vec x) {
    for (int i = 0; i < 3; i++){
        coord_gauss_old[i] = x[i];
    }
     // Update the old coordinates with new ones
};

// Implementation of the Metropolis algorithm for Gaussian distribution
void metropolis_gauss::metropolis_method(Waves_function* WF, vec coordinates) {

    if (WF == nullptr) {
    cerr << "Error: null WF pointer passed to initialize_acceptance!" << endl;
    return;
    }
    if (coordinates.n_elem != 3) {
    cerr << "Error: coordinates vector must have exactly 3 elements!" << endl;
    return;
    }

    // Propose new coordinates using gauss random distribution
    vec coordinate_proposed = RND.random_newcoord_gauss(get_coord_old(), gauss_sigma);
    

    
    //coordinate_proposed.print("coordinate del punto proposto: ");
    
    // Calculate wave function values at current and proposed coordinates
    double x = WF->measure(get_coord_old());
    //cout << "x: " << x << endl;
    double y = WF->measure(coordinate_proposed);
    //cout << "y: " << y << endl;
    double z = RND.Rannyu(0,1); // Generate a gauss random number in [0,1]

    // Metropolis acceptance criteria
    if (pow(y / x, 2) >= 1.) {
        set_coord_new(coordinate_proposed); // Accept new coordinates
        //cout << "coordinate accetate" << endl;
    }
    if (pow(y / x, 2) < 1. && z > pow(y / x, 2)) {
        set_coord_new(get_coord_old()); // Keep old coordinates
        //cout << "coordinate non accetate" << endl;
    }
    if (pow(y / x, 2) < 1. && z <= pow(y / x, 2)){
        set_coord_new(coordinate_proposed); // Accept new coordinates
        //cout << "coordinate accetate" << endl;
    }
    
    };

void metropolis_gauss::initialize_acceptance(Waves_function* WF, vec coordinates,int wf) {
    double counter = 0.;
    bool repeat = false;
    ofstream out_acc;
    out_acc.open("./OUTPUT/accettanza_gauss_"+to_string(wf)+"_.dat");
    vec coordinates_true = coordinates;
    int step = 100000;
    double s = get_sigma(); // Ottieni il sigma iniziale
    double target_acceptance = 0.5; // Target di accettazione del 50%
    double tolerance = 0.0001; // Tolleranza per l'accettazione

    //cout << "prima sigma" << s << endl;
     out_acc << "Accettanza for gaussian method is start:"<< endl <<"sigma"<< " " << "accettanza" << endl;
        do {
        
        counter = 0; //reinitialize the cycle
    
        for (int i = 0; i < step; i++) {
            vec x = get_coord_old();
            //x.print("x: ");
            metropolis_method(WF, x);
            vec y = get_coord_old();
            //y.print("y: ");
            // Confronto tra vettori, if they are the same means that the vector proposed is not accepted
            if (!arma::all(x == y)) {
                counter += 1;
            }
            
        }

        

        double acceptance_rate = counter / (double) step;
        //cout << " sigma: "<< get_sigma() << " " << " accettanza del ciclo: " << acceptance_rate*100 << endl;
        out_acc << get_sigma() << " " << acceptance_rate*100 << endl;
        if (acceptance_rate > target_acceptance + tolerance) {
        s *= 1.0 + 0.1 * (acceptance_rate - target_acceptance); // Aumenta sigma proporzionalmente all'eccesso di accettazione
        }else if (acceptance_rate < target_acceptance - tolerance) {
        s *= 1.0 - 0.1 * (target_acceptance - acceptance_rate); // Diminuisci sigma proporzionalmente al deficit di accettazione
        } 
        if (acceptance_rate > target_acceptance - tolerance && acceptance_rate < target_acceptance + tolerance){
            repeat = true; // Fermati se il tasso è dentro il range accettabile
        }

        set_sigma(s); // Aggiorna sigma
        set_coord_new(coordinates_true);
         


    } while (repeat == false);

    set_coord_new(coordinates_true);
    cout << endl << "sigma accettata per gauss: " << get_sigma() << endl; //finished the acceptance fase I reset the initial coordinates
    out_acc << endl << "ACCETANZA FINISHED!"<< endl;
    out_acc.close();
};

void metropolis_gauss:: print_coordinate_accepted(int wf,double pos){
    vec coordinate = get_coord_old();
    ostringstream oss_x;
    oss_x << setprecision(1) << fixed << pos;
    ofstream out_cor;
    double r = coordinate[0];      // raggio
    double theta = coordinate[1];  // angolo polare(in radianti)
    double phi = coordinate[2];    // angolo azimutale (in radianti)

    // Conversione da coordinate sferiche a cartesiane
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);
    
    out_cor.open("./OUTPUT/Coordinate_accettate_uniforme_"+to_string(wf)+"_"+oss_x.str()+"_.dat", ios:: app);
    out_cor << x << " "<< y << " "<< z << " " << endl;
    
    out_cor.close();
    
};
//////////////////////////////////////////////////////////////////////////
