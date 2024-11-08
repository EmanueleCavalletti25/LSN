#include "tsp.h"

  tsp::~tsp() {
    
  }
// FUNCTION TO INITIALIZE THE PROBLEM
    void tsp :: initialize(int rank){     //function to initialize the system: shape and other parameters
        int p1, p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12; // Read from ./INPUT/Primes a pair of numbers to be used to initialize the RNG
        ifstream Primes("./INPUT/Primes");
        string line;
        if (rank == 0){                   //parallelized the lecture of Primes
        Primes >> p1 >> p2 ;
        }
        if (rank == 1){ 
        getline(Primes,line);                
        Primes >> p3 >> p4 ;
        }
        if (rank == 2){
        getline(Primes,line);    
        getline(Primes,line);    
        Primes >> p5 >> p6 ;
        }
        if (rank == 3){
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        Primes >> p7 >> p8 ;
        }
        if (rank == 4){
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        Primes >> p9 >> p10 ;
        }
        if (rank == 5){
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        getline(Primes,line);    
        Primes >> p11 >> p12 ;
        }
        Primes.close();         
        int seed[4]; // Read the seed of the RNG      //SEED IS EQUAL FOR EVERYONE
        ifstream Seed("./INPUT/seed.in");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        if (rank == 0){_rnd.SetRandom(seed,p1,p2);} //parallelized the initialize for random
        if (rank == 1){_rnd.SetRandom(seed,p3,p4);}
        if (rank == 2){_rnd.SetRandom(seed,p5,p6);}
        if (rank == 3){_rnd.SetRandom(seed,p7,p8);}
        if (rank == 4){_rnd.SetRandom(seed,p9,p10);}
        if (rank == 5){_rnd.SetRandom(seed,p11,p12);}
        
        string filename_input;
        filename_input = "./INPUT/input.dat";
        
        //BEGIN LECTURE
        ifstream inputFile(filename_input);
        if (!inputFile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename_input <<  endl;
        exit(-1);  // Esce dal programma con codice di errore
    }
        ifstream input(filename_input); // Start reading ./INPUT/input.dat
        ofstream coutf;
        string output_f = "./OUTPUT/output.dat";
        coutf.open(output_f);
        string property;
        while ( !input.eof() ){
            input >> property;
            if( property == "SHAPE" ){
            input >> _shape;
        if(_shape > 1){
        cerr << "PROBLEM: unknown _shape type" << endl;
        exit(EXIT_FAILURE);
    }
      if(_shape == 0)      coutf << "TSP WITH CITIES ON CIRCONFERENCE"  << endl;
      else if(_shape == 1) coutf << "TSP WITH CITIES IN SQUARE"         << endl;
     
    } else if( property == "RADIUS" ){
      input >> _r;
      coutf << "RADIUS = " << _r << endl;
    } else if( property == "N_CITY" ){
      input >> _ncity;
      m_dist = arma::zeros(_ncity,_ncity);
      coutf << "NUMBER OF CITY = " << _ncity << endl;
      _cityes.set_size(_ncity);              //set the size of population field
      vec pos_c(2,fill::zeros);
      vec pos_s(2,fill::zeros);

      for(int i=0; i<_ncity; i++){ 
        if (rank == 0){
          if (_shape == 0){pos_c = _rnd.Circonference((double)_r);}
          if (_shape == 1){pos_s = _rnd.Square((double)_r);}
        }
        MPI_Bcast(pos_c.memptr(), pos_c.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(pos_s.memptr(), pos_s.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        _cityes(i).initialize();
        
          if (_shape == 0){_cityes(i).set_pos(pos_c);};
          if (_shape == 1){_cityes(i).set_pos(pos_s);}
        
        
      }
    
    } else if( property == "POPULATION" ){
      input >> _npop;
      coutf << "POPULATION = " << _npop << endl;
    
      
      m_pop = arma::zeros(_npop,_ncity-1); //matrix of population: n_city-1; for better shuffle.
      
      
      
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  
  input.close();
  coutf << "System initialized!" << endl;
  coutf.close();

  /*for(int i = 0 ; i <_ncity; i++){
    cout <<"città: "<< i << " x: " <<_cityes(i).get_pos()[0] << " y:" <<_cityes(i).get_pos()[1] << endl;
    cout << (sqrt(pow(_cityes(i).get_pos()[0],2)+pow(_cityes(i).get_pos()[1],2))) << endl;
  }*/

  //write the output_loss_function
  
  return;
}; 
    
    void tsp :: initialize_prop(int rank){
    
  ofstream out_m;
  ofstream out_p;
  if (_shape == 0){out_m.open("./OUTPUT/Loss_media_circunference_best_"+std::to_string(rank)+"_.dat");}
  if (_shape == 1){out_m.open("./OUTPUT/Loss_media_square_best_"+std::to_string(rank)+"_.dat");}
  if (_shape == 0){out_p.open("./OUTPUT/Loss_media_circunference_best_"+std::to_string(rank)+"_single.dat");}
  if (_shape == 1){out_p.open("./OUTPUT/Loss_media_square_best_"+std::to_string(rank)+"_single.dat");}

  out_m << "BEST LOSS FUNCTION MEAN and DV_STD"<< endl;
  out_p << "BEST LOSS FUNCTION" << endl;

  out_m.close();
  out_p.close();
      //fill the distance matrix: (all the possible distance matrix);
      for(int i = 0; i < _ncity; i++){
        for(int j = 0; j <_ncity; j++){
          m_dist(i,j) = distance(_cityes(i).get_pos(),_cityes(j).get_pos());
          
        }
      }
      //m_dist.print();

      //fill the replica matrix: fix the first one and shuffle the other;
      for(int i = 0; i < _npop; i++){
        for(int j = 0; j <_ncity-1; j++){
          m_pop(i,j) = (int)j+1;
        }
        m_pop.row(i) = arma::shuffle(m_pop.row(i));
        check_vec(m_pop.row(i).t());
      }
     // m_pop.print();



    }; //creation of repliche in the system and matrix distance
    

    //FUNCTION OF GENETIC ALGORITHM//////////////////////////////////////////////////////////////////////////////

    //Calculation of Loss function:: remember that m_pop is (100,33) and not (100,34) because, the first city is fixed.
    double tsp:: Fitness_f_1(vec x){
      double Loss = abs(m_dist(0,x[1]));
      for (int i = 0; i < x.size()-1; i++){
          Loss += abs(m_dist(x[i],x[i+1]));
      }
      Loss += m_dist(x[x.size()-1],0);
      return 1./Loss;
    };


    void tsp:: pair_permutation(vec& x){
        int m = floor(_rnd.Rannyu(0,1)*((double)x.size()));
        int n = floor(_rnd.Rannyu(0,1)*((double)x.size()));
        swap(x,m,n);
      };

    void tsp:: shift_n(vec& x){
      int i = floor(_rnd.Rannyu(0,1)*((double)x.size()));   //start from
      //cout << " i:  " << i << endl; 
      int m = floor(_rnd.Rannyu(0,1)*((double)(x.size()-i-1)))+1; //contiguos cell
      //cout << " m:  " << m << endl; 
      int shift = floor(_rnd.Rannyu(0,1)*((double)(x.size()-i-(m-1)))); //shift
      //out << " shift:  " << shift << endl; 

    for(int s = 0 ; s < shift; s++){
      for (int j = 0; j < m; j++ ){
        swap(x, i + m - j + s , i +(m-1-j)+s);
        }
      }
    };
  
  void tsp:: swap_m_city(vec& x) {
    int i = floor(_rnd.Rannyu() * ((double)x.size()));  // start index 'i'

    // Generate 'm', the length of contiguous block to swap, ensuring 1 <= m < (x.size()-1)/2
    int m = floor(_rnd.Rannyu() * (double)((x.size() - 1) / 2.)) + 1;

    // Ensure the first set of cities is not too long and can fit within the array without overlap
    while (i + m >= x.size() || (i + m > (double)((x.size() - 1) / 2.) && i - m <= 0)) {
        i = abs(floor(_rnd.Rannyu() * ((double)x.size())));  // regenerate 'i' if condition not met
        m = abs(floor(_rnd.Rannyu() * (double)((x.size() - 1) / 2.))) + 1;  // regenerate 'm' if needed
        //cout << "i " << i << endl;
        //cout << "mi blocco quaaa" << endl;
    }

    int f = 0;
    bool ok = false;
    //cout << "i: " << i << endl;
    //cout << "m: " << m << endl;

    // Find a valid second block 'f' that doesn't overlap with the first block (i, i+m)
    do {
        f = abs(floor(_rnd.Rannyu() * ((double)x.size())));  // another start index 'f'

        if (f > i && f < i + (m-1)) {  // 'f' cannot be between 'i' and 'i+m'
            //cout << "entro cond 1" << endl;
            ok = false;
        } else if (f > i + (m-1) && f + (m-1) < x.size()) {  // valid if 'f+m' fits within the array
            //cout << "entro cond 2" << endl;
            ok = true;
        } else if (f > i + (m-1) && f + (m-1) >= x.size()) {  // invalid if 'f+m' exceeds the array
            //cout << "entro cond 3" << endl;
            ok = false;
        } else if (f < i && f + (m-1)< i) {  // valid if 'f+m' is before 'i'
            //cout << "entro cond 4" << endl;
            ok = true;
        } else if (f < i && f + (m-1) > i) {  // invalid if 'f+m' overlaps with 'i'
            //cout << "entro cond 5" << endl;
            ok = false;
        } else if (f == i || f + (m-1) == i + (m-1) || f + (m-1) == i) {  // invalid if blocks overlap exactly
            //cout << "entro cond 6" << endl;
            ok = false;
        } 

        // Ensure if m == (x.size() - 1) / 2 and i == 0, we avoid potential out-of-bound issues
        if (((m-1) == (double)(x.size() - 1) / 2.) && i == 0) {
            //cout << "entro cond 7" << endl;
            i = abs(floor(_rnd.Rannyu() * ((double)x.size())));  // regenerate 'i' if condition not met
            m = abs(floor(_rnd.Rannyu() * (double)((x.size() - 1) / 2.))) + 1;  // regenerate 'm' if needed
            ok = false;
        }

    } while (!ok);

    //cout << "f: " << f << endl;

    // Perform the swap between the two blocks of cities
    for (int l = 0; l < m; l++) {
        swap(x, i + l, f + l);
    }
}

void tsp:: swap_back(vec& x){
      int i = floor(_rnd.Rannyu(0,1)*((double)x.size()));   //start from i
      //cout << " i:  " << i << endl; 
      int m = floor(_rnd.Rannyu(0,1)*((double)(x.size()-i-1)))+1; //contiguos cell
      //cout << " m:  " << m << endl; 
      if (m % 2 ){
        for (int s = 0; s <= m/2; s++){
          swap(x, i+s , i+(m-1)-s);
        }
      }else{
        for (int s = 0; s <= (m-1)/2; s++){
          swap(x,i+s,i+(m-1)-s);
        }
      }
    };

    vec tsp::cumulative_fit( mat popul , double pot) {
    int dim = popul.n_rows; // Population
    vec fit(dim, fill::zeros);
    
    // Pre-calculation of total fitness
    double sum_fit = 0.;
    for (int i = 0; i < dim; ++i) {
        sum_fit += pow(Fitness_f_1(vec(popul.row(i).t())),pot);
    }

    if (sum_fit == 0.0) {
    
    // if fitness is zero, give back an array of uniform probability uniformi to avoid dividing for zero.
        fit.fill(1.0 / dim);
        return cumsum(fit); // Give the cumulative sum.
    }

    // Calcolo dei valori di fitness cumulativo normalizzato
    fit[0] = pow(Fitness_f_1(vec(popul.row(0).t())),pot) / sum_fit;
    for (int i = 1; i < dim; ++i) {
        fit[i] = fit[i-1] + pow(Fitness_f_1(vec(popul.row(i).t())),pot) / sum_fit;
    }

    return fit;
}

    vec tsp:: Loss_vector(mat pop){
      double dim = (double)pop.size()/(double)(_ncity-1);
      vec loss (dim,fill::zeros);

      for (int i = 0; i < dim; i++){
        loss[i] = 1./Fitness_f_1(vec(pop.row(i).t()));
      }
      
      return loss;

    };

   int tsp:: select_idx_fit(vec cum_fit){
    double x = _rnd.Rannyu(0,1);
     for (int i = 0; i < cum_fit.size(); ++i) {
        if (x < cum_fit[i]) {
          return i;  // give the index of x taken;
         }
      }
    
      return cum_fit.size() - 1;  // nothing is taken (impossible); give last element
    
  };
      
void tsp:: mutation (vec& v){
        double x = 0.;
        int counter = 1;
        while (counter < 4){
        x = _rnd.Rannyu();
        if (counter == 1){
        //cout << "sto facendo pair permutation "<< endl;
          if(x < p_m){
            pair_permutation(v);
            check_vec(v);           //necessity to check if there are all the city in the vector;
            counter ++;
            //cout << "fatto pair permutation "<< endl;
            
          }else{
            //cout << " non fatto pair permutation "<< endl;
           counter++;}
        }
        x = _rnd.Rannyu();
        if (counter == 2){
        //cout << "sto facendo shift n "<< endl;
        if(x < p_m){
            shift_n(v);
            check_vec(v);
            counter ++;
            //cout << "fatto shift n "<< endl;
           
          }else{
            //cout << " non fatto shift n "<< endl;
            counter++;}
        }
        x = _rnd.Rannyu();
          if (counter == 3){
          //cout << "sto facendo m city "<< endl;
          if(x < p_m){
            swap_m_city(v);
            check_vec(v);
            counter ++;
            //cout << "fatto m city "<< endl;
            
          }else{
            //cout << "non fatto m city "<< endl;
            counter++;}
        }
        x = _rnd.Rannyu();
          if (counter == 4){
          //cout << "sto facendo swap back "<< endl;
          if(x < p_m){
            swap_back(v);
            check_vec(v);
            counter ++;
            //cout << "fatto swap back"<< endl;
           
          }else{
            //cout << "non fatto swap back"<< endl;
            counter++;}
        }
      }
        

    };

    void tsp:: children_generator(vec& a, vec& b){
    if (a.size()!=b.size()){
      cerr << "vector haven't same dimension: no child_gereator"<< endl;
      exit(-1);
    }
       
       int cut = std::floor(_rnd.Rannyu(0,1)*a.size()); //index cut
       //cout << "cut at " << cut << endl;
      //cout << "cut non floor" << ds << endl;
       
       vec child_a(a.size(),fill::zeros);
       vec child_b(b.size(),fill::zeros);

       for(int j = 0; j <= cut; j++){
        child_a[j]=a[j];
        child_b[j]=b[j];
       }
       //find the lost (index) city for A
       vec idx_a_lost(a.size()-(cut+1), fill::zeros);
       int counter = 0;
       arma::uvec indices_find_a;
       
      for (int l = 1 ; l < _ncity; l++){
        indices_find_a = arma::find(child_a == l);
        if(indices_find_a.is_empty()){
          for (int i = 0; i < b.size();i++){
            if (b[i] == l){
              idx_a_lost[counter] = i;
              counter++;
            }
          }
        }
      }
      
      idx_a_lost=arma::sort(idx_a_lost);
      for (int m = 0 ; m < idx_a_lost.size(); m++){
        child_a[cut+1+m] = b[idx_a_lost[m]];
      }

      //find the lost (index) city for B
      vec idx_b_lost(b.size()-(cut+1), fill::zeros);
      counter = 0;
      arma::uvec indices_find_b;
       
      for (int l = 1 ; l < _ncity; l++){
        indices_find_b = arma::find(child_b == l);
        if(indices_find_b.is_empty()){
          for (int i = 0; i < a.size();i++){
            if (a[i] == l){
              idx_b_lost[counter] = i;
              counter++;
            }
          }
        }
      }
      
      idx_b_lost=arma::sort(idx_b_lost);
      for (int m = 0 ; m < idx_b_lost.size(); m++){
        child_b[cut+1+m] = a[idx_b_lost[m]];
      }

      a = child_a;
      b = child_b;
    };

    void tsp :: crossover (vec& v, vec& w){
      double s = _rnd.Rannyu(0,1);
      if (s < p_c){
        children_generator(v,w);
        //cout << "crossover happened"<< endl;
        check_vec(v);
        check_vec(w);
        }
    };

    
    //FUNCTION TO WRITE OUTPUT

      double tsp::write_loss_function(int rank){
      vec loss_vectore = Loss_vector(m_pop);
      double media_50 = Loss_best_media(loss_vectore);
      vec v = arma::sort(loss_vectore);
      double std_deviation = arma::stddev(v.head(_npop/2));

      ofstream out_m;
      if (_shape == 0){out_m.open("./OUTPUT/Loss_media_circunference_best_"+std::to_string(rank)+"_.dat", ios::app);}
      if (_shape == 1){out_m.open("./OUTPUT/Loss_media_square_best_"+std::to_string(rank)+"_.dat", ios::app);}

      out_m << media_50 << " "<< std_deviation << endl;
      
      out_m.close();

      return media_50;
      
    };

      double tsp:: write_loss_function_single(int rank){
      vec loss_vectore = Loss_vector(m_pop);

      vec v = arma::sort(loss_vectore);
      ofstream out_m;
      if (_shape == 0){out_m.open("./OUTPUT/Loss_media_circunference_best_"+std::to_string(rank)+"_single.dat", ios::app);}
      if (_shape == 1){out_m.open("./OUTPUT/Loss_media_square_best_"+std::to_string(rank)+"_single.dat", ios::app);}

      out_m << v[0] << endl;

      
      
      out_m.close();
      return v[0];

    };

    

    void tsp::write_best_path(mat pop,int rank){
      vec loss_vec = Loss_vector(m_pop); 
      ofstream out_mat;
      
      if (_shape == 0){out_mat.open("./OUTPUT/best_path_circunference_"+std::to_string(rank)+"_.dat");}
      if (_shape == 1){out_mat.open("./OUTPUT/best_path_square_"+std::to_string(rank)+"_.dat");}
      
      out_mat << "X:" << setw(12) << "Y:" << endl;
      
      int best_idx = Loss_best_idx(loss_vec);
    
      vec best_v = pop.row(best_idx).t();
      
      out_mat << _cityes(0).get_pos()[0] << setw(12) << _cityes(0).get_pos()[1] << endl; 
      for (int i = 0 ; i < best_v.size();i++ ){
        int order_city = best_v[i];
        out_mat << _cityes(order_city).get_pos()[0] << setw(12) << _cityes(order_city).get_pos()[1] << endl; 
      }
      
      out_mat << _cityes(0).get_pos()[0] << setw(12) << _cityes(0).get_pos()[1] << endl; 
      out_mat.close();
      };



      void tsp:: write_core_best(int rank,int min_rank, double global_min){
          if (rank==0){
          ofstream out_best;
          if(_shape==0)out_best.open("./OUTPUT/output_best_circunference.dat");
          if(_shape==1)out_best.open("./OUTPUT/output_best_square.dat");
          out_best << "il rank migliore è: " << min_rank << " con loss_function: " << global_min << endl;

           out_best.close();
         }
      };


    //FUNCTION TO ELABORATE DATA///////////////////////////////////////////////////////////////////////////
    double tsp :: distance(vec v, vec w){
        double distance = 0.;
        if(v.size() != w.size()){
            cerr << "DIMENSIONS OF VECTOR DON'T CORRESPOND";
            exit(-1);
        }
        for (int i = 0; i < v.size(); i++){
             distance += pow(v[i] - w[i],2);
            
        }
        //cout << distance << endl;
        distance = sqrt(distance);
        if (_shape == 0){
        if(distance > (double)2*_r){
            cerr << "DIMENSIONS OF VECTOR DON'T CORRESPOND";
            exit(-1);
        }
        }

        return distance;
    };
        //function to check if there is all the city in the row of mpop;
    void tsp :: check_vec(vec a){
      int counter = 1;
      int not_f = 0;
      bool find_el = false;
      while (counter < _ncity){
        for (int i = 0; i < a.size(); i++){
          
          if (a[i] == counter){
           //cout <<"trovato: "<< a[i] << endl; 
            find_el = true;
          } 
          
        }
        if (find_el != true)
          not_f = counter;
        counter++;
        find_el = false;
        }
      if (not_f != 0){
        
        cerr << "ERROR, THERE IS NO ALL CITY IN THE VEC" << endl;
        exit (-1);
      }

    };
    
    void tsp :: swap(vec& v ,int id_x ,int id_y){
      int t = v[id_x];
      v[id_x] = v[id_y];
      v[id_y] = t;
      
    };

    double tsp :: Loss_best_media(vec v){
      if (v.size()< 50){
        cerr << "problema nella loss size"<< endl;
        exit(-1);
      }
      v = arma::sort(v); //i primi 50, ne faccio la media
      double media = 0.;
      double counter = 0.;

      for (int i = 0; i < _npop/2 ; i++){
        media += v[i];  
        counter++;
      }

      return media/(double)counter;

    };

    

    int tsp:: Loss_best_idx(vec v){
       vec v_copy = arma::sort(v); //sort the vector
       double x = v_copy[0];
       arma::uvec best;

       best = arma::find(v==x);
       //cout << "better path idx" << best[0]<< endl;
       //v.print("LOSS_VECTOR");
       if (best.is_empty()){
        cerr << "ERROR, IN FIND BEST FIT" << endl;
        exit(-1);
       }
       return best[0];
    };

    int tsp:: find_sender_idk(int rank, vec rank_n){
        int target_value = rank; // Value to find
        int index = -1; // Initialize index to -1 to indicate not found

      for (int i = 0; i < rank_n.n_elem; ++i) {
          if (rank_n(i) == target_value) {
              index = i; // Store the index if found
              break; // Exit loop once found
        }
      }

      return index;
    
    }

        /*void tsp::migration(int rank, vec ranks, vec shuf_rank) {
    vec Loss_idx_core(ranks.size(), fill::zeros);
    
    // Calcola l'indice della miglior soluzione per ogni rank
    vec Loss_vector_best = Loss_vector(m_pop);
    double best_idx = Loss_best_idx(Loss_vector_best);
    
    // Condividi il migliore con gli altri rank
   
    
    // A questo punto, ogni processo ha il miglior indice del rank shuffle corrispondente
    Loss_idx_core.print("vettore: ");

}*/
    
    
    
    /* void tsp:: m(int rank,int rank_1, int rank_2, int row_to_send, int row_to_receive){
       int cols = m_pop.n_cols;    // obtain number of column
       arma::rowvec buffer(cols);   // temporary vector for exchange

      // Invio della riga da rank a other_rank
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) == rank_1) {
        // Invia la riga row_to_send al processo other_rank
        MPI_Send(vec(m_pop.row(row_to_send)).memptr(), cols, MPI_INT, rank_2, 0, MPI_COMM_WORLD);
        // Ricevi la riga dal processo other_rank
        MPI_Recv(buffer.memptr(), cols, MPI_INT, rank_2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Sostituisci la riga row_to_receive con i dati ricevuti
        m_pop.row(row_to_receive) = buffer;
        //cout << "ho scambiato la best di " << rank_1 << "con la best di " << rank_2 << endl;
    }else if (rank_2) {
        // Ricevi la riga dal processo rank
        MPI_Recv(buffer.memptr(), cols, MPI_INT, rank_1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Invia la riga row_to_send al processo 0
        MPI_Send(vec(m_pop.row(row_to_send)).memptr(), cols, MPI_INT, rank_1, 0, MPI_COMM_WORLD);
        // Sostituisci la riga row_to_receive con i dati ricevuti
        m_pop.row(row_to_receive) = buffer;

        cout << "ho scambiato la best di " << rank_1 << "con la best di " << rank_2 << endl;
    }
    };*/

    
    //FUNCTION GET  SET/////////////////////
    mat tsp :: get_mat_pop(){
      return m_pop;
    };

    void tsp :: set_mat_pop(mat pop_new){
      m_pop = pop_new;
    };

    void tsp::set_row_pop(vec new_row,int idx_row){ //function to modify the row of matrix population
      check_vec(new_row);
      if (new_row.size()!= vec(m_pop.row(idx_row).t()).size()){
        cerr << "Problem in dimension of new row" << endl;
        exit(-1);
      }
      m_pop.row(idx_row)= new_row.t();
    };

    vec tsp:: get_row_pop(int idx_row){
      return vec(m_pop.row(idx_row).t());
    };      
    
    double tsp:: get_pc(){return p_c;}  //get and set new prob of croossover
    void tsp:: set_pc(double p_new){p_c =p_new;};
    
    double tsp::get_pm(){return p_m;}   //get and set new prob of mutation
    void tsp::set_pm(double p_new){p_m =p_new;};