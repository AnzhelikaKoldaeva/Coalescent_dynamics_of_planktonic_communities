//
// Created by paula on 17/11/15.
//


#include "include.h"

using namespace std;

string MODEL_NAME="";
float MAX_MBITS=100;
bool descending_order (int i,int j) { return (i>j); }
struct individuals_in_area{
    list<int> individual_list;
};
struct tree_individual{
    list<tree_individual>::iterator parent;
    unsigned int own_dir,label, intermediate_births,offspring_number;
    bool active,print,mutation,species_assignated;
    double t_birth;
    list<int>::iterator species_abundance;
    float x_initial,y_initial;

    tree_individual(){
        print=false;
        mutation=false;
        species_assignated=false;
        intermediate_births=0;
        active=false;
        offspring_number=0;
        t_birth=0;
    }

    void initialize(){
	//cout<<"Initialize 34"<<endl;
        mutation=false;
        species_assignated=false;
    }
};
struct ancestor_individual{

    unsigned int alive_offspring;
    list <tree_individual> offspring;

    ancestor_individual(){
        alive_offspring=1;
    };
    ~ancestor_individual(){
        //cout << "removing ancestor with " << offspring.size() << " offsprings" << endl;
    };
};
struct dynamics_individual{

    list<ancestor_individual>::iterator ancestor;
    list <tree_individual>::iterator offspring_it;

    virtual void initial_set(gsl_rng *r,double force, double *p= nullptr){};
};
struct flux_individual:public dynamics_individual{
    float x,y;
    float fieldx,fieldy;
    float fieldoldx,fieldoldy;
    int number_cells_x;
    int number_cells_y;
    double Lx,Ly;
    list<individuals_in_area>::iterator current_area;


    void initial_set(gsl_rng *r, ifstream *initial_spatial_distribution_data,double force=0, double *p= nullptr, double Lx_ini=0, double Ly_ini=0, double xini=0,double yini=0){
        //double Ly=p[5];//7.5;
        //double Lx=p[5];//7.5;
        double x_origin=number_cells_x*Lx/2.-Lx/2.;
        double y_origin=number_cells_y*Ly/2.-Ly/2.;

#ifdef GET_INITIAL_SPATIAL_DISTRIBUTION
        double get_trash;
        (*initial_spatial_distribution_data)>>x>>y>>get_trash;
#else
	//double Lx= 7.5;
	//double Ly = 7.5;
        x=gsl_rng_uniform(r)*Lx;    //xini+(gsl_rng_uniform(r)*Lx_ini-Lx_ini/2.)+x_origin;
        y=gsl_rng_uniform(r)*Ly;   //yini+(gsl_rng_uniform(r)*Ly_ini-Ly_ini/2.)+y_origin;
#endif
        advectx(force,0,0,p);
        advecty(force,0,0,p);
        fieldoldx=fieldx;
        fieldoldy=fieldy;
    }
    void advectx(double force=0, double t=0, double Ly_particular=1, double *p= nullptr, double Kxi=0, double Kyi=0, double ti=0){

        #ifdef ADVECTION
        double Ly=p[5];
	      double Lx=p[5];

        double flux_y=((y+Kyi*ti)-floor((y+Kyi*ti)/Ly)*Ly)-Ly/2.;
        double flux_x=((x+Kxi*ti)-floor((x+Kxi*ti)/Lx)*Lx)-Lx/2.;

      	double B0 = p[0]; //1.2;
      	double eps = p[1]; //0.3;
      	double Fi = p[2]; //M_PI/2;
      	double w = p[3]; //0.4;
      	double c = p[4]; //0.12;
        double L = p[5]; //7.5;
        double k = p[6]; //2*M_PI/L;
        double resc_const =  p[7];


      	double Bt = B0 + eps*cos(w*(t+ti)+Fi);
      	double tanh_arg = (flux_y - Bt*cos(k*flux_x))/sqrt(1+pow(k*Bt*sin(k*flux_x),2));
      	double expr = 1+pow(k*Bt*sin(k*flux_x),2);

        fieldx=((1-pow(tanh(tanh_arg),2))/sqrt(expr)-c)*sqrt(resc_const);
        #else
        fieldx=0;
        #endif
	     }
    void advecty(double force=0, double t=0, double Lx_particular=1, double *p= nullptr, double Kxi=0, double Kyi=0, double ti=0){

        #ifdef ADVECTION
        double Ly=p[5];
	      double Lx=p[5];

        double flux_y=((y+Kyi*ti)-floor((y+Kyi*ti)/Ly)*Ly)-Ly/2.;
        double flux_x=((x+Kxi*ti)-floor((x+Kxi*ti)/Lx)*Lx)-Lx/2.;

        double B0 = p[0]; //1.2;
      	double eps = p[1]; //0.3;
      	double Fi = p[2]; //M_PI/2;
      	double w = p[3]; //0.4;
      	double c = p[4]; //0.12;
        double L = p[5]; //7.5;
        double k = p[6]; //2*M_PI/L;
        double resc_const =  p[7];

        double Bt = B0 + eps*cos(w*(t+ti)+Fi);
      	double tanh_arg = (flux_y - Bt*cos(k*flux_x))/sqrt(1+pow(k*Bt*sin(k*flux_x),2));
      	double expr = 1+pow(k*Bt*sin(k*flux_x),2);

        double fieldy= -(1-pow(tanh(tanh_arg),2))*k*Bt*sin(k*flux_x)*(expr-k*k*Bt*cos(k*flux_x)*(flux_y-Bt*cos(k*flux_x)))/(sqrt(expr)*expr)*sqrt(resc_const);

        #else
        fieldy=0;
        #endif
    }
};
struct species_active_individual{
    list <tree_individual>::iterator individual_it;
};
float tree_individual_struct_size=sizeof(tree_individual)*1e-6;


///GENERIC TREE CLASS
class tree{

private:
    bool previous_print;
    unsigned int MAX_NUMBER_OF_INDIVIDUALS,TOTAL_TRIALS;
    vector <double> SR_histogram,SAD_histogram;
    fstream SAD,SR,generation_tree_ACTIVES,generation_tree_FAMILIES,generation_tree_PROPERTIES;
    list <ancestor_individual>::iterator it_ancestor;

    int **common_species;
    fstream SAD_subsample;
    bool **checked_species;


    virtual void initialize_network(){}
    virtual void initialize(gsl_rng * random_variable, double *p= nullptr){

        tree_individual new_offspring_individual;
        flux_individual new_individual;
        ancestor_individual new_ancestor;
        last_label=0;
        final_spatial_data.open(final_spatial_name, ios::in);
        for(int i=0;i<NNN;i++){

            ancestors_list.push_back(new_ancestor);

            //The first individual is "0": parent of ancestors
            new_offspring_individual.label=0;
            ancestors_list.back().offspring.push_back(new_offspring_individual);
            ancestors_list.back().offspring.back().parent=--ancestors_list.back().offspring.end();
            last_label++;
            new_offspring_individual.parent=--ancestors_list.back().offspring.end();
            new_offspring_individual.label=last_label;

            new_individual.initial_set(random_variable,&final_spatial_data,force,p);
            new_offspring_individual.x_initial=new_individual.x;
            new_offspring_individual.y_initial=new_individual.y;
            ancestors_list.back().offspring.push_back(new_offspring_individual);

            new_individual.ancestor=--ancestors_list.end();
            new_individual.offspring_it=--ancestors_list.back().offspring.end();
            active_individuals_vector.push_back(new_individual);
        }
        final_spatial_data.close();

        initialize_network();
        before_individuals=0;
        tnext=0;
    }
    void eliminate_dead_branches(){

        bool continue_while;

        tree_individual *it_offspring_aux;

        for(it_ancestor=ancestors_list.begin();it_ancestor!=ancestors_list.end();++it_ancestor){
            for(it_offspring1=(*it_ancestor).offspring.begin();it_offspring1!=(*it_ancestor).offspring.end();it_offspring1++){
                (*it_offspring1).print=false;
            }
        }
        for(int i=0;i<active_individuals_vector.size();i++){

            it_ancestor=active_individuals_vector[0].ancestor;
            it_offspring1=active_individuals_vector[i].offspring_it;
            if(!((*it_offspring1).print)){
                (*it_offspring1).print=true;

                it_offspring_aux=&(*it_offspring1);
                it_offspring2=(*it_offspring_aux).parent;

                previous_print=false;
                continue_while=true;
                while((continue_while)&&(!previous_print)){
                    previous_print=(*it_offspring2).print;
                    (*it_offspring2).print=true;
                    if((*it_offspring2).label==0) continue_while=false;
                    it_offspring2=(*it_offspring2).parent;
                }
            }
        }
        after_individuals=0;
        for(it_ancestor=ancestors_list.begin();it_ancestor!=ancestors_list.end();++it_ancestor){
            it_offspring1=(*it_ancestor).offspring.begin();
            while(it_offspring1!=(*it_ancestor).offspring.end()){
                if(!(*it_offspring1).print) it_offspring1=(*it_ancestor).offspring.erase(it_offspring1);
                else {
                    it_offspring1++;
                    after_individuals++;
                }
            }
        }
        before_individuals=after_individuals;
    };

protected:
    time_t itime,ftime,TOTAL_itime,TOTAL_ftime;
    unsigned int before_individuals,after_individuals,random_seed,NNN,NNN_forward,last_label;
    double mutation_rate,force,Lx,Ly,dt,tnext;
    char SAD_name[5000],SAD_subsample_name[1280],SR_name[5000],generation_tree_name_ACTIVES[5000],generation_tree_name_FAMILIES[5000],generation_tree_name_PROPERTIES[5000],spacetime_name[5000],fname[5000];
    list <tree_individual>::iterator offspring_begin_it,offspring_end_it,it_offspring1,it_offspring2;
    string model;
    int N_sample,S_sample;
public:
    vector <flux_individual> active_individuals_vector;
    double t;
    int number_of_subareas;
    int number_of_occupied_subareas;
    struct subsample{
        int number_species;
        int number_individuals;
        vector <int> species_abundance_vector;
        bool occupied=false;
    };

    ifstream final_spatial_data;
    char final_spatial_name[300];

    void set_initial_spatial_data(double D_factor,int my_network){
        sprintf(final_spatial_name,"Final_spatial_distribution_Df%g_Mut%g_N%d_Seed%d_Network%d.dat",D_factor,mutation_rate,NNN,random_seed,my_network);
    };

    double mean_number_species, mean_number_species_subarea,mean_beta_diversity,mean_number_individuals_subarea;
    vector <double> SAD_histogram_subsamples_mean;
    vector <double> mean_heterozigosity_diff_pairs;
    vector <double> mean_heterozigosity_all_pairs;
    vector <double> mean_alpha;
    int number_of_spatial_radius;
    vector<double> coalescences_tdistribution,mutations_tdistribution,events_tdistribution;
    long unsigned int tdistribution_points;
    long unsigned int coalescences_total,mutations_total,current_mutations, current_coalescents;
    double sum_coal_prob;
    long unsigned int curr_no_coal_no_mut;
    double tdmax;
    double mean_final_individuals_distance;
    double mean_final_individuals_distance_all;
    double min_final_individuals_distance_all;
    double max_final_individuals_distance_all;

protected:
    flux_individual new_active_individual;
    list <ancestor_individual> ancestors_list;

    tree(int N,double mutprob,int r,int number_of_spatial_radius_in=0,double Lx_in=-1,double Ly_in=-1,double tst_in=1000):NNN(N),mutation_rate(mutprob),random_seed(r),number_of_spatial_radius(number_of_spatial_radius_in),Lx(Lx_in),Ly(Ly_in){
        printed_networks=0;
        network=0;
        TOTAL_TRIALS=0;
        last_label=0;
        force=0;

        MAX_NUMBER_OF_INDIVIDUALS=NNN+1;
        SAD_histogram.assign(MAX_NUMBER_OF_INDIVIDUALS,0);
        itime=time(NULL);
        TOTAL_itime=time(NULL);
        model=MODEL_NAME;

        mean_heterozigosity_diff_pairs.assign(number_of_spatial_radius,0);
        mean_alpha.assign(number_of_spatial_radius,0);
        mean_heterozigosity_all_pairs.assign(number_of_spatial_radius,0);
        SAD_histogram_subsamples_mean.assign(MAX_NUMBER_OF_INDIVIDUALS,0);
        tdmax=1e4;
        tdistribution_points=int(tdmax)+1;
        coalescences_tdistribution.assign(tdistribution_points,0);
        mutations_tdistribution.assign(tdistribution_points,0);
        events_tdistribution.assign(tdistribution_points,0);
        coalescences_total=0;
        mutations_total=0;
        mean_final_individuals_distance=0;
        mean_final_individuals_distance_all=0;
        min_final_individuals_distance_all=0;
        max_final_individuals_distance_all=0;


    }
    virtual void one_step_dynamics(gsl_rng * random_variable){}

public:

    list <int> species_abundance;
    int printed_networks,network;
    void time_increment(gsl_rng * random_variable){
        dt=1./active_individuals_vector.size();
    };
    virtual void generate_network(gsl_rng * random_variable, double *mean_final_population=0,double *p= nullptr){
	//cout<<"Begin gen_netw 353"<<endl;
        bool overtime=false;
        double max_time=100000;
        int current_time=0;
        int max_size=4000;

	//cout<<"Begin initial"<<endl;
        initialize(random_variable);
	//cout<<"Finish init"<<endl;
        cout<<"CREATING NETWORK "<<network+1<<endl;
        t=0.0;
        current_mutations=0;
        time_increment(random_variable);
        TOTAL_ftime=time(NULL);
        //cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
        do{
            t+=dt;
            ///TIME UPDATE
            ftime=time(NULL);
            TOTAL_ftime=time(NULL);
            if(difftime(ftime, itime)>=(0.01)*60){
                //cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
                itime=time(NULL);
            }
	    cout<<"t = "<<t<<endl;
            one_step_dynamics(random_variable);
            if(before_individuals*tree_individual_struct_size>MAX_MBITS) eliminate_dead_branches();
        }while((ancestors_list.back().alive_offspring!=active_individuals_vector.size())&&(!overtime));

        eliminate_dead_branches();
        eliminate_dead_branches();

        //cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;

        (*mean_final_population)*=network;
        (*mean_final_population)+=active_individuals_vector.size();

        network++;

        (*mean_final_population)/=network;
    }
    void print_network(){

    }
    void assign_species(gsl_rng * random_variable){

        double mutation_probability;
        it_ancestor=ancestors_list.begin();
        offspring_begin_it=(*it_ancestor).offspring.begin();
        offspring_end_it=(*it_ancestor).offspring.end();

        #ifndef MUTATION_IN_DYNAMICS
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++) (*it_offspring1).initialize();

        ///////////////////////// MUTATIONS ///////////////////////
        ++offspring_begin_it;
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){

            mutation_probability=1-pow((1-mutation_rate),(*it_offspring1).intermediate_births);
            if(gsl_rng_uniform(random_variable)<mutation_probability){
                (*it_offspring1).mutation=true;
            }
        }
        #endif

        ///////////////////////// SPECIES ASSIGNMENT ///////////////////////
        for(int i=0;i<active_individuals_vector.size();i++){

            list <tree_individual>::iterator current_it;
            list<int>::iterator aux_it;
            vector <species_active_individual> to_assign_vector;
            species_active_individual to_assign;

            current_it=active_individuals_vector[i].offspring_it;
            while((!(*current_it).mutation)&&(!(*current_it).species_assignated)&&((*current_it).label!=0))
            {
                to_assign.individual_it=current_it;
                to_assign_vector.push_back(to_assign);
                current_it=(*current_it).parent;
            }

            if((*current_it).species_assignated){    ///IF IT IS ASSIGNED
                for(int j=0;j<to_assign_vector.size();j++){
                    (*to_assign_vector[j].individual_it).species_abundance=(*current_it).species_abundance;
                    (*to_assign_vector[j].individual_it).species_assignated=true;
                }
            }
            else{                                                           ///IF IT IS MUTATED OR IS THE ANCESTOR
                to_assign.individual_it=current_it;
                to_assign_vector.push_back(to_assign);
                species_abundance.push_back(0);
                aux_it=--species_abundance.end();

                for(int j=0;j<to_assign_vector.size();j++){
                    (*to_assign_vector[j].individual_it).species_abundance=aux_it;
                    (*to_assign_vector[j].individual_it).species_assignated=true;
                }
            }
            to_assign_vector.clear();
        }

        for(int i=0;i<active_individuals_vector.size();i++){
            (*(*active_individuals_vector[i].offspring_it).species_abundance)++;
        }
    };
    void get_abundances(int trial,double Lx_initial,double D_factor){

        ofstream SvsN_file;
        char SvsN_name[128];

        sprintf(SvsN_name,"SNL_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_rate,int(NNN_forward),Lx,Lx_initial,random_seed);
        SvsN_file.open(SvsN_name, ios::out | ios::trunc);

        list<int>::iterator it;
        for(it=species_abundance.begin();it!=species_abundance.end();it++) {
            SAD_histogram[(*it)]++;
            N_sample+=(*it);
            if((*it)>0) S_sample++;
        }
        TOTAL_TRIALS++;
        SvsN_file<<Lx_initial<<" "<<N_sample/double(TOTAL_TRIALS)<<" "<<S_sample/double(TOTAL_TRIALS)<<" "<<TOTAL_TRIALS<<endl;
    }

    void get_measures_subsamples(double Lxi,int trial,double interaction_range,bool print_spatial,ofstream *spatial_data){
	//cout<<"2.0"<<endl;
	cout<<"Inter rang in get_meas_subs = "<<interaction_range<<endl;
        double Lxi_local=Lxi;
        double Lx_local=Lxi;
        list<int>::iterator it;
        int number_of_subareas=1;//20;
        vector <subsample> subsample_vector;
        subsample subsample_new_element;
        int current_cell,current_cell2;
        int species_label, species_label2;
        int number_of_species;
        double current_distance,alternative_distance;
        double current_distance_x;
        double current_distance_y;
        double alternative_distance_x;
        double alternative_distance_y;

        vector <bool> checked_species_alpha;
        vector <int> current_alpha;

	//cout<<"2.1"<<endl;
        ///Here I construct the vector of subsamples
        for(int i=0;i<number_of_subareas;i++){
            vector <int> species_abundance_new_list;
            for(int j=0;j<species_abundance.size();j++) species_abundance_new_list.push_back(0);
            subsample_new_element.number_species=0;
            subsample_new_element.number_individuals=0;
            subsample_new_element.occupied=false;
            subsample_new_element.species_abundance_vector=species_abundance_new_list;
            subsample_vector.push_back(subsample_new_element);
        }

	    for(int k=0;k<MAX_NUMBER_OF_INDIVIDUALS;k++){
		    checked_species_alpha.push_back(false);
	    }
        for(int k=0;k<number_of_spatial_radius;k++){
            current_alpha.push_back(0);
        }

	//cout<<"2.2"<<endl;

        ///Here I count the species abundance at each subarea
        int max_cell=0;
        int min_cell=10000;
        for(int i=0;i<active_individuals_vector.size();i++){
            current_cell=0;
            species_label=distance(species_abundance.begin(), (*active_individuals_vector[i].offspring_it).species_abundance);
            subsample_vector[current_cell].species_abundance_vector[species_label]++;
            subsample_vector[current_cell].number_individuals++;
            subsample_vector[current_cell].occupied=true;
            if(current_cell>max_cell) max_cell=current_cell;
            if(current_cell<min_cell) min_cell=current_cell;
        }
        max_cell++;


        ///Here I count the number of species at each subarea
        number_of_occupied_subareas=0;
        for(int i=min_cell;i<max_cell;i++){
            number_of_species=0;
            for(int j=0;j<subsample_vector[i].species_abundance_vector.size();j++){
                if(subsample_vector[i].species_abundance_vector[j]!=0){
			        SAD_histogram_subsamples_mean[subsample_vector[i].species_abundance_vector[j]]++;
			        number_of_species++;
                }
            }
            subsample_vector[i].number_species=number_of_species;
            if(subsample_vector[i].occupied) number_of_occupied_subareas++;
        }

        ///Here I measure beta-diversity and heterozigozity
        int all_pairs=0;
        int equal_pairs=0;
        double mean_final_individuals_distance_aux=0;
        double mean_final_individuals_distance_all_aux=0;
        double min_final_individuals_distance_all_aux=10;
        double max_final_individuals_distance_all_aux=0;

        for(int i=0;i<active_individuals_vector.size();i++){
            species_label=distance(species_abundance.begin(), (*active_individuals_vector[i].offspring_it).species_abundance);
            if(print_spatial) (*spatial_data)<<(*active_individuals_vector[i].offspring_it).x_initial<<" "<<(*active_individuals_vector[i].offspring_it).y_initial<<" "<<species_label<<endl;
            for(int i2=i+1;i2<active_individuals_vector.size();i2++){
                species_label2=distance(species_abundance.begin(), (*active_individuals_vector[i2].offspring_it).species_abundance);
                current_distance_x=fabs((*active_individuals_vector[i].offspring_it).x_initial-(*active_individuals_vector[i2].offspring_it).x_initial);
                current_distance_y=fabs((*active_individuals_vector[i].offspring_it).y_initial-(*active_individuals_vector[i2].offspring_it).y_initial);

                if(current_distance_x>0.5*Lx) current_distance_x-=0.5*Lx;
                if(current_distance_y>0.5*Lx) current_distance_y-=0.5*Lx;

                current_distance=sqrt(pow(current_distance_x,2)+(pow(current_distance_y,2)));

                mean_final_individuals_distance_all_aux+=current_distance;
                if(current_distance>max_final_individuals_distance_all_aux) max_final_individuals_distance_all_aux=current_distance;
                if(current_distance<min_final_individuals_distance_all_aux) min_final_individuals_distance_all_aux=current_distance;
                if(species_label==species_label2){
                    mean_heterozigosity_diff_pairs[int(current_distance/interaction_range)]++;
                    mean_final_individuals_distance_aux+=current_distance;
                    equal_pairs++;
                }
                all_pairs++;
                mean_heterozigosity_all_pairs[int(current_distance/interaction_range)]++;
            }
            if(!checked_species_alpha[species_label]){

                //ORIGIN (0.5,0.5)
		            Lxi_local = Lx;//7.5;
                current_distance=sqrt(pow((*active_individuals_vector[i].offspring_it).x_initial-Lxi_local/2.,2)+pow((*active_individuals_vector[i].offspring_it).y_initial-Lxi_local/2.,2));
                current_alpha[int(current_distance/interaction_range)]++;

            }
	        checked_species_alpha[species_label]=true;
        }


        if(equal_pairs!=0) mean_final_individuals_distance+=mean_final_individuals_distance_aux/double(equal_pairs);
        if(all_pairs!=0){
            mean_final_individuals_distance_all+=mean_final_individuals_distance_all_aux/double(all_pairs);
            min_final_individuals_distance_all+=min_final_individuals_distance_all_aux;
            max_final_individuals_distance_all+=max_final_individuals_distance_all_aux;
        }

        int cummulated_alpha=0;
	    for(int i=0;i<number_of_spatial_radius;i++){
        	cummulated_alpha+=current_alpha[i];
		    mean_alpha[i]+=cummulated_alpha;
        }

        ///Here I average over all subsamples
        double number_pair_subareas=0;
        for(int i=min_cell;i<max_cell;i++){
            mean_number_individuals_subarea+=subsample_vector[i].number_individuals/double(number_of_occupied_subareas);
            mean_number_species_subarea+=subsample_vector[i].number_species/double(number_of_occupied_subareas);
        }
        species_abundance.clear();


    }

    void print_abundances(double D_factor,double Lxi){

        ///////////////////////// SAD ///////////////////////
        sprintf(SAD_subsample_name,"SAD_subsample_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_rate,NNN_forward,Lx,Lxi,random_seed);
        SAD_subsample.open(SAD_subsample_name, ios::out | ios::trunc);
        for(int i=0;i<MAX_NUMBER_OF_INDIVIDUALS;i++){
            if(SAD_histogram_subsamples_mean[i]!=0){
                SAD_subsample<<i<<" "<<SAD_histogram_subsamples_mean[i]/double(TOTAL_TRIALS*number_of_occupied_subareas)<<" "<<network<<" "<<TOTAL_TRIALS<<endl;
            }
        }
        SAD_subsample.close();
    }
    void clear(){
        active_individuals_vector.clear();
	      ancestors_list.clear();
    }
};


///COALESCESNCE TREE CLASS
class coalescing_random_walk_tree:public tree{

public:
    vector <double> n_coalescence_vector;
    unsigned int border_colisions;
private:
    tree_individual new_individual;

    list<individuals_in_area> area_list;
    list<list<individuals_in_area>::iterator> area_list_occupied;

    list<individuals_in_area>::iterator **area;
    //int **area;
    int number_of_areas_x, number_of_areas_y;
    double space_between_areas,delta,sqrtdt;
    int number_of_periodic_cells_x;
    int number_of_periodic_cells_y;
    double *jet_parameters;
    double mutation_probability;
    double diffusion_rate;
    double dt_continuous;
    int particles_limit;
    double Lx_ini,Ly_ini,x_ini,y_ini;
    double ndnd,coalescence_probability,real_coalescence_probability;

    void initialize_area_list(){

        list<individuals_in_area>::iterator aux_it;
        area_list.clear();

        for(int ix=0;ix<number_of_areas_x;ix++){
            for(int iy=0;iy<number_of_areas_y;iy++){

                individuals_in_area aux_area_list_element;
                area_list.push_back(aux_area_list_element);
                aux_it=area_list.end();
                --aux_it;
                area[ix][iy]=aux_it;
            }
        }
    }
    virtual void initialize(gsl_rng * random_variable){

        ancestor_individual new_ancestor;
        list<ancestor_individual>::iterator ancestor;
	      initialize_area_list();
        ancestors_list.push_back(new_ancestor);

        ///The first position of offspring vector (parent_dir=0) is reserved to the individual 0="parent of ancestors"
        new_individual.label=0;
        ancestors_list.back().offspring.push_back(new_individual);

        last_label=0;
        area_list_occupied.clear();
        final_spatial_data.open(final_spatial_name, ios::in);
        for(unsigned int i=0;i<NNN;i++){
            last_label++;
            new_individual.label=last_label;
            new_individual.t_birth=0;


            new_active_individual.Lx=Lx;
            new_active_individual.Ly=Ly;
            new_active_individual.number_cells_x=number_of_periodic_cells_x;
            new_active_individual.number_cells_y=number_of_periodic_cells_y;
            new_active_individual.initial_set(random_variable,&final_spatial_data,0,jet_parameters,Lx_ini,Ly_ini,x_ini,y_ini);

            new_individual.x_initial=new_active_individual.x;
            new_individual.y_initial=new_active_individual.y;

            ancestors_list.back().offspring.push_back(new_individual);
            new_active_individual.offspring_it=--ancestors_list.back().offspring.end();
	          new_active_individual.current_area=area[0][0];//[int(new_active_individual.x/space_between_areas)][int(new_active_individual.y/space_between_areas)];
	          active_individuals_vector.push_back(new_active_individual);
            (*new_active_individual.current_area).individual_list.push_back(active_individuals_vector.size()-1);
	    area_list_occupied.push_back(new_active_individual.current_area);
        }
        final_spatial_data.close();
        offspring_begin_it=ancestors_list.back().offspring.begin();
        offspring_end_it=ancestors_list.back().offspring.end();
    }
    virtual void particular_change_variables(){

        int current_label,new_individual_label;
        unsigned int i=0;

        active_individuals_vector.clear();
        offspring_begin_it=++ancestors_list.back().offspring.begin();
        offspring_end_it=offspring_begin_it;
        advance(offspring_end_it,NNN);

        for (it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            new_active_individual.offspring_it=it_offspring1;
            active_individuals_vector.push_back(new_active_individual);
            i++;
        }

        offspring_begin_it=++ancestors_list.back().offspring.begin();
        offspring_end_it=ancestors_list.back().offspring.end();
    }
    void contour_conditions(int i){

        if(active_individuals_vector[i].x>number_of_periodic_cells_x*Lx) border_colisions++;
        if(active_individuals_vector[i].y>number_of_periodic_cells_y*Ly) border_colisions++;

        if(active_individuals_vector[i].x<0) border_colisions++;
        if(active_individuals_vector[i].y<0) border_colisions++;

        active_individuals_vector[i].x=fmod((active_individuals_vector[i].x+Lx*number_of_periodic_cells_x),Lx*number_of_periodic_cells_x);
        active_individuals_vector[i].y=fmod((active_individuals_vector[i].y+Ly*number_of_periodic_cells_y),Ly*number_of_periodic_cells_y);
    };
    int check_occupied(int isite,list<int>::iterator *coalescent_site=nullptr){

        auto jaux=active_individuals_vector[isite].offspring_it;
        double ix=active_individuals_vector[isite].x;
        double iy=active_individuals_vector[isite].y;
        double jx,jy;
        int number_of_individuals_area=(*active_individuals_vector[isite].current_area).individual_list.size();
        auto jsite_it=(*active_individuals_vector[isite].current_area).individual_list.begin();
        int occupied=0;
        int kxi,kyi,kxj,kyj;

        kxi=(int)floor(ix*delta);
        kyi=(int)floor(iy*delta);
        for(int i=0;i<number_of_individuals_area;i++){
            jx=active_individuals_vector[(*jsite_it)].x;
            jy=active_individuals_vector[(*jsite_it)].y;
            kxj=(int)floor(jx*delta);
            kyj=(int)floor(jy*delta);

            if((kxi==kxj)&&(kyi==kyj)&&(ix!=jx)&&(iy!=jy)){
                (*coalescent_site)=jsite_it;
                occupied++;
            }
            ++jsite_it;
        }
        return occupied;
    }

public:
    coalescing_random_walk_tree(int N,double mutprob,int r,double h,double nc,double DIFF_in,double delta_in,
                                int number_of_spatial_radius_in=0,double Lx_in=-1,double Ly_in=-1,double *extra_parameters=nullptr,
                                double Lx_i=0, double Ly_i=0, double xi=0,double yi=0):tree(N,mutprob,r,number_of_spatial_radius_in,Lx_in,Ly_in){
        Lx=Lx_in;
        Ly=Ly_in;
        dt=1./N;
        dt_continuous=0.001;
        model=MODEL_NAME;
        Lx_ini=Lx_i;
        Ly_ini=Ly_i;
        x_ini=xi;
        y_ini=yi;
        border_colisions=0;
        number_of_periodic_cells_x=nc;
        number_of_periodic_cells_y=nc;
        particles_limit=200;
	      delta=delta_in;
        space_between_areas=Lx;//7.5;   // 1/delta;   //1;
        number_of_areas_x=(number_of_periodic_cells_x*Lx_in/space_between_areas)+1;//(int)(Lx*delta)+1;
        number_of_areas_y=(number_of_periodic_cells_y*Ly_in/space_between_areas)+1;//(int)(Ly*delta)+1;

        diffusion_rate=DIFF_in;
        sqrtdt=sqrt(2.0*dt_continuous*diffusion_rate);
        area=(list<individuals_in_area>::iterator **) malloc(sizeof(list<individuals_in_area>::iterator *)*(number_of_areas_x)); //(int **) malloc(sizeof(int *)*number_of_areas_x);
        for(int i=0;i<(number_of_areas_x);i++) *(area+i)=(list<individuals_in_area>::iterator *)malloc(sizeof(list<individuals_in_area>::iterator)*(number_of_areas_y));  //{*(area+1)=new list<individuals_in_area>::iterator[number_of_areas_y];} // *(area+i)=(list<individuals_in_area>::iterator *)malloc(sizeof(list<individuals_in_area>::iterator)*(number_of_areas_y)); //(int i=0;i<number_of_areas_x;i++) *(area+i)=(int *)malloc(sizeof(int)*number_of_areas_y);
        jet_parameters=(double *) malloc(sizeof(double *)*8);
        jet_parameters[0]=extra_parameters[0];
        jet_parameters[1]=extra_parameters[1];
        jet_parameters[2]=extra_parameters[2];
        jet_parameters[3]=extra_parameters[3];
        jet_parameters[4]=extra_parameters[4];
        jet_parameters[5]=extra_parameters[5];
        jet_parameters[6] = extra_parameters[6];//(double)(Lx*Lx/16384.)/(0.05*0.05/2);  //(Lx*Lx/NNN)/(0.05*0.05/2);
        NNN_forward=N;//extra_parameters[6];
        NNN=N;
        mutation_rate=mutprob;
        mutation_probability=dt_continuous*mutation_rate;
        char aux_name[5000];
        ndnd=(double)NNN/(delta*delta*Lx_i*Lx_i);
        coalescence_probability=dt_continuous*1;
        N_sample=0;
        S_sample=0;
    }
    void input_ncoal(double factor=1){
        ifstream extra_data;
        char extra_name[5000];
        sprintf(extra_name,"ncoalescence_N%d_seed%d.dat",NNN,random_seed);
        extra_data.open(extra_name, ios::in);
        double distance,mean_n,min_n,max_n;
        double prev_distance=0;

        while(extra_data>>distance>>mean_n>>min_n>>max_n){
            n_coalescence_vector.push_back(mean_n/factor);
        }
        extra_data.close();
    }
    virtual void generate_network(gsl_rng * random_variable, ofstream& coal_part_data,  double *mean_final_population=0,double *p= nullptr){
	//cout<<"Begin gen_netw 822"<<endl;
#ifdef OUTPUT_TREE
        ofstream actives_t_data;
        char actives_t_name[5000];
        sprintf(actives_t_name,"ACTIVES-t_generation_tree_%s_N%d_Network%d.dat",model.c_str(),NNN,random_seed);
        actives_t_data.open(actives_t_name, ios::out | ios::trunc);
#endif
        itime=time(NULL);
        int iselected;
        time_t itime=time(NULL),ftime;
        int next_particles=NNN;
        initialize(random_variable);


        cout<<"CREATING NETWORK "<<network+1<<" IN "<<model.c_str()<<endl;

        t=dt_continuous;
	cout<<"t = "<<t<<endl;
        current_mutations = 0;
	current_coalescents = 0;
	curr_no_coal_no_mut = 0;
	sum_coal_prob=0;
        double tnext=t+dt_continuous;
        double last_t=t;
        bool flag_lookup=true;
	int number_of_coalescent_particles;
	double time_bef_mut, time_aft_mut, time_bef_coal, time_aft_coal,time_bef_move, time_aft_move;
	double  total_mut_time=0, total_coal_time = 0, total_move_time= 0;

        do{
            vector<int> aux_actives(active_individuals_vector.size());
            iota (std::begin(aux_actives), std::end(aux_actives), 0);

	    int num_coal_part_tot_step = 0;
	    int occupied;
            while(aux_actives.size()>0){
                iselected=aux_actives[gsl_rng_uniform_int(random_variable,aux_actives.size())];
                #ifdef MUTATION_IN_DYNAMICS
	        //list<int>::iterator coalescent_tree_site;
                if(gsl_rng_uniform(random_variable)<1-mutation_probability){
		    //time_bef_coal = time(NULL);
                    coalesce(random_variable,iselected,&aux_actives,coal_part_data);
                }
                else{
		    //time_bef_mut = time(NULL);
                    mutate(iselected,&aux_actives);
		    //time_aft_mut = time(NULL);
		    current_mutations++;

                }
                #else
                coalesce(random_variable,iselected,&aux_actives,&actives_t_data,coal_part_data);
                #endif
            }

            for(int i=0;i<active_individuals_vector.size();i++) move(random_variable,i);
	    t+=dt_continuous;

            ftime=time(NULL);
            TOTAL_ftime=time(NULL);
            if(difftime(ftime, itime)>=(2)*60){
                cout<<"!!t="<<t<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
                //cout<<"!!t="<<t<<", curr mut="<<current_mutations<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
		itime=time(NULL);
		//cout<<"Total coal time = "<<total_coal_time<<"; Total mut time = "<<total_mut_time<<"; Total move time = "<<total_move_time<<endl;
		}

            if(flag_lookup){
                if(active_individuals_vector.size()<particles_limit){
                    //space_between_areas=number_of_periodic_cells_x*Lx;
                    update_areas();
                    flag_lookup=false;
                }
                else update_areas();
            }

        }while(active_individuals_vector.size()>1);

        last_step_dynamics(random_variable);
        particular_change_variables();
        network++;
#ifdef OUTPUT_TREE
        actives_t_data.close();
#endif
    }
    void update_areas(){
        for(auto area_it=area_list_occupied.begin();area_it!=area_list_occupied.end();++area_it){
            (*(*area_it)).individual_list.clear();
        }
        area_list_occupied.clear();
        for(int i=0;i<active_individuals_vector.size();i++){
            active_individuals_vector[i].current_area=area[int(active_individuals_vector[i].x/space_between_areas)][int(active_individuals_vector[i].y/space_between_areas)];
            (*active_individuals_vector[i].current_area).individual_list.push_back(i);
            area_list_occupied.push_back(active_individuals_vector[i].current_area);
        }
    }
    void move(gsl_rng * random_variable, int imove){

        double incr_x,incr_y;
	double K0x, K0y, K1x, K1y, K2x, K2y, K3x, K3y;
	#ifdef ADVECTION
        active_individuals_vector[imove].advectx(force,t,Ly,jet_parameters,0,0,0);
        active_individuals_vector[imove].advecty(force,t,Lx,jet_parameters,0,0,0);
        K0x = active_individuals_vector[imove].fieldx;
        K0y = active_individuals_vector[imove].fieldy;
        active_individuals_vector[imove].advectx(force,t,Ly,jet_parameters,K0x,K0y,dt_continuous/2);
        active_individuals_vector[imove].advecty(force,t,Lx,jet_parameters,K0x,K0y,dt_continuous/2);
        K1x = active_individuals_vector[imove].fieldx;
        K1y = active_individuals_vector[imove].fieldy;

        active_individuals_vector[imove].advectx(force,t,Ly,jet_parameters,K1x,K1y,dt_continuous/2);
        active_individuals_vector[imove].advecty(force,t,Lx,jet_parameters,K1x,K1y,dt_continuous/2);
        K2x = active_individuals_vector[imove].fieldx;
        K2y = active_individuals_vector[imove].fieldy;

        active_individuals_vector[imove].advectx(force,t,Ly,jet_parameters,K2x,K2y,dt_continuous);
        active_individuals_vector[imove].advecty(force,t,Lx,jet_parameters,K2x,K2y,dt_continuous);
        K3x = active_individuals_vector[imove].fieldx;
        K3y = active_individuals_vector[imove].fieldy;

	incr_x =sqrtdt*gsl_ran_gaussian(random_variable,1)+ dt_continuous*(K0x+2*K1x+2*K2x+K3x)/6;
	incr_y =sqrtdt*gsl_ran_gaussian(random_variable,1)+ dt_continuous*(K0y+2*K1y+2*K2y+K3y)/6;

	#else
	incr_x = sqrtdt*gsl_ran_gaussian(random_variable,1);
	incr_y = sqrtdt*gsl_ran_gaussian(random_variable,1);
	#endif

        active_individuals_vector[imove].x+=incr_x;
        active_individuals_vector[imove].y+=incr_y;

        contour_conditions(imove);

        active_individuals_vector[imove].fieldoldx=active_individuals_vector[imove].fieldx;
        active_individuals_vector[imove].fieldoldy=active_individuals_vector[imove].fieldy;
    }
    void mutate(int imove, vector<int> *aux_actives){

        if(t<tdmax){
            mutations_tdistribution[int(t)]++;
            events_tdistribution[int(t)]++;
        }
        mutations_total++;
        current_mutations++;


        list<tree_individual>::iterator it_move=active_individuals_vector[imove].offspring_it;

        (*it_move).mutation=true;

        ///The father is created: a random individual of the ancestor family as it will not affect the results
        if(imove>0) (*it_move).parent=active_individuals_vector[imove-1].offspring_it;
        else (*it_move).parent=active_individuals_vector[imove+1].offspring_it;

        ///Update of active list: remove the one moving and change the (previously) active (son) iterator to the parent
        (*active_individuals_vector[imove].current_area).individual_list.remove(imove);
        (*active_individuals_vector.back().current_area).individual_list.remove(active_individuals_vector.size()-1);
        if(imove!=(active_individuals_vector.size()-1)){
            (*active_individuals_vector.back().current_area).individual_list.push_back(imove);
            active_individuals_vector[imove].current_area=active_individuals_vector.back().current_area;
        }
        active_individuals_vector[imove]=active_individuals_vector.back();
        active_individuals_vector.pop_back();

        (*aux_actives).erase(std::remove((*aux_actives).begin(), (*aux_actives).end(), active_individuals_vector.size()), (*aux_actives).end());

    }
    void coalesce(gsl_rng * random_variable, int imove, vector<int> *aux_actives, ofstream& coal_part_data,ofstream *extra_data= nullptr){

        list<int>::iterator coalescent_tree_site;
        list<tree_individual>::iterator coalescent_active_site_it,it_move;
        vector<flux_individual>::iterator active_it;
        int number_of_coalescent_particles;
        it_move=active_individuals_vector[imove].offspring_it;
        number_of_coalescent_particles=check_occupied(imove,&coalescent_tree_site);

	if(t == 0.001){
	coal_part_data<<number_of_coalescent_particles<<", ";}
        real_coalescence_probability=coalescence_probability*number_of_coalescent_particles;

        if(gsl_rng_uniform(random_variable)<real_coalescence_probability){

            if(t<tdmax){
                coalescences_tdistribution[int(t)]++;
                events_tdistribution[int(t)]++;
            }
            coalescences_total++;
	    current_coalescents++;

            ///COALESCENCE
            coalescent_active_site_it=active_individuals_vector[(*coalescent_tree_site)].offspring_it;
            int actual_coalescent=(*coalescent_tree_site);

            ///The father is created
            last_label++;
            new_individual.label=last_label;
            new_individual.t_birth=t;
            new_individual.x_initial=(*coalescent_active_site_it).x_initial;
            new_individual.y_initial=(*coalescent_active_site_it).y_initial;
            ancestors_list.back().offspring.push_back(new_individual);

#ifdef OUTPUT_TREE
            (*extra_data)<<t<<" "<<(*it_move).label<<" "<<(*coalescent_active_site_it).label<<" "<<last_label<<" "<<active_individuals_vector.size()-1<<endl;
#endif
            ///The parent of the two coalescent particles is the one that have just been created
            (*it_move).parent=--ancestors_list.back().offspring.end();
            (*coalescent_active_site_it).parent=--ancestors_list.back().offspring.end();

            ///Update of properties
            (*it_move).intermediate_births=(t-(*it_move).t_birth);

            ///Update of active list: remove the one moving and change the (previously) active (son) iterator to the parent
            (*active_individuals_vector[imove].current_area).individual_list.remove(imove);
            (*active_individuals_vector.back().current_area).individual_list.remove(active_individuals_vector.size()-1);
            if(imove!=(active_individuals_vector.size()-1)){
                (*active_individuals_vector.back().current_area).individual_list.push_back(imove);
                active_individuals_vector[imove].current_area=active_individuals_vector.back().current_area;
            }
            active_individuals_vector[actual_coalescent].offspring_it=--ancestors_list.back().offspring.end();
            active_individuals_vector[imove]=active_individuals_vector.back();
            active_individuals_vector.pop_back();

            (*aux_actives).erase(std::remove((*aux_actives).begin(), (*aux_actives).end(), active_individuals_vector.size()), (*aux_actives).end());
            (*aux_actives).erase(std::remove((*aux_actives).begin(), (*aux_actives).end(), actual_coalescent), (*aux_actives).end());
        }
        else{
            (*aux_actives).erase(std::remove((*aux_actives).begin(), (*aux_actives).end(), imove), (*aux_actives).end());
	    curr_no_coal_no_mut++;

        }
	//return number_of_coalescent_particles;
    }
    virtual void one_step_dynamics(gsl_rng * random_variable){

        list<int>::iterator coalescent_tree_site;
        list<tree_individual>::iterator coalescent_active_site_it,it_move;
        vector<flux_individual>::iterator active_it;
        int imove;
        bool occupied;

        #ifdef MUTATION_IN_DYNAMICS
        imove=gsl_rng_uniform_int(random_variable,active_individuals_vector.size());
        it_move=active_individuals_vector[imove].offspring_it;
        occupied=check_occupied(imove,&coalescent_tree_site);

        if(occupied){
            ///COALESCENCE

            coalescent_active_site_it=active_individuals_vector[(*coalescent_tree_site)].offspring_it;
            int actual_coalescent=(*coalescent_tree_site);

            ///The father is created
            last_label++;
            new_individual.label=last_label;
            ancestors_list.back().offspring.push_back(new_individual);

            ///The parent of the two coalescent particles is the one that have just been created
            (*it_move).parent=--ancestors_list.back().offspring.end();
            (*coalescent_active_site_it).parent=--ancestors_list.back().offspring.end();

            ///Update of active list: remove the one moving and change the (previously) active (son) iterator to the parent
            (*active_individuals_vector[imove].current_area).individual_list.remove(imove);
            (*active_individuals_vector.back().current_area).individual_list.remove(active_individuals_vector.size()-1);
            if(imove!=(active_individuals_vector.size()-1)){
                (*active_individuals_vector.back().current_area).individual_list.push_back(imove);
                active_individuals_vector[imove].current_area=active_individuals_vector.back().current_area;
            }
            active_individuals_vector[actual_coalescent].offspring_it=--ancestors_list.back().offspring.end();
            active_individuals_vector[imove]=active_individuals_vector.back();
            active_individuals_vector.pop_back();


        }

        #else
        imove=gsl_rng_uniform_int(random_variable,active_individuals_vector.size());
        it_move=active_individuals_vector[imove].offspring_it;
        occupied=check_occupied(imove,&coalescent_tree_site);

        if(occupied){
            ///COALESCENCE

            coalescent_active_site_it=active_individuals_vector[(*coalescent_tree_site)].offspring_it;
            int actual_coalescent=(*coalescent_tree_site);

            ///The father is created
            last_label++;
            new_individual.label=last_label;
            new_individual.t_birth=t;
            ancestors_list.back().offspring.push_back(new_individual);

            ///The parent of the two coalescent particles is the one that have just been created
            (*it_move).parent=--ancestors_list.back().offspring.end();
            (*coalescent_active_site_it).parent=--ancestors_list.back().offspring.end();

            ///Update of properties
            (*it_move).intermediate_births=(t-(*it_move).t_birth);

            ///Update of active list: remove the one moving and change the (previously) active (son) iterator to the parent
            (*active_individuals_vector[imove].current_area).individual_list.remove(imove);
            (*active_individuals_vector.back().current_area).individual_list.remove(active_individuals_vector.size()-1);
            if(imove!=(active_individuals_vector.size()-1)){
                (*active_individuals_vector.back().current_area).individual_list.push_back(imove);
                active_individuals_vector[imove].current_area=active_individuals_vector.back().current_area;
            }
            active_individuals_vector[actual_coalescent].offspring_it=--ancestors_list.back().offspring.end();
            active_individuals_vector[imove]=active_individuals_vector.back();
            active_individuals_vector.pop_back();
        }
        #endif
        ///TIME STEP UPDATE
        time_increment(random_variable);
    }
    void last_step_dynamics(gsl_rng * random_variable){

        ///LAST COALESCENCE
        (*active_individuals_vector[0].offspring_it).parent=offspring_begin_it;
    }
};
