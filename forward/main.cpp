#include "classes.h"

using namespace std;

int main(int argc, char *argv[])
{
    ///////////////////////// PARAMETERS ///////////////////////
    const double pi = 4.0 * atan(1.0);
    int number_of_spatial_radius=0;		//Fluid force
    double NNN=16384;				//Initial number of alive_individuals
    int random_seed=1;
    double mu=1.0;
    double refresh_time=30;			//REAL TIME BETWEEN SAD OUTPUTS (in minutes)
    double mutation_probability=0.001;	//Mutation probability
    double Lx=7.5;				//Linear size x axis
    int NETWORKS_TRIALS=100;
    int NETWORKS_TO_PRINT=1;
    int ABUNDANCES_TRIALS=100;
    int READ_NETWORK=0;			//Set to 0 if not wanted to read the network
    double B0_jet=1.2;
    double epsilon_jet=0.3;
    double w_jet=0.4;
    double L_jet=7.5;
    double c_jet=0.12;
    double k_jet = 2*pi/L_jet;
    double Fi_jet = pi/2;

    double sigma = 0.05;

    double Lxi=L_jet;
    double D_factor=1;
    int NNN_sample=16000;
    double xcell=0.5;
    double ycell=0.5;

    //double resc_const = (D_factor/(NNN))/(3*pow(10, -9));          //sqrt(3)*pow(10, -3);

    //When running the program
    for (int i = 1; i < argc; i++) {
        if(string(argv[i]) == "-N"){
            NNN = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Lx"){
            Lx = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Lxi"){
            Lxi = atof(argv[i + 1]);
        }else if(string(argv[i]) == "-mutation_prob"){
            mutation_probability = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-random_seed"){
            random_seed = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-networks_to_average"){
            NETWORKS_TRIALS = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-networks_to_print"){
            NETWORKS_TO_PRINT = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-read_network"){
            READ_NETWORK = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-trials_per_network"){
            ABUNDANCES_TRIALS = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-rt") {
            refresh_time = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-B0_jet"){
            B0_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-epsilon_jet"){
            epsilon_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-w_jet"){
            w_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-L_jet"){
            L_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-c_jet"){
            c_jet = atof(argv[i + 1]);
        }  else if(string(argv[i]) == "-D_factor"){
            D_factor = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Ns"){
            NNN_sample = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-xcell"){
            xcell = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-ycell"){
            ycell = atof(argv[i + 1]);
        } else if(i%2!=0){
            cout<<endl<<" CAUTION!: "<<string(argv[i])<<" HAVE NOT SET ANY PARAMETER!"<<endl;
        }
    }

    ///////////////////////// VARIABLES ///////////////////////
    NNN_sample=NNN_sample;
    double Ly=Lx;//1;
    time_t itime=time(NULL),ftime,TOTAL_itime=time(NULL),TOTAL_ftime;
    int trial;
    double DIFF= (D_factor*Lx*Lx)/(NNN);                     //Diffusion
    double interaction_range;
    double jet_parameters[8];
    //double g_jet=0.001;//1;
    //double l_jet=2*M_PI*L_jet;
    //double a_jet=1;
    //double v0_jet=1;
    //double k_jet=1;
    jet_parameters[0]=B0_jet;
    jet_parameters[1]=epsilon_jet;
    jet_parameters[2]=Fi_jet;
    jet_parameters[3]=w_jet;
    jet_parameters[4]=c_jet;
    jet_parameters[5]=L_jet;
    jet_parameters[6]=k_jet;
    jet_parameters[7]=(double)(Lx*Lx/16384.)/(sigma*sigma/2);  //rescaling const
    double oneoverd = (int)floor(sqrt(NNN))/double(Lx);
    interaction_range=1./oneoverd; //1./(int)floor(sqrt(NNN))/double(Lx);;//   double(Lx)/(int)floor(sqrt(NNN));//(int)floor(sqrt(NNN))/double(Lx); (int)floor(sqrt(NNN))/double(Lx);
    number_of_spatial_radius=int((Lx+Lx/2)/interaction_range);//int(1.5/interaction_range);
    fluid2D_tree current_tree(NNN,mutation_probability,random_seed,number_of_spatial_radius,Lx,Ly,mu,DIFF,jet_parameters);


    ///////////////////////// INITIALIZATION ///////////////////////
    //INITIALIZE_RANDOM_SEED(rand_seed):
    int rand_seed=random_seed;
    const gsl_rng_type * gsl_rng_T;                 //Random variable
    gsl_rng * random_variable;                		//Random variable
    gsl_rng_env_setup();
    gsl_rng_default_seed=rand_seed;
    gsl_rng_T=gsl_rng_default;
    random_variable=gsl_rng_alloc(gsl_rng_T);
    if(refresh_time>1) refresh_time=(refresh_time+gsl_rng_uniform_int(random_variable,5))*60;

    ofstream spatial_data;
    char spatial_name[300];
    sprintf(spatial_name,"Spatial_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
    bool print_spatial=true;
    spatial_data.open(spatial_name, ios::out | ios::trunc);

    char final_spatial_name[300];
    ofstream final_spatial;

    ofstream final_data;
    char final_name[300];
    sprintf(final_name,"final_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
    final_data.open(final_name, ios::out | ios::trunc);

    ofstream current_events_data;
    char current_events_name[300];
    sprintf(current_events_name,"current_events_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
    current_events_data.open(current_events_name, ios::out | ios::trunc);


    ///////////////////////// NETWORKS LOOP ///////////////////////
    double mean_final_population=0;
    current_tree.network=0;
    double N_mean=0,N2_mean=0;
    double L_mean=0,L2_mean=0;
    while(current_tree.network<NETWORKS_TRIALS){//while(1){//
        current_tree.generate_network(random_variable,&mean_final_population,jet_parameters);
        current_tree.get_actives_index_vector();

        if(ABUNDANCES_TRIALS!=0){
            ///////////////////////// OBTAINING SAD AND SRC ///////////////////////
            trial=0;
            cout<<"OBTAINING SAD AND SRC"<<endl;
            while(trial<ABUNDANCES_TRIALS){
                current_tree.get_sample_index_vector(NNN_sample,random_variable);
                current_tree.assign_species(random_variable);
                current_tree.get_abundances(trial);

                if(current_tree.actives_index_vector.size()>NNN_sample) current_tree.get_measures_nosubsamples(Lxi,trial,interaction_range,print_spatial,&spatial_data,NNN_sample,random_variable);
                print_spatial=false;
                trial++;
                ftime=time(NULL);
                TOTAL_ftime=time(NULL);
                if(difftime(ftime, itime)>refresh_time){
                    current_tree.print_abundances(Lxi,D_factor);
                    cout<<"TRIALS="<<trial<<", Running time="<<difftime(TOTAL_ftime, TOTAL_itime)/60.<<" min"<<endl;
                    itime=time(NULL);
                }
		//cout<<"0"<<endl;
                current_events_data.close();
		//cout<<"1"<<endl;
                current_events_data.open(current_events_name, ios::out | ios::app);
		//cout<<"2"<<endl;
                current_events_data<<current_tree.current_mutations<<" "<<double(trial*current_tree.network)<<endl;
                //cout<<"3"<<endl;
		current_tree.current_mutations=0;
		//cout<<"4"<<endl;
                current_events_data.close();
		//cout<<"End of the abunt. loop"<<endl;
            }
            current_tree.print_abundances(Lxi,D_factor);


        }


        ofstream extra_data;
        char extra_name[300];
        sprintf(extra_name,"SNL_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        extra_data.open(extra_name, ios::out | ios::trunc);
        extra_data<<Lxi<<" ";
        extra_data<<NNN/double(current_tree.number_of_subareas)<<" ";
        extra_data<<current_tree.mean_number_individuals_subarea/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        extra_data<<current_tree.mean_number_species/double(current_tree.number_of_subareas*ABUNDANCES_TRIALS*current_tree.network)<<" ";
        extra_data<<current_tree.mean_number_species_subarea/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        extra_data<<current_tree.number_of_subareas<<" "<<current_tree.number_of_occupied_subareas<<" ";
        extra_data<<ABUNDANCES_TRIALS<<" "<<current_tree.network<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        extra_data.close();


        ofstream heterozigosity_data;
        char heterozigosity_name[300];
        sprintf(heterozigosity_name,"Beta-diversity_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        heterozigosity_data.open(heterozigosity_name, ios::out | ios::trunc);

        for(int i=0;i<current_tree.mean_heterozigosity_all_pairs.size();i++){
            if(current_tree.mean_heterozigosity_all_pairs[i]!=0){
                heterozigosity_data<<i*interaction_range<<" ";
                heterozigosity_data<<current_tree.mean_heterozigosity_diff_pairs[i]/double(current_tree.mean_heterozigosity_all_pairs[i])<<" ";
                heterozigosity_data<<ABUNDANCES_TRIALS<<" "<<current_tree.network<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
            }
        }

        ofstream alpha_data;
        char alpha_name[300];
        sprintf(alpha_name,"Alpha-diversity_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        alpha_data.open(alpha_name, ios::out | ios::trunc);
        for(int i=0;i<current_tree.mean_alpha.size();i++){
            if(current_tree.mean_alpha[i]!=0){
                alpha_data<<i*interaction_range<<" ";
                alpha_data<<current_tree.mean_alpha[i]/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
                alpha_data<<ABUNDANCES_TRIALS<<" "<<current_tree.network<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
            }
        }
        alpha_data.close();

        ofstream tdistribution_data;
        char tdistribution_name[300];
        sprintf(tdistribution_name,"tdistribution_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        tdistribution_data.open(tdistribution_name, ios::out | ios::trunc);

        for(double ipoint=0;ipoint<current_tree.tdistribution_points;ipoint++){
            if((current_tree.coalescences_tdistribution[ipoint]>1)||(current_tree.mutations_tdistribution[ipoint]>1)){
                tdistribution_data<<ipoint<<" "<<current_tree.coalescences_tdistribution[ipoint]/double(current_tree.events_tdistribution[ipoint]);
                tdistribution_data<<" "<<current_tree.mutations_tdistribution[ipoint]/double(current_tree.events_tdistribution[ipoint]);
                tdistribution_data<<" "<<current_tree.events_tdistribution[ipoint]<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
            }
            else tdistribution_data<<ipoint<<" . . "<<current_tree.events_tdistribution[ipoint]<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        }
        tdistribution_data.close();

        ofstream total_events_data;
        char total_events_name[300];
        sprintf(total_events_name,"total_events_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        total_events_data.open(total_events_name, ios::out | ios::trunc);
        total_events_data<<current_tree.coalescences_total/double(ABUNDANCES_TRIALS*current_tree.network);
        total_events_data<<" "<<current_tree.mutations_total/double(ABUNDANCES_TRIALS*current_tree.network)<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        total_events_data.close();



        ofstream final_separation_data;
        char final_separation_name[300];
        sprintf(final_separation_name,"final_separation_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN),Lx,Lxi,random_seed);
        final_separation_data.open(final_separation_name, ios::out | ios::trunc);
        final_separation_data<<current_tree.mean_final_individuals_distance/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.mean_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.min_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.max_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        final_separation_data.close();


        N_mean+=current_tree.get_Nactives();
        L_mean+=current_tree.max_final_individuals_distance_all_network/sqrt(2.);

        final_data<<current_tree.get_Nactives()<<" "<<N_mean/double(current_tree.network)<<" ";
        final_data<<current_tree.max_final_individuals_distance_all_network/sqrt(2.)<<" "<<L_mean/double(current_tree.network)<<" ";
        final_data<<current_tree.network<<endl;


        current_tree.clear();

        TOTAL_ftime=time(NULL);
        cout<<"NETWORK "<<current_tree.network<<" ENDED. TOTAL RUNNING TIME="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl<<endl;


    }

    final_data.close();
    return 0;
}
