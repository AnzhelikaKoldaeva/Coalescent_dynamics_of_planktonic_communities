// 2D particle simulations, with and without flow

#include "classes.h"

using namespace std;

int main(int argc, char *argv[])
{
    ///////////////////////// PARAMETERS ///////////////////////
    const double pi = 4.0 * atan(1.0);

    ///A priori
    int number_of_spatial_radius=0;
    int random_seed=1;
    double refresh_time=30;             //REAL TIME BETWEEN SAD OUTPUTS (in minutes)
    double mutation_probability=0.001;  //Mutation probability
    double Lx=7.5;                        //Linear size x axis

    int NETWORKS_TRIALS=400;
    int NETWORKS_TO_PRINT=1;
    int ABUNDANCES_TRIALS=100;

    double B0_jet=1.2;
    double epsilon_jet=0.3;
    double w_jet=0.4;
    double L_jet=7.5;
    double c_jet=0.12;
    double k_jet = 2*pi/L_jet;
    double Fi_jet = pi/2;

    double sigma = 0.05;

    //double theta_jet=pi/2.;

    double number_cells=1;
    double Ly=7.5;
    double Lxi=Lx;
    double Lyi=Ly;
    double xi=Lx/2.;
    double yi=Ly/2.;
    double factor=1;
    double delta;
    double NNN_forward;
    double D_factor=1;

    //double resc_const = (D_factor/NNN)/(3*pow(10, -9));          //sqrt(3)*pow(10, -3);



    ///When running the program
    for (int i = 1; i < argc; i++) {
        if(string(argv[i]) == "-Ns"){
            NNN_forward = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Lx"){
            Lx = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Ly"){
            Ly = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-Lxi"){
            Lxi = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-xi"){
            xi = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-yi"){
            yi = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-mutation_prob"){
            mutation_probability = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-random_seed"){
            random_seed = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-networks_to_average"){
            NETWORKS_TRIALS = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-networks_to_print"){
            NETWORKS_TO_PRINT = atoi(argv[i + 1]);
        }  else if(string(argv[i]) == "-trials_per_network"){
            ABUNDANCES_TRIALS = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-rt") {
            refresh_time = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-B0_jet"){
            B0_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-epsilon_jet"){
            epsilon_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-w_jet"){
            w_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-l_jet"){
            L_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-c_jet"){
            c_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-theta_jet"){
            k_jet = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-number_cells"){
            number_cells = atoi(argv[i + 1]);
        } else if(string(argv[i]) == "-factor"){
            factor = atof(argv[i + 1]);
        } else if(string(argv[i]) == "-D_factor"){
            D_factor = atof(argv[i + 1]);
        } else if(i%2!=0){
            cout<<endl<<" CAUTION!: "<<string(argv[i])<<" HAVE NOT SET ANY PARAMETER!"<<endl;
        }
    }


    ///////////////////////// VARIABLES ///////////////////////
    //cout<<"Lx is"<<Lx<<"; "<<"Ly is "<<Ly<<endl;
    time_t itime=time(NULL),ftime,TOTAL_itime=time(NULL),TOTAL_ftime;
    int trial;
    #if defined(TWO_DIMENSIONAL)
    double space_between_areas=Lx/10.;
    Lyi=Lxi;
    delta=(int)floor(sqrt(16384.))/Lx;   //double oneoverd = (int)floor(sqrt(NNN))/double(Lx);
    //cout<<"delta is"<<delta<<endl;
    double DIFF= (D_factor*Lx*Lx)/double(16384.); //3*pow(10,-9)/pow(resc_const,2);
    //cout<<"DIFF is"<<DIFF<<endl;
    double NNN=NNN_forward;
    #elif defined(ONE_DIMENSIONAL)
    double space_between_areas=Lx/10.;
    Lyi=Lxi;
    delta=(int)floor(NNN_forward);
    double DIFF=(D_factor*Lx*Lx)/double(NNN_forward*NNN_forward);
    double NNN=NNN_forward*double(Lxi*Lxi);
    #endif

    #ifdef MUTATION_IN_DYNAMICS
    ABUNDANCES_TRIALS=1;
    #endif

    double jet_parameters[8];
    //double g_jet=0.001;//1;
    //double l_jet=2*M_PI*L_jet;
    //double a_jet=1;
    //double v0_jet=1;
    //double k_jet=1;
    //w_jet=30;
    jet_parameters[0]=B0_jet;
    jet_parameters[1]=epsilon_jet;
    jet_parameters[2]=Fi_jet;
    jet_parameters[3]=w_jet;
    jet_parameters[4]=c_jet;
    jet_parameters[5]=L_jet;
    jet_parameters[6]=k_jet;
    jet_parameters[7]=(double)(Lx*Lx/16384.)/(sigma*sigma/2);  //rescaling const
    double interaction_range=1./delta;
    //cout<<"Inter range is"<<interaction_range<<endl;
    number_of_spatial_radius=int((Lx+Lx/2)/interaction_range)+1;//int(1.5/interaction_range)+1;
    //cout<<"Numb of spat rad = "<<number_of_spatial_radius<<endl;
    coalescing_random_walk_tree current_tree(NNN,mutation_probability,random_seed,space_between_areas,number_cells,DIFF,delta,number_of_spatial_radius,Lx,Ly,jet_parameters,Lxi,Lyi,xi,yi);
    cout<<endl<<endl;
    cout<<"*******************************************************************************************"<<endl;
    cout<<"***                     YOU HAVE SELECTED DUAL VOTER MODEL IN A FLUX                    ***"<<endl;
    cout<<"*******************************************************************************************"<<endl;
    cout<<endl<<endl;

    ///////////////////////// INITIALIZATION ///////////////////////
    //INITIALIZE_RANDOM_SEED(rand_seed):
    //cout<<"Begin to init rand numb"<<endl;
    int rand_seed=random_seed;
    const gsl_rng_type * gsl_rng_T;	//Random variable
    gsl_rng * random_variable;	//Random variable
    gsl_rng_env_setup();
    gsl_rng_default_seed=rand_seed;
    gsl_rng_T=gsl_rng_default;
    random_variable=gsl_rng_alloc(gsl_rng_T);
    if(refresh_time>1) refresh_time=(refresh_time+gsl_rng_uniform_int(random_variable,5))*60;
    //cout<<"End to rnd numb init"<<endl;

    ofstream spatial_data;
    char spatial_name[300];
    sprintf(spatial_name,"Spatial_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
    bool print_spatial=true;
    spatial_data.open(spatial_name, ios::out | ios::trunc);

    ofstream current_events_data;
    char current_events_name[300];
    sprintf(current_events_name,"current_events_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
    current_events_data.open(current_events_name, ios::out | ios::trunc);


    ofstream coal_part_data;
    char coal_part_name[300];
    sprintf(coal_part_name,"Numb_Coal_Part_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
    coal_part_data.open(coal_part_name, ios::out | ios::trunc);


    ///////////////////////// NETWORKS LOOP ///////////////////////
    double mean_final_population=0;
    current_tree.network=0;
    current_tree.input_ncoal(factor);
    while(current_tree.network<NETWORKS_TRIALS){

#ifdef GET_INITIAL_SPATIAL_DISTRIBUTION
        int my_network=floor(current_tree.network/100.)+1;
        current_tree.set_initial_spatial_data(D_factor,my_network);
#endif
	//cout<<"Begin generate_network"<<endl;
        current_tree.generate_network(random_variable, coal_part_data, &mean_final_population);
	//cout<<"End gener netw"<<endl;

        #ifdef PRINT_SAD
        if(ABUNDANCES_TRIALS!=0){
            ///////////////////////// OBTAINING SAD AND SRC ///////////////////////
            trial=0;
            cout<<"OBTAINING SAD AND SRC"<<endl;
            while(trial<ABUNDANCES_TRIALS){
		cout<<"while loop"<<endl;
                current_tree.assign_species(random_variable);
		cout<<"1"<<endl;
                current_tree.get_abundances(trial,Lxi,D_factor);
		cout<<"2"<<endl;
		current_tree.get_measures_subsamples(Lxi,trial,interaction_range,print_spatial,&spatial_data);
		cout<<"3"<<endl;
        	current_tree.print_abundances(D_factor,Lxi);
		cout<<"4"<<endl;
		print_spatial=false;
                cout<<"5"<<endl;
		trial++;
            }
        }

        ofstream alpha_data;
        char alpha_name[300];
        sprintf(alpha_name,"Alpha-diversity_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
        alpha_data.open(alpha_name, ios::out | ios::trunc);
        for(int i=0;i<current_tree.mean_alpha.size();i++){
            if(current_tree.mean_alpha[i]!=0){
                alpha_data<<i*interaction_range<<" ";
                alpha_data<<current_tree.mean_alpha[i]/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
                alpha_data<<ABUNDANCES_TRIALS<<" "<<current_tree.network<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
            }
        }
        alpha_data.close();

        ofstream heterozigosity_data;
        char heterozigosity_name[300];
        sprintf(heterozigosity_name,"Beta-diversity_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
        heterozigosity_data.open(heterozigosity_name, ios::out | ios::trunc);
        for(int i=0;i<current_tree.mean_heterozigosity_all_pairs.size();i++){
            if(current_tree.mean_heterozigosity_all_pairs[i]!=0){
                heterozigosity_data<<i*interaction_range<<" ";
                heterozigosity_data<<current_tree.mean_heterozigosity_diff_pairs[i]/double(current_tree.mean_heterozigosity_all_pairs[i])<<" ";
                heterozigosity_data<<ABUNDANCES_TRIALS<<" "<<current_tree.network<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
            }
        }
        heterozigosity_data.close();

        ofstream tdistribution_data;
        char tdistribution_name[300];
        sprintf(tdistribution_name,"tdistribution_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
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
        sprintf(total_events_name,"total_events_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
        total_events_data.open(total_events_name, ios::out | ios::trunc);
        total_events_data<<current_tree.coalescences_total/double(ABUNDANCES_TRIALS*current_tree.network);
        total_events_data<<" "<<current_tree.mutations_total/double(ABUNDANCES_TRIALS*current_tree.network)<<endl;
        total_events_data.close();

        current_events_data.close();
        current_events_data.open(current_events_name, ios::out | ios::app);
        current_events_data<<current_tree.current_mutations<<" "<<double(ABUNDANCES_TRIALS*current_tree.network)<<" "<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        //current_tree.current_mutations=0;
        current_events_data.close();

        ofstream final_separation_data;
        char final_separation_name[300];
        sprintf(final_separation_name,"final_separation_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_probability,int(NNN_forward),Lx,Lxi,random_seed);
        final_separation_data.open(final_separation_name, ios::out | ios::trunc);
        final_separation_data<<current_tree.mean_final_individuals_distance/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.mean_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.min_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<current_tree.max_final_individuals_distance_all/double(ABUNDANCES_TRIALS*current_tree.network)<<" ";
        final_separation_data<<ABUNDANCES_TRIALS*current_tree.network<<endl;
        final_separation_data.close();

	//cout<<"Closed final separation"<<endl;

        #endif
	//cout<<"After endif"<<endl;

        current_tree.clear();

	//cout<<"After current"<<endl;
        TOTAL_ftime=time(NULL);

	//cout<<"After total"<<endl;

	cout<<"Mean number of border colisions "<<current_tree.border_colisions/current_tree.network<<" (periodice cells="<<number_cells<<")"<<endl;
        cout<<"NETWORK "<<current_tree.network<<" ENDED. TOTAL RUNNING TIME="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl<<endl;
    }

    return 0;
}
