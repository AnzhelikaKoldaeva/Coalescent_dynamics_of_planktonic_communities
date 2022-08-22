//
// Created by paula on 17/11/15.
//

#ifndef TREE_CLASSES_H
#define TREE_CLASSES_H

#include "include.h"



const double DT=0.001;//0.001;
const double TWO_PI=2*M_PI;
float MAX_MBITS=100;
//string MODEL_NAME="";  /*mistake here? --- added it from the classes.h - backward code*/

using namespace std;

bool IsnotO(int i) { return (i!=0); }
bool descending_order (int i,int j) { return (i>j); }
void process_mem_usage(double& vm_usage, double& resident_set) {
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage     = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned int vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}


struct network_element{
    vector <unsigned int> neighbors;
    bool occupied;
    unsigned int vector_position;
    network_element(){ occupied=1;}
};
struct network_class{

    vector <network_element> sites_vector;

    void link_1D_periodic(int N){

        network_element network_element_new;
        network_element_new.neighbors.push_back(0);
        network_element_new.neighbors.push_back(0);
        sites_vector.clear();

        network_element_new.neighbors[0]=N-1;
        network_element_new.neighbors[1]=1;
        sites_vector.push_back(network_element_new);
        for(int i=1;i<N-1;i++){
            network_element_new.neighbors[0]=i+1;
            network_element_new.neighbors[1]=i-1;
            sites_vector.push_back(network_element_new);
        }
        network_element_new.neighbors[0]=N-2;
        network_element_new.neighbors[1]=0;
        sites_vector.push_back(network_element_new);

    };

    void link_2D_periodic(int N){

        int L=sqrt(N);
        network_element network_element_new;
        network_element_new.neighbors.push_back(0);
        network_element_new.neighbors.push_back(0);
        network_element_new.neighbors.push_back(0);
        network_element_new.neighbors.push_back(0);

        sites_vector.clear();
        for(int isite=0;isite<L;isite++){
            for(int jsite=0;jsite<L;jsite++){
                network_element_new.neighbors[0]=((isite+1+L)%L)*L+jsite;
                network_element_new.neighbors[1]=((isite-1+L)%L)*L+jsite;
                network_element_new.neighbors[2]=isite*L+(jsite+1+L)%L;
                network_element_new.neighbors[3]=isite*L+(jsite-1+L)%L;
                sites_vector.push_back(network_element_new);
            }
        }
    }
};
struct tree_individual{
    list<tree_individual>::iterator parent;
    unsigned int own_dir,label, intermediate_births,intermediate_births2,offspring_number;
    bool active,print,mutation,species_assignated;
    //double t_birth;
    list<int>::iterator species_abundance;

    tree_individual(){
        print=false;
        mutation=false;
        species_assignated=false;
        intermediate_births=0;
        intermediate_births2=0;
        active=false;
        offspring_number=0;
        //t_birth=0;
    }

    void initialize(){
        mutation=false;
        species_assignated=false;
    }
};
struct position_element{
    double prop;
    tree_individual *dir;

    bool operator<(const  position_element & other){ return prop < other.prop; }
};
struct ascending_prop{
    inline bool operator() (const position_element& struct1, const position_element& struct2){
        return (struct1.prop < struct2.prop);
    }
};
struct ancestor_individual{

    unsigned int alive_offspring;
    list <tree_individual> offspring;

    ancestor_individual(){
        alive_offspring=1;
    };

};
struct dynamics_individual{

    list<ancestor_individual>::iterator ancestor;
    list <tree_individual>::iterator offspring_it;
};
struct ascending_label{
    inline bool operator() (const dynamics_individual& struct1, const dynamics_individual& struct2){
        return ((*struct1.offspring_it).label < (*struct2.offspring_it).label);
    }
};
struct descending_label{
    inline bool operator() (const dynamics_individual& struct1, const dynamics_individual& struct2){
        return ((*struct1.offspring_it).label > (*struct2.offspring_it).label);
    }
};
struct flux_individual:public dynamics_individual{
    double x,y;
    double fieldx,fieldy;
    double fieldoldx,fieldoldy;

    void initial_set(gsl_rng *r,double force, double *p= nullptr){


#ifdef JET_FLOW
	double Lx = 7.5;
	double Ly = 7.5;
        x=gsl_rng_uniform(r)*Lx;
        y=gsl_rng_uniform(r)*Ly;
        advectx(force,0,0,p,0,0,0);
        advecty(force,0,0,p,0,0,0);
#else
	double Lx = 7.5;
	double Ly = 7.5;
        x=gsl_rng_uniform(r)*Lx;
        y=gsl_rng_uniform(r)*Ly;
        advectx(force);
        advecty(force);
#endif
        fieldoldx=fieldx;
        fieldoldy=fieldy;
    }
    void advectx(double force=0, double t=0, double Ly_particular=1, double *p= nullptr, double Kxi=0, double Kyi=0, double ti=0){

        #ifdef ADVECTION
        //double Ly=7.5, Lx=7.5;
        double flux_y=((y+Kyi*ti)-floor((y+Kyi*ti)/Ly)*Ly)-Ly/2.;
        double flux_x=((x+Kxi*ti)-floor((x+Kxi*ti)/Lx)*Lx)-Lx/2.;
        //double gamma=p[0];
        //double l=p[1];
        //double alpha=p[2];
        //double v0=p[3];
        //double k1=p[4];
        //double w1=p[5];
        //double arg=2*M_PI/l;
        //double resc_const = p[6];//(D_factor/(NNN))/pow(10, -9);
        //double fieldx0=gamma/l*(sinh(flux_y*arg)/(cosh(flux_y*arg)-alpha*cos(flux_x*arg)));
        //double fieldx_pert=v0*sin(2*k*flux_x-2*w*t);
        //fieldx_pert=0;
        //fieldx=(fieldx0+fieldx_pert);

        //double B0 = 1.2;
      	//double eps =0.3;// 0.3;
      	//double Fi = M_PI/2;
      	//double w = 0.4;
      	//double c = 0.12;
        //double L = 7.5;
        //double k = 2*M_PI/L;

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
	//cout<<"(t, Bt, tanh_arg,  resc_const, x, y, fieldx) "<<t<<" "<<Bt<<" "<<tanh_arg<<" "<<resc_const<<" "<<flux_x<<" "<<flux_y<<" "<<fieldx<<endl;

	//cout<<"fieldx is calculated"<<endl;
    }
    void advecty(double force=0, double t=0, double Lx_particular=1, double *p= nullptr, double Kxi=0, double Kyi=0, double ti=0){

        #ifdef ADVECTION
        //double Ly=7.5, Lx=7.5;
        double flux_y=((y+Kyi*ti)-floor((y+Kyi*ti)/Ly)*Ly)-Ly/2.;
        double flux_x=((x+Kxi*ti)-floor((x+Kxi*ti)/Lx)*Lx)-Lx/2.;
	//cout<<"(x, y) is "<<x<<", "<<y<<"; (flux_x, flux_y) is "<<flux_x<<"; "<<flux_y<<endl;
        //double gamma=p[0];
        //double l=p[1];
        //double alpha=p[2];
        //double v0=p[3];
        //double k1=p[4];
        //double w1=p[5];
        //double arg=2*M_PI/l;

        //double B0 = 1.2;
      	//double eps = 0.3;//0.3;
      	//double Fi = M_PI/2;
      	//double w = 0.4;
      	//double c = 0.12;
        //double L = 7.5;
        //double k = 2*M_PI/L;
        //double resc_const = p[6];

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

        //double fieldy0=-gamma*alpha/l*(sin(flux_x*arg)/(cosh(flux_y*arg)-alpha*cos(flux_x*arg)));
        fieldy= -(1-pow(tanh(tanh_arg),2))*k*Bt*sin(k*flux_x)*(expr-k*k*Bt*cos(k*flux_x)*(flux_y-Bt*cos(k*flux_x)))/(sqrt(expr)*expr)*sqrt(resc_const);
	//fieldy= -(1-pow(tanh(tanh_arg),2)*k*Bt*sin(k*flux_x)*(expr-k*k*Bt*cos(k*flux_x)*(flux_y-Bt*cos(k*flux_x)))/(sqrt(expr)*expr))*sqrt(resc_const);
        #else
        fieldy=0;
        #endif
	//cout<<"fieldy is calculated"<<endl;
        //cout<<"(t, Bt, tanh_arg, resc_const, x, y, fieldy) "<<t<<" "<<Bt<<" "<<tanh_arg<<" "<<resc_const<<" "<<flux_x<<" "<<flux_y<<" "<<fieldy<<endl;
    }
};
struct network_individual:public dynamics_individual{
    vector<network_element>::iterator network_site;
};
struct species_active_individual{
    list <tree_individual>::iterator individual_it;
};
float tree_individual_struct_size=sizeof(tree_individual)*1e-6;

class tree{

private:
    bool previous_print;
    unsigned int MAX_NUMBER_OF_INDIVIDUALS,TOTAL_TRIALS;
    vector <double> SR_histogram,SAD_histogram;
    list <int> species_abundance;
    fstream SAD,SAD_subsample,SR,generation_tree_ACTIVES,generation_tree_FAMILIES,generation_tree_PROPERTIES;
    time_t itime,ftime,TOTAL_itime,TOTAL_ftime, time_one_netw;
    list <ancestor_individual>::iterator it_ancestor;

    virtual void initialize_network(){}
    virtual void initialize(gsl_rng * random_variable, double *p= nullptr){

        tree_individual new_offspring_individual;

#ifdef FLUX
        flux_individual new_individual;
#elif defined(NETWORK)
        network_individual new_individual;
        #else
            dynamics_individual new_individual;
#endif

        ancestor_individual new_ancestor;

        last_label=0;
        for(int i=0;i<NNN;i++){

            ancestors_list.push_back(new_ancestor);
            //The first individual is "0": parent of ancestors
            new_offspring_individual.label=0;
            ancestors_list.back().offspring.push_back(new_offspring_individual);
            ancestors_list.back().offspring.back().parent=--ancestors_list.back().offspring.end();

            last_label++;
            new_offspring_individual.parent=--ancestors_list.back().offspring.end();
            new_offspring_individual.label=last_label;
            ancestors_list.back().offspring.push_back(new_offspring_individual);

            new_individual.initial_set(random_variable,force,p);
            new_individual.ancestor=--ancestors_list.end();
            new_individual.offspring_it=--ancestors_list.back().offspring.end();
            active_individuals_vector.push_back(new_individual);
        }

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
    unsigned int before_individuals,after_individuals,random_seed,NNN,last_label;
    double mutation_rate,force,Lx,Ly,t,dt,tnext;
    char SAD_name[1280],SAD_subsample_name[1280],SR_name[1280],generation_tree_name_ACTIVES[1280],generation_tree_name_FAMILIES[1280],generation_tree_name_PROPERTIES[1280],spacetime_name[1280],fname[1280];
    list <tree_individual>::iterator offspring_begin_it,offspring_end_it,it_offspring1,it_offspring2;
    string model;
    vector <flux_individual> active_individuals_vector;
    flux_individual new_active_individual;
    list <ancestor_individual> ancestors_list;

    tree(int N,double mutprob,int r,double number_of_spatial_radius_in=0,double Lx_in=-1,double Ly_in=-1):NNN(N),mutation_rate(mutprob),random_seed(r),number_of_spatial_radius(number_of_spatial_radius_in),Lx(Lx_in),Ly(Ly_in){
        printed_networks=0;
        network=0;
        TOTAL_TRIALS=0;
        last_label=0;

        MAX_NUMBER_OF_INDIVIDUALS=NNN+1;
        SR_histogram.assign(MAX_NUMBER_OF_INDIVIDUALS,0);
        SAD_histogram.assign(MAX_NUMBER_OF_INDIVIDUALS,0);
        SAD_histogram_subsamples_mean.assign(MAX_NUMBER_OF_INDIVIDUALS,0);
        itime=time(NULL);
        TOTAL_itime=time(NULL);
        model=MODEL_NAME;

        mean_number_individuals_subarea=0;
        mean_number_species_subarea=0;
        mean_number_species=0;

        mean_heterozigosity_diff_pairs.assign(number_of_spatial_radius,0);
        mean_alpha.assign(number_of_spatial_radius,0);
        mean_heterozigosity_all_pairs.assign(number_of_spatial_radius,0);

        tdmax=1e4;
        tdistribution_points=int(tdmax)+1;
        coalescences_tdistribution.assign(tdistribution_points,0);
        mutations_tdistribution.assign(tdistribution_points,0);
        events_tdistribution.assign(tdistribution_points,0);
        coalescences_total=0;
        mutations_total=0;
        mean_interacting_individuals=0;
        mean_final_individuals_distance=0;
        mean_final_individuals_distance_all=0;
        min_final_individuals_distance_all=0;
        max_final_individuals_distance_all=0;
    }
    virtual void one_step_dynamics(gsl_rng * random_variable){}

public:
    int number_of_subareas;
    int number_of_occupied_subareas;
    struct subsample{
        int number_species;
        int number_individuals;
        vector <int> species_abundance_vector;
        bool occupied=false;
    };

    double mean_number_species, mean_number_species_subarea, mean_number_individuals_subarea,mean_beta_diversity;
    vector <double> SAD_histogram_subsamples_mean;
    vector <double> mean_heterozigosity_diff_pairs;
    vector <double> mean_heterozigosity_all_pairs;
    vector <double> mean_alpha;
    int number_of_spatial_radius;
    vector<double> coalescences_tdistribution,mutations_tdistribution,events_tdistribution;
    long unsigned int tdistribution_points;
    long unsigned int coalescences_total,mutations_total,current_mutations;
    double tdmax;
    double mean_final_individuals_distance;
    double mean_final_individuals_distance_all;
    double min_final_individuals_distance_all;
    double max_final_individuals_distance_all;
    double max_final_individuals_distance_all_network;
    double mean_interacting_individuals;
    vector <int> actives_index_vector;
    vector <int> subsample_index;
    vector <int> subsample_index_muted;
    vector <int> t_muted_species,i_muted_species;
    int printed_networks,network;
    void remove_tree_individual(dynamics_individual *individual){
        (*(*individual).ancestor).alive_offspring--;
        if((*(*individual).ancestor).alive_offspring<1){
            before_individuals-=(*(*individual).ancestor).offspring.size();
            (*individual).ancestor=ancestors_list.erase((*individual).ancestor);
        }
    };
    void add_tree_individual(dynamics_individual *son,dynamics_individual *father){

        tree_individual new_offspring;
        (*son)=(*father);
        last_label++;
        new_offspring.parent=(*father).offspring_it;
        new_offspring.label=last_label;
        (*(*father).ancestor).offspring.push_back(new_offspring);
        (*son).offspring_it=--(*(*father).ancestor).offspring.end();
        before_individuals++;
        (*(*son).ancestor).alive_offspring++;
    }
    void remove_intermediate_individuals(){

        list <tree_individual>::iterator previous_it,current_it;
        ///FIRSTLY WE COUNT OFFSPRING OF EACH INDIVIDUAL
        it_ancestor=active_individuals_vector[0].ancestor;
        offspring_begin_it=++(*it_ancestor).offspring.begin();
        offspring_end_it=(*it_ancestor).offspring.end();
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            (*it_offspring1).print=false;
            (*(*it_offspring1).parent).offspring_number++;
        }

        ///SECONDLY WE ASSIGN PARENTS AND INTERMEDIATE BIRTHS
        //cout<<"ACTIVES: ";
        for(int i=0;i<active_individuals_vector.size();i++){
            (*active_individuals_vector[i].offspring_it).active=true;
            (*active_individuals_vector[i].offspring_it).print=true;
            //cout<<(*active_individuals_vector[i].offspring_it).label<<" ";
        }
        //cout<<endl;

        offspring_begin_it=++(*it_ancestor).offspring.begin();
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            previous_it=it_offspring1;
            current_it=(*it_offspring1).parent;

            while(!(((*current_it).offspring_number>1)||((*current_it).active)||((*current_it).label==0))) {
                previous_it=current_it;
                current_it=(*current_it).parent;
            }

            (*it_offspring1).parent=current_it;
            (*it_offspring1).intermediate_births+=(*previous_it).intermediate_births;
            (*it_offspring1).intermediate_births2+=(*previous_it).intermediate_births2;
            (*it_offspring1).intermediate_births++;
            (*current_it).print=true;
        }

#ifdef OUTPUT_NETWORKS
        offspring_begin_it=++(*it_ancestor).offspring.begin();
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            cout<<"AFTER SIMPLIFICATION: "<<(*it_offspring1).label<<" "<<(*(*it_offspring1).parent).label<<" "<<(*it_offspring1).intermediate_births<<endl;
        }
#endif

        ///FINALLY WE REMOVE INTERMEDIATE VALUES
        it_offspring1=++(*it_ancestor).offspring.begin();
        while(it_offspring1!=(*it_ancestor).offspring.end()){
            if(!(*it_offspring1).print) it_offspring1=(*it_ancestor).offspring.erase(it_offspring1);
            else it_offspring1++;
        }

#ifdef OUTPUT_NETWORKS
        int j=0;
        offspring_begin_it=++(*it_ancestor).offspring.begin();
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            cout<<"AFTER REMOVING: "<<(*it_offspring1).label<<" "<<(*(*it_offspring1).parent).label<<" "<<(*it_offspring1).intermediate_births<<endl;
            j++;
        }
        cout<<"AFTER REMOVING INDIVIDUALS N="<<j+1<<" "<<NNN<<endl;
#endif

    }
    void time_increment(gsl_rng * random_variable){dt=DT;};
    int get_Nactives(){ return active_individuals_vector.size();}
    virtual void generate_network(gsl_rng * random_variable, double *mean_final_population=0,double *p= nullptr){
        bool overtime=false;
        double max_time=100000;
        int current_time=0;
        int max_size=4000;


        initialize(random_variable,p);
        cout<<"CREATING NETWORK "<<network+1<<" IN "<<model.c_str()<<endl;
        t=0.0;
        time_increment(random_variable);
        TOTAL_ftime=time(NULL);
	time_one_netw = time(NULL);
        cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
        do{
	    //cout<<"Begin the cycle; t = "<<t<<endl;
            t+=dt;
            ///TIME UPDATE
            ftime=time(NULL);
            TOTAL_ftime=time(NULL);
            if(difftime(ftime, itime)>=(2)*60){
                cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
                itime=time(NULL);
                //cout<<current_time<<endl;
            }
	    //cout<<"t = "<<t<<endl;
            one_step_dynamics(random_variable);
            if(before_individuals*tree_individual_struct_size>MAX_MBITS) eliminate_dead_branches();
	    //cout<<"list.back "<<ancestors_list.back().alive_offspring<<"; vector.size "<<active_individuals_vector.size()<<endl;
	    //if(difftime(TOTAL_ftime,time_one_netw)/60. > 120) {
	    //overtime = true;
	    //cout<<"Overtime!"<<endl;
	    //}
        }while((ancestors_list.back().alive_offspring!=active_individuals_vector.size())&&(!overtime));

        eliminate_dead_branches();
        if(overtime){
            auto it_ancestor_begin=ancestors_list.begin();
            auto it_ancestor_end=ancestors_list.end();
            --it_ancestor_end;
            auto it_common_ancestor=it_ancestor_end;
            tree_individual new_offspring,other_ancestor_tree_ind;
            list <tree_individual>::iterator it_common_ancestor_tree_ind,it_other_ancestor_tree_ind,new_parent_it,new_it;
            int position,last_position;
            vector <int> ancestor_size;
            int i=0;
            it_common_ancestor_tree_ind=(*it_common_ancestor).offspring.begin();
            it_common_ancestor_tree_ind++;

            last_position=0;
            i=0;
            for(auto it_ancestor=ancestors_list.begin();it_ancestor!=it_ancestor_end;++it_ancestor){

                last_position=(*it_common_ancestor).offspring.size();
                ancestor_size.push_back(last_position);

                it_other_ancestor_tree_ind=(*it_ancestor).offspring.begin();
                it_other_ancestor_tree_ind++;

                other_ancestor_tree_ind.intermediate_births=1;
                other_ancestor_tree_ind.parent=(*it_common_ancestor).offspring.begin();
                other_ancestor_tree_ind.label=(*it_other_ancestor_tree_ind).label;
                (*it_common_ancestor).offspring.push_back(other_ancestor_tree_ind);
                offspring_begin_it=(*it_ancestor).offspring.begin();
                ++offspring_begin_it;
                ++offspring_begin_it;
                offspring_end_it=(*it_ancestor).offspring.end();

                for (it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
                    new_offspring.intermediate_births=(*it_offspring1).intermediate_births+1;
                    position=distance((*it_ancestor).offspring.begin(),(*it_offspring1).parent);
                    new_parent_it=(*it_common_ancestor).offspring.begin();
                    advance(new_parent_it,last_position+position-1);
                    new_offspring.parent=new_parent_it;
                    new_offspring.label=(*it_offspring1).label;
                    (*it_common_ancestor).offspring.push_back(new_offspring);
                }
                i++;
            }

            int ancestor_number;
            for(int i=0;i<active_individuals_vector.size();i++){
                if(active_individuals_vector[i].ancestor!=it_common_ancestor){

                    ancestor_number=distance(ancestors_list.begin(),active_individuals_vector[i].ancestor);
                    new_it=(*it_common_ancestor).offspring.begin();
                    position=distance((*active_individuals_vector[i].ancestor).offspring.begin(),active_individuals_vector[i].offspring_it);
                    advance(new_it,ancestor_size[ancestor_number]+position-1);
                    active_individuals_vector[i].offspring_it=new_it;
                    active_individuals_vector[i].ancestor=it_common_ancestor;
                }
            }
            int extra_ancestors=ancestors_list.size()-1;
            for(int i=0;i<extra_ancestors;i++) ancestors_list.pop_front();
            cout<<"TREE RECONSTRUCTED!"<<endl;

            offspring_begin_it=++(*it_common_ancestor).offspring.begin();
            offspring_end_it=(*it_common_ancestor).offspring.end();
        }
        eliminate_dead_branches();
        remove_intermediate_individuals();
        cout<<"t="<<t<<", Number of ancestors="<<ancestors_list.size()<<", Number of alive_individuals="<<active_individuals_vector.size()<<", Running time="<<(difftime(TOTAL_ftime, TOTAL_itime))/60.<<" min"<<endl;
        (*mean_final_population)*=network;
        (*mean_final_population)+=active_individuals_vector.size();
        network++;
        (*mean_final_population)/=network;
    }
    void assign_species(gsl_rng * random_variable){

        double mutation_probability;
        it_ancestor=ancestors_list.begin();
        offspring_begin_it=(*it_ancestor).offspring.begin();
        offspring_end_it=(*it_ancestor).offspring.end();
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++) (*it_offspring1).initialize();

        ///////////////////////// MUTATIONS ///////////////////////
        ++offspring_begin_it;
        double t_local;
        for(it_offspring1=offspring_begin_it;it_offspring1!=offspring_end_it;it_offspring1++){
            mutation_probability=1-pow((1-mutation_rate),(*it_offspring1).intermediate_births);
            if(gsl_rng_uniform(random_variable)<mutation_probability) (*it_offspring1).mutation=true;
        }

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
            if((*current_it).species_assignated){///IF IT IS ASSIGNED
                for(int j=0;j<to_assign_vector.size();j++){
                    (*to_assign_vector[j].individual_it).species_abundance=(*current_it).species_abundance;
                    (*to_assign_vector[j].individual_it).species_assignated=true;
                }
            }
            else{///IF IT IS MUTATED OR IS THE ANCESTOR
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
        for(int i=0;i<active_individuals_vector.size();i++) (*(*active_individuals_vector[i].offspring_it).species_abundance)++;
    };
    void get_abundances(int trial){

        list<int>::iterator it;
        for(it=species_abundance.begin();it!=species_abundance.end();it++) {
            SAD_histogram[(*it)]++;
            if((*it)!=0) mean_number_species++;
        }
        TOTAL_TRIALS++;
    }
    void get_actives_index_vector(){
        actives_index_vector.clear();
        for(int i=0;i<active_individuals_vector.size();i++){
            actives_index_vector.push_back(i);
        }
    }
    void get_sample_index_vector(int NNN_sample,gsl_rng * random_variable){
        subsample_index.clear();
        subsample_index_muted.clear();
        subsample_index=actives_index_vector;
        subsample_index_muted=subsample_index;
    }
    void get_measures_nosubsamples(double Lxi,int trial,double interaction_range,bool print_spatial,ofstream *spatial_data,int NNN_sample,gsl_rng * random_variable){

        list<int>::iterator it;
        number_of_subareas=1;
        vector <subsample> subsample_vector;
        subsample subsample_new_element;
        int current_cell;
        int species_label, species_label2;
        int number_of_species;
        double current_distance,alternative_distance;
        double current_distance_x;
        double current_distance_y;
        double alternative_distance_x;
        double alternative_distance_y;

        int middle_cell1=int(0.5*Lx/double(1))+Lx/double(0.5)*int(0.5*Lx/double(1)); //int(0.5/double(1))+Lx/double(0.5)*int(0.5/double(1));
        vector <int> small_index;
        bool *checked_species_alpha;
        int *current_alpha;
        int number_of_spatial_radius=int(1.5*Lx/interaction_range);

        checked_species_alpha=(bool *) malloc(MAX_NUMBER_OF_INDIVIDUALS*sizeof(bool));
        current_alpha=(int *) malloc(number_of_spatial_radius*sizeof(int));
        for(int k=0;k<MAX_NUMBER_OF_INDIVIDUALS;k++) checked_species_alpha[k]=false;
        for(int k=0;k<number_of_spatial_radius;k++) current_alpha[k]=0;

        ///Here I construct the vector of subsamples
        vector <int> species_abundance_new_list;
        vector <int> abundance_sample;
        for(int j=0;j<species_abundance.size();j++) species_abundance_new_list.push_back(0);
        subsample_new_element.number_species=0;
        subsample_new_element.number_individuals=0;
        subsample_new_element.occupied=false;
        subsample_new_element.species_abundance_vector=species_abundance_new_list;
        subsample_vector.push_back(subsample_new_element);

        ///Here I count the species abundance at each subarea
        int max_species_label=0;
        if(subsample_index.size()>NNN_sample){

            int random_index;
            int diff_N=subsample_index.size()-NNN_sample;
            for(int i=0;i<diff_N;i++){
                random_index=gsl_rng_uniform_int(random_variable,subsample_index.size());
                subsample_index.erase (subsample_index.begin()+random_index);
            }
            subsample_index_muted=subsample_index;

            max_species_label=0;
            for(int i_pre=0;i_pre<subsample_index.size();i_pre++){
                int i=subsample_index[i_pre];
                species_label=distance(species_abundance.begin(), (*active_individuals_vector[i].offspring_it).species_abundance);
                if(max_species_label<species_label) max_species_label=species_label;
                subsample_index_muted[i_pre]=species_label;
                small_index.push_back(i);
                subsample_vector[0].number_individuals++;
                subsample_vector[0].species_abundance_vector[species_label]++;
            }

            i_muted_species.clear();
            for(int i=0;i<max_species_label;i++) i_muted_species.push_back(-1);
            for(int i_pre=0;i_pre<small_index.size();i_pre++){
                coalescences_total++;
                i_muted_species[subsample_index_muted[i_pre]]=i_pre;
            }
            for(int ispecies=0;ispecies<max_species_label;ispecies++){
                if(i_muted_species[ispecies]!=-1){
                    mutations_total++;
                    current_mutations++;
                    coalescences_total--;
                }
            }

            ///Here I measure beta-diversity and heterozigozity of the two middle cells
            int all_pairs=0;
            int equal_pairs=0;
            double mean_final_individuals_distance_aux=0;
            double mean_final_individuals_distance_all_aux=0;
            double max_final_individuals_distance_all_aux=0;
            double min_final_individuals_distance_all_aux=10;
            for(int j=0;j<small_index.size();j++){
                int i=small_index[j];
                species_label=subsample_index_muted[j];
                if(print_spatial) (*spatial_data)<<active_individuals_vector[i].x<<" "<<active_individuals_vector[i].y<<" "<<species_label<<endl;
                for(int j2=j+1;j2<small_index.size();j2++){
                    int i2=small_index[j2];
                    species_label2=subsample_index_muted[j2];
                    current_distance_x=fabs(active_individuals_vector[i].x-active_individuals_vector[i2].x);
                    current_distance_y=fabs(active_individuals_vector[i].y-active_individuals_vector[i2].y);
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
                    current_distance=sqrt(pow(active_individuals_vector[i].x-0.5*Lx,2)+pow(active_individuals_vector[i].y-0.5*Lx,2));
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
            max_final_individuals_distance_all_network=max_final_individuals_distance_all_aux;
            int cummulated_alpha=0;
            for(int i=0;i<number_of_spatial_radius;i++){
                cummulated_alpha+=current_alpha[i];
                mean_alpha[i]+=cummulated_alpha;
            }

            ///Here I count the number of species at each subarea and average over all subsamples
            number_of_occupied_subareas=1;
            number_of_species=0;
            for(int j=0;j<subsample_vector[0].species_abundance_vector.size();j++){
                if(subsample_vector[0].species_abundance_vector[j]!=0){
                    SAD_histogram_subsamples_mean[subsample_vector[0].species_abundance_vector[j]]++;
                    number_of_species++;
                }
            }
            subsample_vector[0].number_species=number_of_species;

            mean_number_individuals_subarea+=subsample_vector[0].number_individuals/double(number_of_occupied_subareas);
            mean_number_species_subarea+=subsample_vector[0].number_species/double(number_of_occupied_subareas);
        }

        species_abundance.clear();
        free(current_alpha);
    }
    void get_measures_subsamples(double Lxi,int trial,double interaction_range,bool print_spatial,ofstream *spatial_data,int NNN_sample,gsl_rng * random_variable){

        list<int>::iterator it;
        number_of_subareas=1;
        vector <subsample> subsample_vector;
        subsample subsample_new_element;
        int current_cell;
        int species_label, species_label2;
        int number_of_species;
        double current_distance,alternative_distance;
        double current_distance_x;
        double current_distance_y;
        double alternative_distance_x;
        double alternative_distance_y;

        int middle_cell1=int(0.25*Lx/double(0.5))+Lx/double(0.5)*int(0.25*Lx/double(0.5)); // middle_cell1=int(0.25/double(0.5))+Lx/double(0.5)*int(0.25/double(0.5));

        vector <int> small_index;

        bool *checked_species_alpha;
        int *current_alpha;
        int number_of_spatial_radius=int(1.5*Lx/interaction_range);

        checked_species_alpha=(bool *) malloc(MAX_NUMBER_OF_INDIVIDUALS*sizeof(bool));
        current_alpha=(int *) malloc(number_of_spatial_radius*sizeof(int));
        for(int k=0;k<MAX_NUMBER_OF_INDIVIDUALS;k++) checked_species_alpha[k]=false;
        for(int k=0;k<number_of_spatial_radius;k++) current_alpha[k]=0;

        ///Here I construct the vector of subsamples
        vector <int> species_abundance_new_list;
        vector <int> abundance_sample;
        for(int j=0;j<species_abundance.size();j++) species_abundance_new_list.push_back(0);
        subsample_new_element.number_species=0;
        subsample_new_element.number_individuals=0;
        subsample_new_element.occupied=false;
        subsample_new_element.species_abundance_vector=species_abundance_new_list;
        subsample_vector.push_back(subsample_new_element);

        ///Here I count the species abundance at each subarea
        int max_species_label=0;
        for(int i_pre=0;i_pre<subsample_index.size();i_pre++){
            int i=subsample_index[i_pre];
            current_cell=int(active_individuals_vector[i].x/double(0.5))+Lx/double(0.5)*int(active_individuals_vector[i].y/double(0.5));
            if(current_cell==middle_cell1){
                small_index.push_back(i);
            }
        }

        if(small_index.size()>NNN_sample){
            subsample_index.clear();
            subsample_index_muted.clear();
            subsample_index=small_index;

            int random_index;
            int diff_N=small_index.size()-NNN_sample;
            for(int i=0;i<diff_N;i++){
                random_index=gsl_rng_uniform_int(random_variable,subsample_index.size());
                subsample_index.erase (subsample_index.begin()+random_index);
            }

            subsample_index_muted=subsample_index;
            small_index.clear();
            max_species_label=0;
            for(int i_pre=0;i_pre<subsample_index.size();i_pre++){
                int i=subsample_index[i_pre];
                current_cell=int(active_individuals_vector[i].x/double(0.5))+Lx/double(0.5)*int(active_individuals_vector[i].y/double(0.5));
                species_label=distance(species_abundance.begin(), (*active_individuals_vector[i].offspring_it).species_abundance);
                if(max_species_label<species_label) max_species_label=species_label;
                subsample_index_muted[i_pre]=species_label;
                if(current_cell==middle_cell1){
                    small_index.push_back(i);
                    subsample_vector[0].number_individuals++;
                    subsample_vector[0].species_abundance_vector[species_label]++;
                }

            }
            i_muted_species.clear();
            for(int i=0;i<max_species_label;i++) i_muted_species.push_back(-1);
            for(int i_pre=0;i_pre<small_index.size();i_pre++){
                coalescences_total++;
                i_muted_species[subsample_index_muted[i_pre]]=i_pre;
            }
            for(int ispecies=0;ispecies<max_species_label;ispecies++){
                if(i_muted_species[ispecies]!=-1){
                    mutations_total++;
                    current_mutations++;
                    coalescences_total--;
                }
            }

            ///Here I measure beta-diversity and heterozigozity of the two middle cells
            int all_pairs=0;
            int equal_pairs=0;
            double mean_final_individuals_distance_aux=0;
            double mean_final_individuals_distance_all_aux=0;
            double max_final_individuals_distance_all_aux=0;
            double min_final_individuals_distance_all_aux=10;
            for(int j=0;j<small_index.size();j++){
                int i=small_index[j];
                species_label=subsample_index_muted[j];
                if(print_spatial) (*spatial_data)<<active_individuals_vector[i].x<<" "<<active_individuals_vector[i].y<<" "<<species_label<<endl;
                for(int j2=j+1;j2<small_index.size();j2++){
                    int i2=small_index[j2];
                    species_label2=subsample_index_muted[j2];
                    current_distance_x=fabs(active_individuals_vector[i].x-active_individuals_vector[i2].x);
                    current_distance_y=fabs(active_individuals_vector[i].y-active_individuals_vector[i2].y);
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
                    //ORIGIN (0.25,0.25)
                    current_distance=sqrt(pow(active_individuals_vector[i].x-0.25*Lx,2)+pow(active_individuals_vector[i].y-0.25*Lx,2));
		    //sqrt(pow(active_individuals_vector[i].x-0.25,2)+pow(active_individuals_vector[i].y-0.25,2));
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
            max_final_individuals_distance_all_network=max_final_individuals_distance_all_aux;
            int cummulated_alpha=0;
            for(int i=0;i<number_of_spatial_radius;i++){
                cummulated_alpha+=current_alpha[i];
                mean_alpha[i]+=cummulated_alpha;
            }

            ///Here I count the number of species at each subarea and average over all subsamples
            number_of_occupied_subareas=1;
            number_of_species=0;
            for(int j=0;j<subsample_vector[0].species_abundance_vector.size();j++){
                if(subsample_vector[0].species_abundance_vector[j]!=0){
                    SAD_histogram_subsamples_mean[subsample_vector[0].species_abundance_vector[j]]++;
                    number_of_species++;
                }
            }
            subsample_vector[0].number_species=number_of_species;

            mean_number_individuals_subarea+=subsample_vector[0].number_individuals/double(number_of_occupied_subareas);
            mean_number_species_subarea+=subsample_vector[0].number_species/double(number_of_occupied_subareas);
        }

        species_abundance.clear();
        free(current_alpha);
    }
    void print_abundances(double Lxi,double D_factor){

        sprintf(SAD_name,"SAD_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_rate,NNN,Lx,Lxi,random_seed);
        sprintf(SAD_subsample_name,"SAD_subsample_Df%g_Mut%g_N%d_Lx%g_Lxi%g_Seed%d.dat",D_factor,mutation_rate,NNN,Lx,Lxi,random_seed);
        ///////////////////////// SAD ///////////////////////
        SAD.open(SAD_name, ios::out | ios::trunc);
        for(int i=0;i<MAX_NUMBER_OF_INDIVIDUALS;i++) if(SAD_histogram[i]!=0) SAD<<i<<" "<<SAD_histogram[i]/double(TOTAL_TRIALS)<<" "<<network<<" "<<TOTAL_TRIALS<<endl;
        SAD.close();
        ///////////////////////// SAD ///////////////////////
        SAD_subsample.open(SAD_subsample_name, ios::out | ios::trunc);
        for(int i=0;i<MAX_NUMBER_OF_INDIVIDUALS;i++) if(SAD_histogram_subsamples_mean[i]!=0) SAD_subsample<<i<<" "<<SAD_histogram_subsamples_mean[i]/double(TOTAL_TRIALS*number_of_occupied_subareas)<<" "<<network<<" "<<TOTAL_TRIALS<<endl;
        SAD_subsample.close();
    }
    void clear(){
        active_individuals_vector.clear();
        ancestors_list.clear();
    }
    void print_spatial_distribution(ofstream *file){

        for(int i=0;i<active_individuals_vector.size();i++) {
            (*file)<<active_individuals_vector[i].x<<" "<<active_individuals_vector[i].y<<" "<<t<<endl;
        }
        (*file)<<endl<<endl;
        (*file).close();
    }
};
class fluid2D_tree:public tree{

public:

    ///PARAMETERS
    int kx,ky;
    double Lx,Ly,ratebirth,sqrtdt,death,ndnd,DIFF,Ly_half,y,delta,oneoverd;
    int max_x,max_y;
    tree_individual new_individual;

    ///VARIABLES
    int **bindistr;
    double ratedeath;
    double *jet_parameters;



    fluid2D_tree(int N,double mutprob,int r,double number_of_spatial_radius_in,double Lx_in,double Ly_in,double mu,double DIFF_in,double *extra_parameters=nullptr):tree(N,mutprob,r,number_of_spatial_radius_in,Lx_in,Ly_in){
	cout<<"Starting to construct a tree"<<endl;
        dt=DT;
        Lx=Lx_in;
        Ly=Lx;
        Ly_half=Ly/2.;
        DIFF=DIFF_in;
        oneoverd=(int)floor(sqrt(NNN))/double(Lx);//(int)floor(sqrt(NNN))/double(Lx);
        cout<<"oneoverd is "<<oneoverd<<endl;
	max_x=oneoverd*Lx+1;
        max_y=oneoverd*Ly+1;
	cout<<"max_x is "<<max_x<<endl;
	cout<<"max_y is "<<max_y<<endl;
        bindistr=(int **) malloc(sizeof(int *)*max_x);
        for(int i=0;i<max_x;i++) *(bindistr+i)=(int *)malloc(sizeof(int)*max_y);
        ndnd=(double)NNN/(Lx*Lx*oneoverd*oneoverd);
        death=mu*dt;// /ndnd;
        ratebirth=mu*dt;
	cout<<"DIFF_in = "<<DIFF_in<<endl;
        sqrtdt=sqrt(2.0*dt*DIFF_in);
	cout<<"sqrtrdt = "<<sqrtdt<<endl;
        sprintf(fname,"fixation_%s_Mut%g_N%d_Lx%g_Ly%g_Force%g_Seed%d.dat",model.c_str(),mutation_rate,NNN,Lx,Ly,force,random_seed);
        sprintf(spacetime_name,"spacetime_%s_Mut%g_N%d_Lx%g_Ly%g_Force%g_Seed%d.dat",model.c_str(),mutation_rate,NNN,Lx,Ly,force,random_seed);

        jet_parameters=(double *) malloc(sizeof(double *)*8);
        jet_parameters[0]=extra_parameters[0];
        jet_parameters[1]=extra_parameters[1];
        jet_parameters[2]=extra_parameters[2];
        jet_parameters[3]=extra_parameters[3];
        jet_parameters[4]=extra_parameters[4];
        jet_parameters[5]=extra_parameters[5];
        jet_parameters[6] = extra_parameters[6];//

        //jet_parameters[6]=(Lx*Lx/NNN)/(0.05*0.05/2); //(3*pow(10, -6));// -9));

        //resc_const = (D_factor/(NNN))/pow(10, -9);

        delta=1./oneoverd;///extra_parameters[6];   double interaction_range=double(Lx)/(int)floor(sqrt(NNN));
	cout<<"Interaction range for the tree is "<<delta<<endl;


    }
    void contour_conditions(int i){
        do{
            active_individuals_vector[i].x=fmod((active_individuals_vector[i].x+Lx),Lx);
            active_individuals_vector[i].y=fmod((active_individuals_vector[i].y+Ly),Ly);

        }while((active_individuals_vector[i].x<0)||(active_individuals_vector[i].y<0)||(active_individuals_vector[i].x>Lx)||(active_individuals_vector[i].y>Ly));
    };
    virtual void one_step_dynamics (gsl_rng * random_variable){
	//cout<<"Starting the one step dynamics; size of the loop is "<<active_individuals_vector.size()<<endl;
        int current_actives,index_to_kill;
        vector <int> index_death;

        for(int i=0;i<active_individuals_vector.size();i++){
	    #ifdef ADVECTION
	    double K0x, K0y, K1x, K1y, K2x, K2y, K3x, K3y;
	    //cout<<"Constructing the act indiv vector; i = "<<i<<endl;
            active_individuals_vector[i].advectx(force,t,Ly,jet_parameters,0,0,0);
            active_individuals_vector[i].advecty(force,t,Lx,jet_parameters,0,0,0);
	    //cout<<"calculat the increment"<<endl;
	    K0x = active_individuals_vector[i].fieldx;
	    K0y = active_individuals_vector[i].fieldy;

            active_individuals_vector[i].advectx(force,t,Ly,jet_parameters,K0x,K0y,dt/2);
            active_individuals_vector[i].advecty(force,t,Lx,jet_parameters,K0x,K0y,dt/2);
            K1x = active_individuals_vector[i].fieldx;
            K1y = active_individuals_vector[i].fieldy;

            active_individuals_vector[i].advectx(force,t,Ly,jet_parameters,K1x,K1y,dt/2);
            active_individuals_vector[i].advecty(force,t,Lx,jet_parameters,K1x,K1y,dt/2);
            K2x = active_individuals_vector[i].fieldx;
            K2y = active_individuals_vector[i].fieldy;

            active_individuals_vector[i].advectx(force,t,Ly,jet_parameters,K2x,K2y,dt);
            active_individuals_vector[i].advecty(force,t,Lx,jet_parameters,K2x,K2y,dt);
            K3x = active_individuals_vector[i].fieldx;
            K3y = active_individuals_vector[i].fieldy;

            active_individuals_vector[i].x+=sqrtdt*gsl_ran_gaussian(random_variable,1)+dt*(K0x+2*K1x+2*K2x+K3x)/6;
            active_individuals_vector[i].y+=sqrtdt*gsl_ran_gaussian(random_variable,1)+dt*(K0y+2*K1y+2*K2y+K3y)/6;
	    #else
            active_individuals_vector[i].x+=sqrtdt*gsl_ran_gaussian(random_variable,1);
            active_individuals_vector[i].y+=sqrtdt*gsl_ran_gaussian(random_variable,1);
	    #endif
            contour_conditions(i);
	    //cout<<"x increment is "<<sqrtdt*gsl_ran_gaussian(random_variable,1)+dt*(1.5*active_individuals_vector[i].fieldx-0.5*active_individuals_vector[i].fieldoldx)<<endl;
	    //cout<<"y increment is "<<sqrtdt*gsl_ran_gaussian(random_variable,1)+dt*(1.5*active_individuals_vector[i].fieldy-0.5*active_individuals_vector[i].fieldoldy)<<endl;
	    //cout<<"starting to check the contour condit"<<endl;
	    //cout<<"(x,y) = "<<active_individuals_vector[i].x<<", "<<active_individuals_vector[i].y<<endl;
	    //cout<<"act individ vect .x "<<active_individuals_vector[i].x<<endl;
	    //cout<<"act individ vect .y"<<active_individuals_vector[i].y<<endl;
	    //contour_conditions(i);
	    //cout<<active_individuals_vector[i].x<<", "<<active_individuals_vector[i].y<<endl;

	    //cout<<"finished contour condit"<<endl;
            active_individuals_vector[i].fieldoldx=active_individuals_vector[i].fieldx;
            active_individuals_vector[i].fieldoldy=active_individuals_vector[i].fieldy;
	    //cout<<"Finish th eloop step"<<endl;
        }


        for(kx=0;kx<max_x;kx++) for(ky=0;ky<max_y;ky++) bindistr[kx][ky]=0;
        for(int i=0;i<active_individuals_vector.size();i++) bindistr[(int)floor(active_individuals_vector[i].x*oneoverd)][(int)floor(active_individuals_vector[i].y*oneoverd)]++;
	//for(kx=0;kx<max_x;kx++) for(ky=0;ky<max_y;ky++) cout<<"the matr is "<<bindistr[kx][ky]<<" "<<endl;
	//for(kx=0;kx<max_x;kx++) cout<<"the matr is "<<bindistr[kx]<<" "<<endl;

        ///DEATH PROCESS
        current_actives=active_individuals_vector.size();
        for(int i=0;i<current_actives;i++){
            kx=(int)floor(active_individuals_vector[i].x*oneoverd);
            ky=(int)floor(active_individuals_vector[i].y*oneoverd);
	    //cout<<"x = "<<active_individuals_vector[i].x<<"; y = "<<active_individuals_vector[i].y<<endl;
            //cout<<"kx = "<<kx<<"; ky = "<<ky<<"; of "<<max_x<<endl;
	    ratedeath=death*(bindistr[kx][ky]-1);
	    //cout<<"death = "<<death<<"; bindistr = "<<bindistr[kx][ky]<<endl;
            if(gsl_rng_uniform(random_variable)<ratedeath) index_death.push_back(i);
        }

        for(int i=index_death.size()-1;i>=0;i--){
            index_to_kill=index_death[i];
            remove_tree_individual(&active_individuals_vector[index_to_kill]);
            active_individuals_vector[index_to_kill]=active_individuals_vector.back();
            active_individuals_vector.pop_back();
        }
	//cout<<index_death.size()<<" individ. died; death rate is "<<ratedeath<<endl;

        current_actives=active_individuals_vector.size();
        ///REPRODUCTION PROCESS
	int count_repr = 0;
        for(int i=0;i<current_actives;i++){
            if(gsl_rng_uniform(random_variable)<ratebirth){
		count_repr = count_repr+1;
		//cout<<count_repr<<" indiv were born"<<endl;
                active_individuals_vector.push_back(active_individuals_vector[i]);
                active_individuals_vector.back().x+=delta*(gsl_rng_uniform(random_variable)-1/2.);
                active_individuals_vector.back().y+=delta*(gsl_rng_uniform(random_variable)-1/2.);
                add_tree_individual(&active_individuals_vector.back(),&active_individuals_vector[i]);
            }
	}
	//cout<<count_repr<<" individ were born; birth rate is"<<ratebirth<<endl;
    }
};






#endif //TREE_CLASSES_H
