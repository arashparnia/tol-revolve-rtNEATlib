#include <iostream>
//#include <iomanip>
#include "String"
#include "neat.h"
#include "nnode.h"
#include "trait.h"
#include "link.h"
#include "population.h"
//#include "experiments.h"

//Perform evolution on spline, for gens generations
//NEAT::Population *spline_test(int gens);

using namespace std;

using namespace NEAT;



int spline_epoch(Population *pop, int generation, char *filename, int &winnernum, int &winnergenes, int &winnernodes);

bool spline_evaluate(Organism *org);

int spline_realtime_loop(NEAT::Population *pop);





int main(int argc, char *argv[]) {


    //***********RANDOM SETUP***************//
    /* Seed the random-number generator with current time so that
        the numbers will be different every time we run.    */
    srand((unsigned) time(NULL));

    //Load in the params
    //load_neat_params("./spline_param.ne",true);

    //cout<<"loaded parameters from ./spline_param.ne "<<endl;

    NEAT::time_alive_minimum = 10;

    NEAT::trait_param_mut_prob = 0;
    NEAT::trait_mutation_power = 0; // Power of mutation on a signle trait param
    NEAT::linktrait_mut_sig = 0; // Amount that mutation_num changes for a trait change inside a link
    NEAT::nodetrait_mut_sig = 0; // Amount a mutation_num changes on a link connecting a node that changed its trait

    NEAT::weight_mut_power = 0.01; // The power of a linkweight mutation
    NEAT::recur_prob = 0; // Prob. that a link mutation which doesn't have to be recurrent will be made recurrent
    NEAT::disjoint_coeff = 0;
    NEAT::excess_coeff = 0;
    NEAT::mutdiff_coeff = 0;
    NEAT::compat_threshold = 0; //Modify compat thresh to control speciation
    NEAT::age_significance = 0; // How much does age matter?
    NEAT::survival_thresh = 0.1; // Percent of ave fitness for survival

    NEAT::mutate_only_prob = 0; // Prob. of a non-mating reproduction
    NEAT::mutate_random_trait_prob = 0;
    NEAT::mutate_link_trait_prob = 0;
    NEAT::mutate_node_trait_prob = 0;
    NEAT::mutate_link_weights_prob = 0.05;
    NEAT::mutate_toggle_enable_prob = 0.0;
    NEAT::mutate_gene_reenable_prob = 0.0;
    NEAT::mutate_add_node_prob = 0;
    NEAT::mutate_add_link_prob = 0;
    NEAT::interspecies_mate_rate = 0; // Prob. of a mate being outside species
    NEAT::mate_multipoint_prob = 1;
    NEAT::mate_multipoint_avg_prob = 1;
    NEAT::mate_singlepoint_prob = 1;
    NEAT::mate_only_prob = 0; // Prob. of mating without mutation
    NEAT::recur_only_prob = 0;  // Probability of forcing selection of ONLY links that are naturally recurrent
    NEAT::pop_size = 50;  // Size of population
    NEAT::dropoff_age = 10;  // Age where Species starts to be penalized
    NEAT::newlink_tries = 0;  // Number of tries mutate_add_link will attempt to find an open link
    NEAT::print_every = 1; // Tells to print population to file every n generations
    NEAT::babies_stolen = 0; // The number of babies to siphen off to the champions
    NEAT::num_runs = 10;


    Population *pop=0;
    Genome *start_genome;
    start_genome = new Genome(1, 20, 0, 0);
    cout << "Spawning Population off Genome" << endl;

    pop = new Population(start_genome, NEAT::pop_size);

//    for (NEAT::Organism *o : pop->organisms) {
//
//        cout << "CHECK " << o->net->numlinks << endl;
//
//    }



//    pop = spline_test(100);  // generation  experiment



    cout << " population size " << pop->organisms.size() << endl;
    cout << "REAL TIME" << endl;

    spline_realtime_loop(pop);



//    if (pop)
//        delete pop;
    return (0);
}



int spline_realtime_loop(Population *pop) {



    vector<Organism *>::iterator curorg;
    vector<Species *>::iterator curspecies;

    vector<Species *>::iterator curspec; //used in printing out debug info

    vector<Species *> sorted_species;  //Species sorted by max fit org in Species

    int pause;
    bool win = false;

    double champ_fitness;
    Organism *champ;



    //Real-time evolution variables
    int offspring_count;
    Organism *new_org;


    //We try to keep the number of species constant at this number
    int num_species_target = NEAT::pop_size / 5;

    //This is where we determine the frequency of compatibility threshold adjustment
    int compat_adjust_frequency = NEAT::pop_size / 10;
    if (compat_adjust_frequency < 1)
        compat_adjust_frequency = 1;

    //Initially, we evaluate the whole population
    //Evaluate each organism on a test
    for (curorg = (pop->organisms).begin(); curorg != (pop->organisms).end(); ++curorg) {

        //shouldn't happen
        if (((*curorg)->gnome) == 0) {
            cout << "ERROR EMPTY GEMOME!" << endl;
            cin >> pause;
        }

        if (spline_evaluate((*curorg))) win = true;

    }

    //Get ready for real-time loop

    //Rank all the organisms from best to worst in each species
    pop->rank_within_species();

    //Assign each species an average fitness
    //This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
    pop->estimate_all_averages();


    //Now create offspring one at a time, testing each offspring,
    // and replacing the worst with the new offspring if its better
    for (offspring_count = 0; offspring_count < 20000; offspring_count++) {


        //Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
        //and reassign the population to new species
        if (offspring_count % compat_adjust_frequency == 0) {

            int num_species = pop->species.size();
            double compat_mod = 0.1;  //Modify compat thresh to control speciation

            // This tinkers with the compatibility threshold
            if (num_species < num_species_target) {
                NEAT::compat_threshold -= compat_mod;
            } else if (num_species > num_species_target)
                NEAT::compat_threshold += compat_mod;

            if (NEAT::compat_threshold < 0.3)
                NEAT::compat_threshold = 0.3;

            cout << "compat_thresh = " << NEAT::compat_threshold << endl;

            //Go through entire population, reassigning organisms to new species
            for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
                pop->reassign_species(*curorg);
            }
        }


        //For printing only
        for (curspec = (pop->species).begin(); curspec != (pop->species).end(); curspec++) {
            cout << "Species " << (*curspec)->id << " size" << (*curspec)->organisms.size() << " average= "
                 << (*curspec)->average_est << endl;
        }

        cout << "Pop size: " << pop->organisms.size() << endl;

        //Here we call two rtNEAT calls:
        //1) choose_parent_species() decides which species should produce the next offspring
        //2) reproduce_one(...) creates a single offspring fromt the chosen species
        new_org = (pop->choose_parent_species())->reproduce_one(offspring_count, pop, pop->species);

        //Now we evaluate the new individual
        //Note that in a true real-time simulation, evaluation would be happening to all individuals at all times.
        //That is, this call would not appear here in a true online simulation.
        cout << "Evaluating new baby: " << endl;
        if (spline_evaluate(new_org)) win = true;

        if (win) {
            cout << "WINNER" << endl;
            pop->print_to_file_by_species("spline_rt_winpop");

            cout << endl << "WINNER GENE " << endl;
            for (NEAT::Gene *g : new_org->gnome->genes) {
                cout << " innovation " << g->innovation_num;
                cout << " link " << g->lnk->weight;
                cout << endl;
            }

            //    for (NEAT::Trait *t : org->gnome->traits) {
            //        for (int i = 0; i < 10; i++) {
            //            cout << " traits " << t->params[i];
            //            cout << endl;
            //        }
            //    }

            break;
        }

        //Now we reestimate the baby's species' fitness
        new_org->species->estimate_average();

        //Remove the worst organism
        pop->remove_worst();

    }

    return 0;
}



int spline_epoch(Population *pop, int generation, char *filename, int &winnernum, int &winnergenes, int &winnernodes) {
vector<Organism *>::iterator curorg;
vector<Species *>::iterator curspecies;

bool win = false;

//Evaluate each organism on a test
for (curorg = (pop->organisms).begin(); curorg != (pop->organisms).end(); ++curorg) {
if (spline_evaluate(*curorg)) {
win = true;
winnernum = (*curorg)->gnome->genome_id;
winnergenes = (*curorg)->gnome->extrons();
winnernodes = ((*curorg)->gnome->nodes).size();
if (winnernodes == 5) {
//You could dump out optimal genomes here if desired
//(*curorg)->gnome->print_to_filename("xor_optimal");
cout<<"OPTIMAL GENE"<<endl;
cout<<(*curorg)->gnome->genome_id<<endl;
}
}
}

//Average and max their fitnesses for dumping to file and snapshot
for (curspecies = (pop->species).begin(); curspecies != (pop->species).end(); ++curspecies) {

//This experiment control routine issues commands to collect ave
//and max fitness, as opposed to having the snapshot do it,
//because this allows flexibility in terms of what time
//to observe fitnesses at

(*curspecies)->compute_average_fitness();
(*curspecies)->compute_max_fitness();
}

//Take a snapshot of the population, so that it can be
//visualized later on
//if ((generation%1)==0)
//  pop->snapshot();

//Only print to file every print_every generations
if (win ||
((generation % (NEAT::print_every)) == 0))
pop->print_to_file_by_species(filename);


if (win) {
for (curorg = (pop->organisms).begin(); curorg != (pop->organisms).end(); ++curorg) {
if ((*curorg)->winner) {
cout << "WINNER IS #" << ((*curorg)->gnome)->genome_id << endl;
//Prints the winner to file
//IMPORTANT: This causes generational file output!
print_Genome_tofile((*curorg)->gnome, "winner");
}
}

}

pop->epoch(generation);

if (win) return 1;
else return 0;

}

bool spline_evaluate(Organism *org) {
    cout << endl << "evaluating organism from generation " << org->generation << endl;

    // PRINTING STUFF
    Network *net;
    double in[1] = {1};
    double out[20] ;
    double this_out;
    bool success = false;
    int net_depth; //The max depth of the network to be activated
    int relax; //Activates until relaxation

    net=org->net;
    net->flush();
    net->load_sensors(in);
    cout << "network has " << net->numnodes << " nodes and " << net->numlinks <<" links" << endl;
    //Relax net and get output
    success=net->activate();
    net->show_activation();
    for(int i=0;i<20;i++) {
        out[i] = net->outputs[i]->output;
//        out[count]=(*(net->outputs.begin()))->activation;
    }
//    cout << "neural net output: " << endl;
//    for (int i = 0 ; i<20;i++){
//        cout  << out[i] << " , ";
//    }
//    cout << endl;


    double sum;

    for (NEAT::Gene *g : org->gnome->genes) {
//        cout <<  g->innovation_num;
        cout << g->lnk->weight << " , ";

        sum +=  g->lnk->weight;
    }
    cout << endl;

    //----------------------------------------------------------------------------------------------------

    float fittness;
//    do {
//        fittness = rand();//static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//    }while (fittness < 0.9 && fittness > 0.1);

//    float fittness;
//    cout << endl <<"FITTNESS: ";
//    cin >> fittness;
    fittness = sum;


    double errorsum;



    if (fittness >= survival_thresh) success = true; else success = false;

    if (success) {
        errorsum = (1 - fittness);
        org->fitness = fittness;
        org->error = errorsum;
    } else {
        //The network is flawed (shouldnt happen)
        errorsum = 999.0;
        org->fitness = 0.001;
    }


    cout << " the fittness is" << fittness << endl;

    if (errorsum < 1-survival_thresh) {
        org->winner = true;
        return true;
    } else {
        org->winner = false;
        return false;
    }



}




//single test
//Perform evolution on , for gens generations
NEAT::Population *spline_test(int gens) {
    Population *pop = 0;

    //memset(pop,0,sizeof(pop));

    Genome *start_genome;
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;

    int evals[NEAT::num_runs];  //Hold records for each run
    int genes[NEAT::num_runs];
    int nodes[NEAT::num_runs];
    int winnernum;
    int winnergenes;
    int winnernodes;
    //For averaging
    int totalevals = 0;
    int totalgenes = 0;
    int totalnodes = 0;
    int expcount;
    int samples;  //For averaging

    memset(evals, 0, NEAT::num_runs * sizeof(int));
    memset(genes, 0, NEAT::num_runs * sizeof(int));
    memset(nodes, 0, NEAT::num_runs * sizeof(int));



    /*
     * reading genome1 from file
    ifstream iFile("xorstartgenes",ios::in);

    cout<<"START XOR TEST"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    //start_genome = new Genome(int new_id,int i, int o, int n,int nmax, bool r, double linkprob);

    iFile.close();
    */


//    //manual deceleration of genome
//
    std::vector<Trait*> traits; //parameter conglomerations
    std::vector<NNode*> nnodes; //List of NNodes for the Network
    std::vector<Link*> links; //List of innovation-tracking genes
//
//
//
    NEAT::Trait *trait1 = new Trait(0,0,0,0,0,0,0,0,0,0);
    traits.push_back(trait1);
//
//    /*
//     * enum nodetype {
//		NEURON = 0,
//		SENSOR = 1
//	};
//
//	enum nodeplace {
//		HIDDEN = 0,
//		INPUT = 1,
//		OUTPUT = 2,
//		BIAS = 3
//	};
//
//	enum functype {
//		SIGMOID = 0
//	};
//     // A NODE is either a NEURON or a SENSOR.
//	//   - If it's a sensor, it can be loaded with a value for output
//	//   - If it's a neuron, it has a list of its incoming input signals (List<Link> is used)
//     *NNode(nodetype ntype,int nodeid, nodeplace placement);
//     */
//
//    //input nodes
    NEAT::NNode *nnodeIn1;
    nnodeIn1 = new NEAT::NNode(NEURON, 1,INPUT);
//
//    //hidden nodes
//
//    //output nodes
    NNode *nnodeOut1;
    nnodeOut1 = new NEAT::NNode(SENSOR, 1,OUTPUT);


    //links

    Link *link1 = new Link(0,nnodeIn1,nnodeOut1,false);
    links.push_back(link1);

    //Genome(int id, std::vector<Trait*> t, std::vector<NNode*> n, std::vector<Link*> links);

//    start_genome = new Genome(1,traits,nnodes,links);


    //Special constructor that creates a Genome of 3 possible types:
    //0 - Fully linked, no hidden nodes
    //1 - Fully linked, one hidden node splitting each link
    //2 - Fully connected with a hidden layer, recurrent
    //num_hidden is only used in type 2
    //Genome(int num_in,int num_out,int num_hidden,int type);

    start_genome = new Genome(1, 20, 0, 0);
//    Network *nn =start_genome->phenotype;

//    cout << nn->numlinks <<endl;
//
//    for (int i = 0; i < 20; ++i) {
//        start_genome->genes[i]->lnk->weight = 0;
//
//    }
//    start_genome->nodes[5]->activation = 0;


    cout << "starrting gene is " <<endl;
    for (NEAT::Gene *g : start_genome->genes) {
        cout << " innovation " << g->innovation_num;
        cout << " link " << g->lnk->weight;
        cout << endl;
    }




    for (expcount = 0; expcount < NEAT::num_runs; expcount++) {
        //Spawn the Population
        cout << "Spawning Population off Genome" << endl;

        pop = new Population(start_genome, NEAT::pop_size);

        for (NEAT::Organism *o : pop->organisms) {

            cout << "CHECK " << o->net->numlinks << endl;

        }


        cout << "Verifying Spawned Pop" << endl;
        pop->verify();

        for (gen = 1; gen <= gens; gen++) {
            cout << "Epoch " << gen << endl;

            //This is how to make a custom filename
//            fnamebuf = new ostringstream();
//            (*fnamebuf) << "gen_" << gen << ends;  //needs end marker

#ifndef NO_SCREEN_OUT
//            cout << "name of fname: " << fnamebuf->str() << endl;
#endif


            char temp[50];
            sprintf(temp, "gen_%d", gen);

            //Check for success
            if (spline_epoch(pop, gen, temp, winnernum, winnergenes, winnernodes)) {
                //Collect Stats on end of experiment
                evals[expcount] = NEAT::pop_size * (gen - 1) + winnernum;
                genes[expcount] = winnergenes;
                nodes[expcount] = winnernodes;
                gen = gens;
            }
            //Clear output filename
//            fnamebuf->clear();
//            delete fnamebuf;
        }

        if (expcount < NEAT::num_runs - 1) delete pop;

    }


    //Average and print stats
    cout << "Nodes: " << endl;
    for (expcount = 0; expcount < NEAT::num_runs; expcount++) {
        cout << nodes[expcount] << endl;
        totalnodes += nodes[expcount];
    }

    cout << "Genes: " << endl;
    for (expcount = 0; expcount < NEAT::num_runs; expcount++) {
        cout << genes[expcount] << endl;
        totalgenes += genes[expcount];
    }

    cout << "Evals " << endl;
    samples = 0;
    for (expcount = 0; expcount < NEAT::num_runs; expcount++) {
        cout << evals[expcount] << endl;
        if (evals[expcount] > 0) {
            totalevals += evals[expcount];
            samples++;
        }
    }

    cout << "Failures: " << (NEAT::num_runs - samples) << " out of " << NEAT::num_runs << " runs" << endl;
    cout << "Average Nodes: " << (samples > 0 ? (double) totalnodes / samples : 0) << endl;
    cout << "Average Genes: " << (samples > 0 ? (double) totalgenes / samples : 0) << endl;
    cout << "Average Evals: " << (samples > 0 ? (double) totalevals / samples : 0) << endl;

    return pop;
}

