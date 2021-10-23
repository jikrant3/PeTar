#pragma once

//! IO parameters for Galpy manager
class IOParamsDisk{
public:
    IOParamsContainer input_par_store;
    IOParams<double> gravitational_constant;
    IOParams<double> semi_in;
    IOParams<double> semi_out;
    IOParams<double> epsilon;
    IOParams<double> r0;
    IOParams<double> tau_dep;
    IOParams<double> zk;
    IOParams<double> imc_mass;
    IOParams<long long int> power_order;
    IOParams<long long int> center_id;

    bool print_flag;
    
    IOParamsDisk(): input_par_store(),
                    gravitational_constant(input_par_store, 1.0, "disk-G","gravitational constant"),
                    semi_in (input_par_store, 1e-5, "disk-ain", "inner semi-major axis"),
                    semi_out(input_par_store, 1e-2, "disk-ain", "inner semi-major axis"),
                    epsilon (input_par_store, 1.0, "disk-epsilon", "epsilon 0"),
                    r0      (input_par_store, 1e-2, "disk-R0", "R0"),
                    tau_dep (input_par_store, 1.0, "disk-tau-dep", "tau_dep"),
                    zk(input_par_store, 1.094, "disk-zk", "Z_k"),
                    imc_mass(input_par_store, 100, "disk-imc-mass", "mass cutoff to find IMC"),
                    power_order(input_par_store, 20, "disk-order", "power index order of force"),
                    center_id(input_par_store, 1, "disk-center-id", "ID of centeral particle (SMBH)"),
                    print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int disk_flag=-1;
        const struct option long_options[] = {
            {gravitational_constant.key, required_argument, &disk_flag, 0}, 
            {semi_in.key,  required_argument, &disk_flag, 1}, 
            {semi_out.key, required_argument, &disk_flag, 2}, 
            {epsilon.key,  required_argument, &disk_flag, 3}, 
            {r0.key,       required_argument, &disk_flag, 4}, 
            {tau_dep.key,  required_argument, &disk_flag, 5}, 
            {power_order.key,  required_argument, &disk_flag, 6}, 
            {center_id.key,    required_argument, &disk_flag, 7}, 
            {imc_mass.key,     required_argument, &disk_flag, 8}, 
            {zk.key,           required_argument, &disk_flag, 9}, 
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;

        while ((copt = getopt_long(argc, argv, "-p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (disk_flag) {
                case 0:
                    gravitational_constant.value = atof(optarg);
                    if (print_flag) gravitational_constant.print(std::cout);
                    opt_used+= 2;
                    break;
                case 1:
                    semi_in.value = atof(optarg);
                    if (print_flag) semi_in.print(std::cout);
                    opt_used+= 2;
                    break;
                case 2:
                    semi_out.value = atof(optarg);
                    if (print_flag) semi_out.print(std::cout);
                    opt_used+= 2;
                    break;
                case 3:
                    epsilon.value = atof(optarg);
                    if (print_flag) epsilon.print(std::cout);
                    opt_used+= 2;
                    break;
                case 4:
                    r0.value = atof(optarg);
                    if (print_flag) r0.print(std::cout);
                    opt_used+= 2;
                    break;
                case 5:
                    tau_dep.value = atof(optarg);
                    if (print_flag) tau_dep.print(std::cout);
                    opt_used+= 2;
                    break;
                case 6:
                    power_order.value = atoi(optarg);
                    if (print_flag) power_order.print(std::cout);
                    opt_used+= 2;
                    break;
                case 7:
                    center_id.value = atoi(optarg);
                    if (print_flag) center_id.print(std::cout);
                    opt_used+= 2;
                    break;
                case 8:
                    imc_mass.value = atof(optarg);
                    if (print_flag) imc_mass.print(std::cout);
                    opt_used+= 2;
                    break;
                case 9:
                    zk.value = atof(optarg);
                    if (print_flag) zk.print(std::cout);
                    opt_used+= 2;
                    break;
                default:
                    break;
                }
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par+".disk"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fgalpy_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fgalpy_par.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used+=2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"Disk options:"<<std::endl;
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of Galpy -----\n";

        return opt_used;
    }    
            
};

//! Disk pontetial
class DiskPotentialManager{
public:
    double time;
    double gravitational_constant;
    double semi_in;
    double semi_out;
    double epsilon0;
    double r0;
    double tau_dep;
    double imc_mass;
    double zk;
    int power_order;
    int center_id;

    double epsilon_prefix;

    //! initialization function
    /*!
      @param[in] _input: input parameters 
      @param[in] _time: current system time
      @param[in] _print_flag: if true, printing information to std::cout
     */
    void initial(const IOParamsDisk& _input, const double _time, const bool _print_flag=false) {
        // 
        time = _time;
        gravitational_constant = _input.gravitational_constant.value;
        semi_in = _input.semi_in.value;
        semi_out = _input.semi_out.value;
        epsilon0 = _input.epsilon.value;
        r0 = _input.r0.value;
        tau_dep = _input.tau_dep.value;
        power_order = _input.power_order.value;
        center_id = _input.center_id.value;
        imc_mass = _input.imc_mass.value;
        zk = _input.zk.value;

        if (_print_flag) {
            std::cout<<"Disk potential, time: "<<time
                     <<" G: "<<gravitational_constant
                     <<" a_in: "<<semi_in
                     <<" a_out: "<<semi_out
                     <<" epsilon_0: "<<epsilon0
                     <<" R0: "<<r0
                     <<" tau_dep: "<<tau_dep
                     <<" power_order: "<<power_order
                     <<" c.m. id: "<<center_id
                     <<" IMC_mass: "<<imc_mass
                     <<" Zk: "<<zk
                     <<std::endl;
        }

        updateParameters(time);
    }

    //! update epsilon coefficient at given time
    /*! 4*pi*G*E0*R0^(3/2)*exp(-t/tau_dep)
      @param[in] time: time in input unit
    */
    void updateParameters(const double time) {
        const double pi = 3.141592653589793;
        epsilon_prefix = 4*pi*gravitational_constant*epsilon0 * std::pow(r0, 1.5) * std::exp(-time/tau_dep);
    }

    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] mass: mass of particle
      @param[in] pos: position of particles referring to SMBH position [input unit];
     */
    void calcAccPot(double* acc, double& pot, double& mass, const double* pos) {
        double r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        double coefficient = epsilon_prefix*std::pow(r, 1.5);
        
        if (mass>imc_mass) { //IMC
            double a = -coefficient*zk;
            acc[0] = a*pos[0]/r;
            acc[1] = a*pos[1]/r;
            acc[2] = a*pos[2]/r;

            pot = 2*coefficient*r*zk;
        }
        else { //Star
            double i_fact = 1;    // i!
            double i2_fact = 1;   // (2*i)!
            double two_power = 1;  // 2^(2*i)

            double a = std::pow(semi_in/r, 0.5); // i=0 no coeff
            pot = std::pow(r/semi_out, 0.5) + std::pow(semi_in/r, 0.5); // i=0 no coeff

            for (int i=1; i<=power_order; i++) {
                i_fact *= i_fact*i;
                i2_fact *= i2_fact*i*2;
                two_power *= 4;

                double Ai = i2_fact/(two_power * i_fact*i_fact);
                Ai = Ai*Ai;
                double ra_left = std::pow(r/semi_out, 2*i+0.5);
                double ra_right = std::pow(semi_in/r, 2*i+0.5);
                a += Ai/(4*i+1) * (2*i* ra_left - (2*i+1)*ra_right);
                pot += Ai/(4*i+1) * (ra_left + ra_right);
            }
            
            a *= coefficient;
            pot *= coefficient;
            
            acc[0] = a*pos[0]/r;
            acc[1] = a*pos[1]/r;
            acc[2] = a*pos[2]/r;

        }

    }
};
