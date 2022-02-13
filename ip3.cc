#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <random>

//Equation for gate variables:
inline double x_inf(double v,double V12,double k) { return 1./(1.+exp(-(v-V12)/k)); }
inline double tau_inf(double tau,double v,double tauV12,double tauk) { return tau/(cosh((v-tauV12)/tauk)); }
inline double Gr(double v) { return 4.34e-5*exp(-0.0211539274*v); }//paper has error 4.34e5 should be e-5
inline double Gd1(double v) { return 0.075+0.043*tanh((v+20)-20); }
inline double Gv(double v) { return (10.6408- 14.6408*exp(-v/42.7671))/v; }


//NaF
const double mV12 = -43.8, mk = 6;
const double mtauV12 = -43.8, mtauk = 14;
const double hV12 = -67.5, hk = -11.8;
const double htau = 8.46, htauk = 12.8;
double gNaF = 150, htauV12 = -67.5, mtau = 0.25;

//NaP
const double mpV12 = -47.1, mpk = 3.1; 
const double mptau = 1, mptauV12 = -47.1, mptauk = 6.2;
double hpV12 = -60, hpk = -9;
const double hptau = 5000, hptauV12 = -60;
double GNAPBLK = 1, hptauk = 9;

//CaL
const double mcV12 = -27.5, mck = 5.7;
const double mctau = 0.5;
const double hcV12 = -52.4, hck = -5.2;
const double hctau = 18;

double Caout=4, alphaCa=2.5e-5,tauCa=500, Cain0=0.0000000001, Cav=0; 
double gCaL=0;
double Kpump=1e-3,Vpump=Kpump/tauCa;

//////////////////////////////////////////////////
//IP3 variables
//////////////////////////////////////////////////
//turn ER ca dynamics on (1) or off (0)
double ER_on_off = 1;

//m_ER (activation) parameters for J_ER_in
double IP3=0.0015,ki=0.001,kaa=0.0001;

// h_ER (inactivation) parameters for of J_ER_in
double lp=0.1,kd=0.0002;

//Current to concentration conversion factors (fca and Aca are the same and already declared above with alphaCa)
double fca=alphaCa,Aca=alphaCa,sigma_ca=0.185;

//Pump strength (Already declared above at 50ms)
//double tauca=500;

//SERCA pump variables: strength (gSERCA mM/ms) and activation (kSERCA mM)
double gSERCA=0.45,kSERCA=0.00005;

//J_ER_in leak conductance (gL_ER) and Maximal conductance (gER)
double gL_ER=0.1,gER=31000*2.5;

///////////////////////////////////////////////////
///////////////////////////////////////////////////

//////////////////////////////////////////////////
//Kout diffustion and glia variables
//////////////////////////////////////////////////
//#Diffusion
double taudiff=750;

//#Glia #Gmax is in units of mM/mS
double zg=10,kf=6.25,Gmax=0.5;

//#intra to extracellular ratio?
double alpha_kout=0.000105;
///////////////////////////////////////////////////
///////////////////////////////////////////////////


//Kdr
const double   nAk = 5;
const double nB = 0.17, nBk = 40;
double nA = 0.011, nAV = 44, nBV = 49;
double gKdr = 220; 
//CAN
double K_CAN = 0.74*1e-3, nc = 0.97, sigma_CAN=-0.05e-3;
double ECAN=0;
const double hcanC12 = .00015, hcank = -.0002;
const double hcantau = 100;

//ChR2
double Echr2 = 0;
double gma_chr2=0.1, eps1=0.8535, eps2=0.14;
double sig_ret=10e-20, Gd2=0.05, tau_chr2=1.3;

//Arch
double Earch = -145;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Synaptic Depression parameters: Each neuron has it's own D_syn value which represents the current level of depression for that neurons synapes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double f_D = 0.2, D_0 = 1, tau_D = 1000; //f_D is the relative amount that the synapse is depressed with each spike, D_0 is the starting value for depression, and tau_D is the recovery rate from depression
					      //parameters are estimated from Kottick and Del Negro 2015
double pD = 0.0, pDr = 1-pD;//(not used right now)pD is the proportion not able to depress
double D_on_off = 1; //Turns depression on (1) or off (0)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double c=36;
double NPM = 100;//number of pacemakers or Sst-

//concentrations
double Nain=15, Naout = 120;
double Kin = 125, kbath=4;                    

double ENa=26.54*log(Naout/Nain);

const double Esyn=-10;


double rnd() {return double(rand())/RAND_MAX; }



double tausyn=5;

double IAPP = 0;
double stim_ON=0, stim_W=0, cells[2]={75,85};

//Ca2+ squarewave 
double square_ON_OFF = 0;//Turning the ca2+ square wave transient on/off: 0=off 1=on
double Iapp_ca=0;

//tonic inhib
double ginh = 0;
double Einh=-70;

//OPIOIDS
double gmu=0;//conductance of mu-opioid activated potassium channel
double musyndep=1;//depression of synaptic strength by opioids

//Values for leak dist. exp((kbath-A_gl)/B_gl);
double A_gl = 3.425;
double B_gl = 4.05;




//////////////////////////////////////
//synaptic weights
//////////////////////////////////////
//nonPM to nonPM(2->2)
double Wnpm2npm = 0.0063;//Weight of nonPM to nonPM
double Pnpm2npm = 0.02;//Prob of nonPM to nonPM
//PM to nonPM (1->2)
double Wpm2npm = 0.000175;//Weight of PM to nonPM
double Ppm2npm = 0.3;//0.25;//Prob of PM to nonPM
//nonPM to PM(2->1)
double Wnpm2pm = 0.05*5;//Weight of PM to nonPM
//double Pnpm2pm = prob;//This is defined below because the parameter prob is passed through the code.
//////////////////////////////////////
//synaptic weights
//////////////////////////////////////

//////////////////////////////////////
//gL, gNaP dist. parameters
//////////////////////////////////////
double rho=0.8;//covariance
double sL_pct=0.05;//gleak dist pct of mean
//////////////////////////////////////
//gL, gNaP dist. parameters
//////////////////////////////////////


class Neuron
{
	public:

// parameters	

	double gCa,gNaP,gCAN,gChR2,gArch;
	double g_L, ID,stimIDs;

// dynamical variables

	double v;
	double m,h; // INaF
	double mp,hp; // INaP
	double mc,hc; // ICa
	double n; // IKdr
	double Cain;
	double Cav;
	double gsyn; 
	double Nain;
	double hcan;
	double OP1,OP2,CL1, CL2, Pchr2; //IChR2
	double Carch,Oarch,Iarch;//arch channel states
	double m_ER,h_ER,J_ER_in, J_ER_out, caER, catot;
	double EK,Kout;
	double D_syn;
	double E_leak;

	Neuron();
	void init(double,double,double,double,double,double);
	void step(double dt,double drive,double fcan, double fsynCa, double irr);
	int spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr);
};

struct connection
{
	int source,target;
	double weight;
};

class Population
{
	public:

	Neuron* net;
	int size;

	connection *w;
	int sex;

	double vth;
	
	Population(int n,double gca0,double gca1,double gl0,double gl1,double gnap0,double gnap1,
		double gcan0,double gcan1,double gchr20,double gchr21,double garch0,double garch1,double prob,double wght, double stimnum)
	{
		vth=-35;
		size=n;
		net=new Neuron [size];
		std::default_random_engine generator;//random number generator for gnap/gleak dist. needs to be before loop
		for(int i=0;i<size;i++)
		{	
			
                     	net[i].ID=i;
                     	net[i].stimIDs=0;
                     	net[i].catot=0.001;
                     	net[i].Cain=0.0000001;
			net[i].gCa=gca0+(gca1-gca0)*rnd();
			net[i].gCAN=gcan0+(gcan1-gcan0)*rnd();
			net[i].gNaP=gnap0+(gnap1-gnap0)*rnd();
			net[i].g_L=gl0+(gl1-gl0)*rnd();
			net[i].gChR2=0.0;//rnd();
			net[i].gArch=0.0;//rnd();//Distrabution of I_Arch conductance
			net[i].gChR2=gchr20+(gchr21-gchr20)*rnd();//rnd();
			

			
			//assign gnap,gleak and gcan values to rhythm generating population
			if (i<NPM && size>2){		
			double mu_nap=3.33;
			double s_nap=0.75;
			std::normal_distribution<double> dist_gnap_PM (mu_nap,s_nap);
			net[i].gNaP=dist_gnap_PM(generator);		
			double M_gl=exp((kbath-A_gl)/B_gl);
			double s_gl=sL_pct*M_gl;
			double mu_gl=M_gl+rho*(s_gl/s_nap)*(net[i].gNaP-mu_nap);
			double sigma_gl=sqrt((1-rho*rho)*s_gl*s_gl);
			
			std::normal_distribution<double> dist_gl_PM (mu_gl,sigma_gl);
			net[i].g_L=dist_gl_PM(generator);	
			net[i].gCAN=0.0;
			}
			
			//assign gnap,gleak and gcan values to pattern generating population
			if (i>=NPM && size>2){
			net[i].g_L=gl0+(gl1-gl0)*rnd();
			std::normal_distribution<double> dist_gnap_nPM (net[i].g_L*0.28,net[i].g_L*0.03);
			double mu_nap=1.5;
			double s_nap=0.25;
			std::normal_distribution<double> dist_gnap_PM (mu_nap,s_nap);
			net[i].gNaP=dist_gnap_PM(generator);
			
			double M_gl=exp((kbath-A_gl)/B_gl);
			double s_gl=0.025*M_gl;
			double mu_gl=M_gl+rho*(s_gl/s_nap)*(net[i].gNaP-mu_nap);
			double sigma_gl=sqrt((1-rho*rho)*s_gl*s_gl);
			std::normal_distribution<double> dist_gl_PM (mu_gl,sigma_gl);
			net[i].g_L=dist_gl_PM(generator);
			
			std::normal_distribution<double> dist_gcan(gcan0,gcan1);
			net[i].gCAN=dist_gcan(generator);
			}
			
			//2 NEURON NETWORK
			if (i==0 && size==2)
			{
			net[i].gCAN=0.0;
			net[i].gNaP=3.33;
			net[i].g_L==gl0+(gl1-gl0)*rnd();
			}
			if (i==1 && size==2)
			{
			net[i].gCAN=1.5;
			net[i].gNaP=1.5;
			net[i].g_L==gl0+(gl1-gl0)*rnd();
			}
			
		}
		
		w=new connection [size*size];
		sex=0;
		for(int i=0;i<size;i++) for(int j=0;j<size;j++)
		{
			/////////////////////
			/////////////////////
			//Neurons 0-24=inhib;25-99=PM; 100-end=NPM
			////////////////////
			////////////////////
			
			double conprob = prob;
			double tmpweight = wght;
			
			
			
			//nonPM to PM
			double Pnpm2pm = prob;//Prob of PM to nonPM
			
			if(size>2){
			//nonPM to nonPM
			if(i>=NPM && j>=NPM) {tmpweight = Wnpm2npm; conprob=Pnpm2npm;}
			//PM to nonPM
			if(i<NPM && j>=NPM) {tmpweight = Wpm2npm; conprob=Ppm2npm;}
			//nonPM to PM
			if(i>=NPM && j<NPM) {tmpweight = Wnpm2pm; conprob=Pnpm2pm; }
			//PM to PM
			if(i<NPM && j<NPM) {tmpweight = wght; conprob=prob;}
			}

			//std::default_random_engine generator;
			if(i==j || rnd()>conprob) continue;
			w[sex].source=i;
			w[sex].target=j;
			w[sex].weight=tmpweight*rnd();
			if(size==2){w[sex].weight=tmpweight;}
			sex++;
		}
	}

	~Population() { delete w; delete net; }

	int step(double dt,double drive,int* spk,double fcan,double fw, double fsynCa, double irr)
	{
		int sp=0;
		for(int i=0;i<size;i++) { spk[i]=net[i].spike(vth,dt,drive,fcan,fsynCa,irr); sp+=spk[i]; }
		for(int i=0;i<sex;i++){	 
					 if(spk[w[i].source]&& w[i].source<NPM) {net[w[i].target].gsyn+=w[i].weight*fw*musyndep*net[w[i].source].D_syn;}//if source is a PM neuron (includes opioid reduction of synapses)
					 else if(spk[w[i].source])              {net[w[i].target].gsyn+=w[i].weight*fw*net[w[i].source].D_syn;}//if source is a follower neuron
					 
					}
		//If synaptic depression is on
		if(D_on_off==1){
				for(int i=0;i<size;i++) if(spk[i]) net[i].D_syn += -f_D*net[i].D_syn; //if neuron i spikes then update D_syn
				}
		return sp;
	}
};

//initial conditons
Neuron::Neuron()
{
		v=-60+5*rnd(); m=.1; h=.1; mp=.1; hp=.4; n=.1; Cain=.000001; Cav=5*1e-5; mc=.1; hc=.1;
		gCa=.003;
		gNaP=5;
		gCAN=1;
		gChR2=0.4;
		gArch=0.4;
		E_leak=-68; g_L=2.5;
		hcan=1;
		OP1 =0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
		Carch = 0.9; Oarch = 0.0;
		m_ER=0.1; h_ER=0.99;
		Kout = 4;
		D_syn=D_0;
		caER=0.002; //catot=0.0005;
}

int Neuron::spike(double vth,double dt,double drive,double fcan, double fsynCa, double irr)
{
		double vpre=v;
		step(dt,drive,fcan,fsynCa, irr);
		return (vpre<vth && v>=vth);
}

void Neuron::init(double gca,double glk,double gnap,double gcan,double gchr2,double garch)
{
	gCa=gca;
	g_L=glk;
	gNaP=gnap;
	gCAN=gcan;
	gChR2=gchr2;
	gArch=garch;
	m=0.1*rnd();
	h=0.1*rnd();
	mp=0.1*rnd();
	hp=0.1*rnd();
	mc=0.1*rnd();
	hc=rnd();
	n=rnd();
	Cain=0.000001;
	caER=0.002; 
	Cav=5e-5*rnd();
	OP1=0.00; OP2=0.00; CL1=0.99; CL2=0.01; Pchr2=0.1;
	Carch = 0.9; Oarch = 0.0;
	
}

int flag_hcan=0;

void Neuron::step(double dt,double drive,double fcan,double fsynCa, double irr)
{
	double ECa=26.54*log(Caout/Cain)/2;
	EK=26.54*log(Kout/Kin);
	double pna = 1, pk=42;
	E_leak = -26.54*log((pna*Nain + pk*Kin)/(pna*Naout + pk*Kout));
	double INaF = gNaF*m*m*m*h;
	double INaP = gNaP*mp*hp*GNAPBLK;
	double ICa = gCa*mc*hc+gCaL;
	double IChR2 = gChR2*Gv(v)*(OP1+gma_chr2*OP2)*(v-Echr2);
	double IArch = gArch*Oarch*(v-Earch);
	double k0 = v+nAV;
	double k1 = nA*k0/(1-exp(-k0/nAk));
	double k2 = nB*exp(-(v+nBV)/nBk);
	double taun_inf = 1/(k1+k2);
	double n_inf = k1/(k1+k2);
	double IKdr = gKdr*n*n*n*n;
	double ICAN = gCAN*fcan/(1.+pow(K_CAN/Cain,nc))*hcan; // Tporikova	
	double Istim = 0;
	
	
	if(stimIDs==1){Istim =stim_W;}

	double Iinh=ginh*(v-Einh);
	if(ID<=NPM){Iinh=0;}
	
	//opioid DEP POTASSIUM CURRENT Imu
	double Imu=0;
	if(ID<NPM){Imu=gmu*(v-EK);}
	

	
	v+=(-INaF*(v-ENa)-INaP*(v-ENa)-IKdr*(v-EK)-ICAN*(v-ECAN)-ICa*(v-ECa)-IChR2-IArch-g_L*(v-E_leak)-(drive+gsyn)*(v-Esyn)-Iinh - Imu + IAPP + Istim)/c*dt;

	//INa
	m+=(x_inf(v,mV12,mk)-m)*(1.-exp(-dt/tau_inf(mtau,v,mtauV12,mtauk)));
	h+=(x_inf(v,hV12,hk)-h)*(1.-exp(-dt/tau_inf(htau,v,htauV12,htauk)));
	//INaP
	mp+=(x_inf(v,mpV12,mpk)-mp)*(1.-exp(-dt/tau_inf(mptau,v,mptauV12,mptauk)));
	hp+=(x_inf(v,hpV12,hpk)-hp)*(1.-exp(-dt/tau_inf(hptau,v,hptauV12,hptauk)));
	//Ica
	mc+=(x_inf(v,mcV12,mck)-mc)*(1.-exp(-dt/mctau));
	hc+=(x_inf(v,hcV12,hck)-hc)*(1.-exp(-dt/hctau));
	//IK
	n+=(n_inf-n)*(1.-exp(-dt/taun_inf));
	//IP3
	m_ER=(Cain*IP3)/((Cain+kaa)*(IP3+ki));
	h_ER=h_ER+(lp*(kd-(Cain+kd)*h_ER))*dt;
	J_ER_in =(gL_ER + gER*m_ER*m_ER*m_ER*h_ER*h_ER*h_ER)*(caER-Cain);
	J_ER_out = gSERCA*Cain*Cain/(kSERCA*kSERCA+Cain*Cain);
	
	if(ER_on_off==0){m_ER=0;J_ER_in=0;J_ER_out=0;}

	//Ca2+ dynamics
	Cain+=(fca*(J_ER_in-J_ER_out) - alphaCa*ICa*(v-ECa) - alphaCa*(gsyn)*fsynCa*(v-ECa) + 0*alphaCa*fsynCa*Istim +alphaCa*Iapp_ca + (Cain0-Cain)/tauCa)*dt;
	catot+=(-alphaCa*ICa*(v-ECa) - alphaCa*(gsyn)*fsynCa*(v-ECa) + 0*alphaCa*fsynCa*Istim  +alphaCa*Iapp_ca + (Cain0-Cain)/tauCa)*dt;
	caER=(catot-Cain)/sigma_ca;
	
	//Make sure that concentrations cant go below 0;
	if(Cain<0){Cain=0;}
	if(catot<0){catot=0;}
	
	Kout=kbath;
	gsyn*=exp(-dt/tausyn);
	D_syn+=((D_0-D_syn)/tau_D)*dt;//Exponential decay back to D_syn=1.
}

using namespace std;

ostream& operator <<(ostream& os,Neuron& N)
{
	return (os<<N.v<<'\t'<<N.m<<'\t'<<N.h<<'\t'<<N.mp<<'\t'<<N.hp<<'\t'<<N.n<<'\t'<<N.Cain);
}

istream& operator >>(istream& is,Neuron& N)
{
    return (is>>N.v>>N.m>>N.h>>N.mp>>N.hp>>N.n>>N.Cain);
}

ostream& operator <<(ostream& os,Population& p)
{
	for(int i=0;i<p.size;i++) os<<p.net[i]<<'\t';
	return os;
}

istream& operator >>(istream& is,Population& p)
{
	for(int i=0;i<p.size;i++) is>>p.net[i];
	return is;
}

int main(int argc,char** argv)
{
	double dt=.025,T=10000;
	double DT=20;
	int size=100;
	double gca[2]={0.0000065,0.0000065},gleak[2]={2.75,2.75},gnap[2]={3,5},gcan[2]={2,1},prob=.1,wght=.1;
	double dr[2]={0,1},fc[2]={1,1},fcw[2]={1,1}, fsca[2]={0.01,0.01}, tauDrug = -5000, pctblk = 0, irrpower[2]={0.0,0.0}, gchr2[2] = {0.0,0.0}, garch[2] = {0.0,0.0};
	double pulseON=0, pulsefreq=0.0, pulsedur=100,pfmax =.25, pfmin=.25, dl=0.1, numdl=0, stimnum= 100, seed=1, max_gnap_blk=0.5, APPramp[2] ={0.0,0.0};
	double stim_time=2000,strength = 25, Nstim=0,stimcount=0, stimseed =0, GINHIB=0,stim_delay=0, stim_advance=0;
	double caTon=250, caP=5000, caToff=caP-caTon, caAmp=0,CRC=0;
	int ini=0, block_type = -9, Nholostim = 0;
	char fname[256]="dat";
	int output=0;
	for(int i=1;i<argc;i++)
	{
		//Flags for changing model parameters
		if(strcmp(argv[i],"-DT")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-i")==0) ini=1;
		else if(strcmp(argv[i],"-hcan")==0) flag_hcan=1;
		else if(strcmp(argv[i],"-o")==0) output=1;
		else if(strcmp(argv[i],"-dt")==0) dt=atof(argv[++i]);
		else if(strcmp(argv[i],"-DT")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-s")==0) size=atoi(argv[++i]);
		else if(strcmp(argv[i],"-T")==0) T=atof(argv[++i])*1000;
		else if(strcmp(argv[i],"-d")==0) { dr[0]=atof(argv[++i]); dr[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-w")==0) wght=atof(argv[++i]);
		else if(strcmp(argv[i],"-tgs")==0) tausyn=atof(argv[++i]);
		else if(strcmp(argv[i],"-nap")==0) { gnap[0]=atof(argv[++i]); gnap[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-na")==0) gNaF=atof(argv[++i]);
		else if(strcmp(argv[i],"-kdr")==0) gKdr=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca")==0) { gca[0]=atof(argv[++i]); gca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-gleak")==0) { gleak[0]=atof(argv[++i]); gleak[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-cal")==0) gCaL=atof(argv[++i]);
		else if(strcmp(argv[i],"-tca")==0) tauCa=atof(argv[++i]);
		else if(strcmp(argv[i],"-kp")==0) Kpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-vp")==0) Vpump=atof(argv[++i]);
		else if(strcmp(argv[i],"-kbath")==0) kbath=atof(argv[++i]);
		else if(strcmp(argv[i],"-prob")==0) prob=atof(argv[++i]);
		else if(strcmp(argv[i],"-can")==0) { gcan[0]=atof(argv[++i]); gcan[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fc")==0) { fc[0]=atof(argv[++i]); fc[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fw")==0) { fcw[0]=atof(argv[++i]); fcw[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-fsca")==0) { fsca[0]=atof(argv[++i]); fsca[1]=atof(argv[++i]); }
		else if(strcmp(argv[i],"-ca0")==0) Cain0=atof(argv[++i]);
		else if(strcmp(argv[i],"-exp")==0) block_type=atof(argv[++i]); //1=can 2=TRPC3 
		else if(strcmp(argv[i],"-tauDrg")==0) tauDrug=atof(argv[++i]);//tau drug
		else if(strcmp(argv[i],"-sd")==0) seed=atof(argv[++i]);
		else if(strcmp(argv[i],"-scalegNaK")==0) {gNaF=gNaF*atof(argv[++i]); gKdr=gKdr*atof(argv[++i]);}
		else if(strcmp(argv[i],"-na_h12")==0) htauV12=atof(argv[++i]);
		else if(strcmp(argv[i],"-mtau")==0) mtau=atof(argv[++i]);
		else if(strcmp(argv[i],"-ktau")==0) nA=atof(argv[++i]);
		else if(strcmp(argv[i],"-Cm")==0) c=atof(argv[++i]);
		else if(strcmp(argv[i],"-K_v12")==0) {nAV=atof(argv[++i]); nBV=atof(argv[++i]);}
		else if(strcmp(argv[i],"-hp12")==0) hpV12=atof(argv[++i]);
		else if(strcmp(argv[i],"-iapp")==0) {APPramp[0]=atof(argv[++i]); APPramp[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-hptauk")==0) hptauk=atof(argv[++i]);
		else if(strcmp(argv[i],"-ttxblk")==0) GNAPBLK=atof(argv[++i]);
		else if(strcmp(argv[i],"-dsyn")==0) D_on_off=atof(argv[++i]);
		else if(strcmp(argv[i],"-ERoff")==0) ER_on_off=atof(argv[++i]);
		else if(strcmp(argv[i],"-CRC")==0) CRC=atof(argv[++i]);
		else if(strcmp(argv[i],"-spotstimIDs")==0) {cells[0]=atof(argv[++i]); cells[1]=atof(argv[++i]);}
		else if(strcmp(argv[i],"-spotstimW")==0) strength=atof(argv[++i]);
		else if(strcmp(argv[i],"-spotstimT")==0) stim_time=atof(argv[++i]);
		else if(strcmp(argv[i],"-NUMspotstim")==0) Nstim=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqONOFF")==0) square_ON_OFF=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqtONtOFF")==0) {caTon=atof(argv[++i]); caToff=atof(argv[++i]);}
		else if(strcmp(argv[i],"-ca_sqAMP")==0) caAmp=atof(argv[++i]);
		else if(strcmp(argv[i],"-ca_sqP")==0) caP=atof(argv[++i]);
		else if(strcmp(argv[i],"-gmu")==0) gmu=atof(argv[++i]);
		else if(strcmp(argv[i],"-Agl")==0) A_gl=atof(argv[++i]);
		else if(strcmp(argv[i],"-Bgl")==0) B_gl=atof(argv[++i]);
		else if(strcmp(argv[i],"-w12")==0) Wpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-p22")==0) Pnpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-w22")==0) Wnpm2npm=atof(argv[++i]);
		else if(strcmp(argv[i],"-w21")==0) Wnpm2pm=atof(argv[++i]);
		else if(strcmp(argv[i],"-syn_step")==0) f_D=atof(argv[++i]);
		else if(strcmp(argv[i],"-tau_syn")==0) tau_D=atof(argv[++i]);
		else if(strcmp(argv[i],"-GSERCA")==0) gSERCA=atof(argv[++i]);
		else if(strcmp(argv[i],"-GER")==0) gER=atof(argv[++i]);
		else if(strcmp(argv[i],"-gldist")==0) sL_pct=atof(argv[++i]);
		else if(strcmp(argv[i],"-GINHIB")==0) GINHIB=atof(argv[++i]);
		else if(strcmp(argv[i],"-kSERCA")==0) kSERCA=atof(argv[++i]);
		else if(strcmp(argv[i],"-gL_ER")==0) gL_ER=atof(argv[++i]);
		else if(strcmp(argv[i],"-Nholo")==0) Nholostim=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimseed")==0) stimseed=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimdelay")==0) stim_delay=atof(argv[++i]);
		else if(strcmp(argv[i],"-stimadvance")==0) stim_advance=atof(argv[++i]);
		else if(strcmp(argv[i],"-musyndep")==0) musyndep=atof(argv[++i]);//#should be between 1 and 0.5
		else strcpy(fname,argv[i]);
	}






	if(seed !=0) srand (seed);
	int freq=int(DT/dt);


	pulseON = (1000/pulsefreq);
	Population pop(size,gca[0],gca[1],gleak[0],gleak[1],gnap[0],gnap[1],gcan[0],gcan[1],gchr2[0],gchr2[1],garch[0],garch[1],prob,wght,stimnum);
    	if(ini) { ifstream is("ini"); is>>pop; }
	ofstream out(fname);
	
	//////////////
	//////////////
	int sp=0;	
	int h[size];
	double last_sp_t[size];
	for(int i=0;i<size;i++) last_sp_t[i]=0;
	double irr=0;
        double PM_sp=0,nPM_sp=0;
        double t_last_stim=-5000;
        double t_Ca_last=0;
        double t_ca_on=0;
        double t_after_B;//time after last burst
        //////////////
	//////////////
        
	caToff=caP-caTon;
	//set seed for stimulation
	if(stimseed !=0) srand (stimseed);
	int count=0;
	while(count<Nholostim){
				pop.net[(rand() % size)].stimIDs=1;
				count=0;
				//count the number of stimulated neurons
				for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){count++;}}
				}
	cerr<<Nholostim<<'\t'<<stimseed<<'\t'<<endl;
	for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){cerr<<i<<'\t'<<pop.net[i].stimIDs<<'\t'<<endl;}}
	//set seed back to what it was
	if(seed !=0) srand (seed);

	//write gLeak,gNaP,CellID,kbath values to .ca file
	for(int i=0;i<size;i++)
	{
	cerr<<pop.net[i].g_L<<'\t'<<pop.net[i].gNaP<<'\t'<<i<<'\t'<<kbath<<'\t'<<endl;
	}


	///////////////////////////////////////////////////////
	//Main Simulation Loop
	//////////////////////////////////////////////////////
	for(double t=0;t<=T;t+=dt)
	{
		//Parameters that can vary during simulations
		double drive=dr[0]+t/T*(dr[1]-dr[0]);//synaptic drive
		double fcan= fc[0]+t/T*(fc[1]-fc[0]);//Can conducatance
		double fw=fcw[1];//-t/T*(fcw[1]-fcw[0]);//Synaptic weights
		double fsynCa=fsca[0]+t/T*(fsca[1]-fsca[0]);//Psynca
		IAPP=APPramp[0]+(t)/(T)*(APPramp[1]-APPramp[0]);//Ramp in applied current
		
		//Ca square wave 
		if(square_ON_OFF==1)
			{	
			Iapp_ca=0;
			if(t_Ca_last>caToff){Iapp_ca=caAmp; t_ca_on+=dt;}
			if(t_ca_on>caTon){t_Ca_last=0; t_ca_on=0; Iapp_ca=0;}
			t_Ca_last+=dt;
			}
		
		//spot stimulation
		if(t_last_stim>stim_time && t_after_B>=stim_delay && t_after_B<(stim_delay+dt) && stimcount<Nstim){
							     stim_ON=1;stim_W=strength;stimcount+=1;t_last_stim=0;
							     stim_delay = stim_delay+stim_advance;
//							     for(int i=0;i<(cells[1]-cells[0]);i++) cerr<<(t/1000)<<'\t'<<cells[0]+i<<'\t'<<-1<<'\t'<<endl;
							     for(int i=0;i<size;i++){if(pop.net[i].stimIDs==1){cerr<<(t/1000)<<'\t'<<i<<'\t'<<-1<<'\t'<<endl;}}
							    }					    
		t_after_B+=dt;
		stim_W*=exp(-dt/100);
		t_last_stim+=dt;
		//END spot stimulation


		//Array for cells that spike during timestep
		int spk[size];
	

		sp+=pop.step(dt,drive,spk,fcan,fw,fsynCa,irr);

		for(int i=0;i<size;i++)
		{
			h[i]+=spk[i];
			if(h[i]>=1)
			{
				//cout<<(t/1000)<<'\t'<<i<<'\t'<<endl;//uncomment this line to write spike times to .hst file for plotting raster plot
				
			}

		}
		

		
		for(int i=0;i<size;i++) h[i]=0;
		


		if(int(t/dt)%freq==0) 
		{	

			if(size>2){out<<(t/1000)<<'\t'<<sp/(.001*DT*size)<<'\t'<<endl;}
			if(size==2){out<<(t/1000)<<'\t'<<pop.net[0].v<<'\t'<<pop.net[1].v<<'\t'<<IAPP<<'\t'<<endl;}
			if(size==1){out<<(t/1000)<<'\t'<<pop.net[0].v<<'\t'<<Iapp_ca<<'\t'<<pop.net[0].Cain<<'\t'<<pop.net[0].caER<<'\t'<<kbath<<'\t'<<endl;}
			sp=0;
	
		}

	}
   	ofstream os("ini"); os<<pop;
	return 0;
}

