import java.io.*;

public class smolenltp 
/* To compile and run this file, save as smolenltp.java. Then with java compiler installed (e.g. JDK from Oracle), type “javac smolenltp.java” to compile and “java smolenltp” to run. Output will be ASCII 2-column text files listed below. */

 {  public static void main (String args[]) throws IOException {

// THESE FILES OUTPUT MODEL VARIABLES OF INTEREST

    PrintWriter out1 = new PrintWriter(new FileWriter("CACYT.txt"));
    PrintWriter out2 = new PrintWriter(new FileWriter("RAF.txt"));
    PrintWriter out3 = new PrintWriter(new FileWriter("MKKPP.txt"));
    PrintWriter out4 = new PrintWriter(new FileWriter("ERKACT.txt"));
    PrintWriter out5 = new PrintWriter(new FileWriter("PCREB.txt"));
    PrintWriter out6 = new PrintWriter(new FileWriter("GPROD.txt"));
    PrintWriter out7 = new PrintWriter(new FileWriter("PKAACT.txt"));
    PrintWriter out8 = new PrintWriter(new FileWriter("TAGSYN.txt"));
    PrintWriter out9 = new PrintWriter(new FileWriter("WSYN.txt"));
    PrintWriter out10 = new PrintWriter(new FileWriter("PERK.txt"));
    PrintWriter out11 = new PrintWriter(new FileWriter("CAMKIV.txt"));
    PrintWriter out12 = new PrintWriter(new FileWriter("CAMKII.txt"));
    PrintWriter out13 = new PrintWriter(new FileWriter("TAGP1.txt"));
    PrintWriter out14 = new PrintWriter(new FileWriter("ERKNUC.txt"));
    PrintWriter out15 = new PrintWriter(new FileWriter("CAMKK.txt"));
    PrintWriter out16 = new PrintWriter(new FileWriter("CAMP.txt"));
    PrintWriter out17 = new PrintWriter(new FileWriter("TAGP2.txt"));
    PrintWriter out18 = new PrintWriter(new FileWriter("RAFP.txt"));
    PrintWriter out19 = new PrintWriter(new FileWriter("CANUC.txt"));
    PrintWriter out20 = new PrintWriter(new FileWriter("TAGP3.txt"));
    PrintWriter out21 = new PrintWriter(new FileWriter("PREW.txt"));
    PrintWriter out22 = new PrintWriter(new FileWriter("ERKSOMACT.txt"));

    double dt=0.0002; // high-res timestep (min)
    double delta=0.01; // low-res timestep (min)

    double recstart=3950.0; // Time to start writing data
    double stime=4000.0; // Time of stimulus onset

    double isi=5.0;      // Spacing between tetanic stimuli

//Time to end simulation and stop writing data
    double recend=4450.0; 

//    double recend=4000.0+24.0*60.0+2.0*isi; 

    int recintvl=10; // (delta) Intervals to record

    double tref; // Time since stimulus onset




//PROGRAM VARIABLES – nomenclature as in Smolen et al 2006.

// Activities of kinases other than MAPK cascades

//CAMKIV
	double ck4act;
//CAMKK
	double ckkact;
	double ckkbas;
// CAMKII
	double ck2act;
// camp and PKA
	double camp;
	double pkaact;

// MAPK CASCADE
	double raf;
	double rafp;
	double mkk;
	double mkkp;
	double mkkpp;
	double erk;
	double erkp;
	double erkpp;
	double erkact;
	double erks;
	double erkps;
	double erkpps;
	double erkacts;
	double erknuc;

/* Sites of translation factor phosphorylation, and translated protein GPROD. */
	double pcreb;
	double perk;
	double gprod;

// SYNAPTIC VARIABLES
	double tagp1;
	double tagp2;
	double tagp3;
	double tltp;
	double tagsyn;

/* “prew” is variable needed to keep W from going too high with multiple stimuli. In Smolen et al 2006 this variable is just denoted “P”. */
	double prew;

/* Synaptic weight variable wsyn. In Smolen et al. 2006 this variable is just denoted “W”. */
	double wsyn;

	
// TIME DERIVATIVES

    	double dv1dt;
    	double dv2dt;
    	double dv3dt;
    	double dv4dt;
    	double dv5dt;
    	double dv6dt;
    	double dv7dt;
    	double dv8dt;
    	double dv9dt;
    	double dv10dt;
    	double dv11dt;
    	double dv12dt;
    	double dv13dt;
    	double dv14dt;
	double dv15dt;
	double dv16dt;
	double dv17dt;
	double dv18dt;
	double dv19dt;
	double dv20dt;
	double dv21dt;

/* Cytoplasmic and nuclear Ca concentrations, followed by power functions used to describe activation of CaMKII, CaMKK, and CAMKIV by Ca/CaM, followed by camp Hill function for PKA activation. Explicit formulas for functions given below. */
	double cacyt;
	double canuc;
	double powca;
	double powkc;
	double powcan;
	double powkcn;
	double termcamp;

// WORK IN CONCENTRATION UNITS OF uM, TIME UNITS OF MIN

// MODEL PARAMETERS FOLLOW:

// BASAL PKA AND CA, STIMULUS AMPLITUDES. In the following, CA2 refers to nuclear Ca.
// Stim amps also drive PKA elevations, by driving step CAMP increases
// THE AMPTETSTIM, AMPTHETASTIM, AMPKANDSTIM refer to brief step activations of Raf,following stimuli.

	double Cabas=0.04;
	double Cabasnuc=Cabas;
	double campbas=0.05;
	double AMPTETCA=1.0;
	double AMPTHETACA=1.0;
	double AMPKANDCA=0.4;
	double AMPCHEMCA=0.1;

	double AMPTETCA2=0.5;
	double AMPTHETACA2=0.5;
	double AMPKANDCA2=0.18;
	double AMPCHEMCA2=0.1;

	double AMPTETCMP=0.15;
	double AMPTHETACMP=0.35;
	double AMPKANDCMP=0.15;
	double AMPCHEMCMP=0.4;

	double AMPTETSTIM=0.15;
	double AMPTHETASTIM=0.4; // 3 TIMES TET TO COMP FOR 3 TETS but only one THETA
	double AMPKANDSTIM=0.15;  // SHOULD BE NO MORE THAN TET

	double AMPCHEMSTIM=0.3;

	double ampstim=0.0;

// RATE CONSTANTS AND MICHAELIS CONSTANTS FOR CAMKIV AND CAMKII 
	double kfck4=10.0;
	double tauck4=20.0;
	double kfckk=2.5;
	double tauckk=0.2;
	double kfck2=200.0;
	double tauck2=1.0;
	double Kck2=0.7;
	double Kck2n=0.3;

//CONSERVED ENZYME AMOUNTS
	double raftot=0.25;
	double mkktot=0.25;
	double erktot=0.25;

// MAPK rate and Michaelis constants
	double kbraf=0.12;
	double kfmkk=0.6;
	double kbmkk=0.025;
	double Kmkk=0.25;
	double kferk=0.52;
	double kberk=0.025;
	double Kmk=0.25;
	double kfbasraf=0.0075;
	double knucin=100.0;
	double knucout=2.5;

// Maximal and basal transcription rates for GENEP1, decay of GENEP1
	double ktrans1=1.0;
	double vtrbas1=0.0004;
	double taug1=100.0;

// "Binding-Dissoc" constants describing interactions of
// TFs with promoters. Separate constants for CREB, CREB2 
// (i.e. Kerk1), CBP (i.e. KCBP).
	double Kcreb=1.0;
	double Kerk1=1.0;

// Rate constants and other parameters
// describing PKA activation and synaptic tagging

	double Kcamp=0.5;

	double kphos1=0.05;
	double kdeph1=0.02;
	double kphos2=2.0;
	double kdeph2=0.2;
	double kphos3=0.06;
	double kdeph3=0.017;
	double taupka=15.0;

// TF PHOS AND DEPHOS RATE PARAMS
	double kphos4=0.12;
	double kdeph4=0.03;
	double kphos5=4.0;
	double kdeph5=0.1;

// Inhibition factors for enzymes and GeneProd1. Usually set to 1 below, values lower than 1 mean application of inhibitors

	double inherk;
	double inhmkk;
	double inhck4;
	double inhck2;
	double inhgp;
	double inhpka;

// Rate constants and other parameters for synaptic 
// weight changes, that is, changes in variable Wsyn

	double kltp=2.0;
	double tausyn=140000.0;
	double Kprew=0.03;
	double vprew=0.0003;
	double tauprew=1000.0;

// TIME STUFF, COUNTERS, VALUE ARRAY

    double time=0.0; // time (minutes)
    double timewrite;

    int i, k, j, l; // counters

/* values array holds dynamic values of model variables with differential equations */
    double[] values=new double[23];

// VARIABLE INITIALIZATION. I like to use small but nonzero initial values to avoid extremely small numbers in internal computations.

	values[1]=0.001;
	values[2]=0.001;
	values[3]=0.001;
	values[4]=0.001;

/* For MAPK cascade variables initial values have to be set to be fractions of the total amounts of the respective kinases. Conservation condition equations below will ensure the remaining portions of total kinases are allocated to the variables that do not have their own differential equations */
	values[5]=0.5*raftot;
	values[6]=0.3*mkktot;
	values[7]=0.4*mkktot;
	values[8]=0.3*erktot;
	values[9]=0.4*erktot;
	values[10]=0.3*erktot;
	values[11]=0.4*erktot;
	values[12]=0.5*values[11];

	values[13]=0.001;
	values[14]=0.001;
	values[15]=0.001;
	values[16]=0.001;
	values[17]=0.001;
	values[18]=0.01;

/* Variables 19 and 20 (synaptic weight and limiting protein P) are slowest in the model. So use initial values that are relatively close to their equilibrium basal values */
	values[19]=0.1;
	values[20]=vprew*tauprew;

// MAIN LOOP (LARGER TIMESTEP)

        k=1;
        do {  

// INNER SIMULATION LOOP

            j=1;
            do {

		tref=time-stime;
		cacyt=Cabas;
		canuc=Cabasnuc;
		ampstim=0.0;
		camp=campbas;

		ckkact=values[1];
		ck4act=values[2];
		ck2act=values[3];
		pkaact=values[4];

		raf=values[5];
		mkk=values[6];
		mkkpp=values[7];
		erk=values[8];
		erkpp=values[9];

		erks=values[10];
		erkpps=values[11];
		erknuc=values[12];

		tagp1=values[13];
		tagp2=values[14];
		tagp3=values[15];

		pcreb=values[16];
		perk=values[17];
		gprod=values[18];

		wsyn=values[19];
		prew=values[20];

// SIMPLE SQUARE WAVES FOR CYTOPLASMIC CA
// DURING THE TETANIC, THETA BURST, AND KANDEL PROTOCOLS. 

// Tetanic protocol

		if (tref > 0.0 && tref < 0.05)
		  {
		  cacyt=AMPTETCA;
		  canuc=AMPTETCA2;
		  }
		if (tref > 0.0 && tref < 1.0)
		  {
		  camp=AMPTETCMP;
		  }
		if (tref > 0.0 && tref < 1.0)
		  {
		  ampstim=AMPTETSTIM;
		  }
		if (tref > (0.0+isi) && tref < (0.0+isi+0.05))
		  {
		  cacyt=AMPTETCA;
		  canuc=AMPTETCA2;
		  }
		if (tref > (0.0+isi) && tref < (0.0+isi+1.0))
		  {
		  camp=AMPTETCMP;
		  }
		if (tref > (0.0+isi) && tref < (0.0+isi+1.0))
		  {
		  ampstim=AMPTETSTIM;
		  }
		if (tref > (0.0+2.0*isi) && tref < (0.0+2.0*isi+0.05))
		  {
		  cacyt=AMPTETCA;
		  canuc=AMPTETCA2;
		  }
		if (tref > (0.0+2.0*isi) && tref < (0.0+2.0*isi+1.0))
		  {
		  camp=AMPTETCMP;
		  }
		if (tref > (0.0+2.0*isi) && tref < (0.0+2.0*isi+1.0))
		  {
		  ampstim=AMPTETSTIM;
		  }

// Theta burst protocol. To simulate theta burst LTP instead of tetanic LTP, uncomment this protocol and comment out the tetanic one.
/*
		if (tref > 0.0 && tref < 0.084)
		  {
		  cacyt=AMPTHETACA;
		  canuc=AMPTHETACA2;
		  }
		if (tref > 0.0 && tref < 1.0)
		  {
		  camp=AMPTHETACMP;
		  }
		if (tref > 0.0 && tref < 1.0)
		  {
		  ampstim=AMPTHETASTIM;
		  }
*/

// Chem LTP
/*

// Activates ERK but not PI3K (maybe lack of PI3K activation explains
// why no ELTP results from Chem treatment).
		if (tref > 0.0 && tref < 30.0)
		  {
		  ampstim=AMPCHEMSTIM;
		  camp=AMPCHEMCMP;
		  }

// Have to assume Chem treatment increases neuron excitability some and
// lets in some extra Ca, otherwise I dont get decent LLTP.
		if (tref > 0.0 && tref < 30.0)
		  {
		  cacyt=AMPCHEMCA;
		  canuc=AMPCHEMCA2;
		  }
*/

// INHIBITOR APPLICATIONS. FIRST SET THE BASELINE NONINHIBITED FACTOR VALUES
// TO 1.

		inherk=1.0;
		inhck4=1.0;
		inhck2=1.0;
		inhgp=1.0;
		inhpka=1.0;
		inhmkk=1.0;

/* Following group of statements allows for inhibitor to be applied for a specific interval (adjust numbers in parentheses to adjust interval). Note, the inhibitor “inhgp” inhibits protein synthesis, others inhibit enzymes. */
		if (tref > -0.001 && tref < 60.0)
		  {
//			inherk=0.1;
//			inhck4=0.1;
//			inhck2=0.1;
//			inhgp=0.4;
//			inhmkk=0.1;
//			inhpka=0.1;
		  }

/* For inhibitors present “throughout experiment” use these statements instead. Here inhibitor application starts at time of stimulus. Model cannot properly simulate inhibitor present much earlier than stimulus, because the model does not represent dynamics of basal synaptic weight maintenance. */
		if (tref > -0.001)
		  {
//			inherk=0.1;
//			inhck4=0.1;
//			inhck2=0.1;
//			inhgp=0.4;
//			inhpka=0.1;
		  }

// BEGINNING OF ACTUAL CODE FOR UPDATING DYNAMIC VARIABLES. 

/* Explicit formulas for power functions used in CaM kinase activation and for the Hill function describing PKA activation by Camp.
Powers of Ca concentrations are used in Hill functions in the appropriate differential equations */
		powca=cacyt*cacyt*cacyt*cacyt;
		powkc=Kck2*Kck2*Kck2*Kck2;
		powcan=canuc*canuc*canuc*canuc;
		powkcn=Kck2n*Kck2n*Kck2n*Kck2n;
		termcamp=camp*camp/(camp*camp+Kcamp*Kcamp);

		rafp=raftot-raf;
		mkkp=mkktot-mkk-mkkpp;
		erkp=erktot-erk-erkpp;
		erkps=erktot-erks-erkpps-erknuc;

/* Activation of CAM kinase kinase by nuclear Ca (CaM dynamics not modeled explicitly, only power of Ca concentration is used */

		dv1dt = kfckk*(powcan/(powcan+powkcn)) 
                        - ckkact/tauckk;

// Activation of CAMKIV by nuclear calcium and by CAMKK (whose difeq is below). 
// Nonlinear (power hill function) due to Ca interactions with calmodulin.

		dv2dt = kfck4*(powcan/(powcan+powkcn))*ckkact 
                       - ck4act/tauck4;

// Activation of CAMKII by a power of cytoplasmic calcium, hill function,
// uses fourth powers.

		dv3dt = kfck2*(powca/(powca+powkc)) - ck2act/tauck2;

		dv4dt = (termcamp - pkaact)/taupka;

// Synaptic and somatic MAPK cascades 
		dv5dt = -(ampstim+kfbasraf)*raf+kbraf*rafp;

		dv6dt = -inhmkk*kfmkk*rafp*mkk/(mkk+Kmkk)+kbmkk*mkkp/(mkkp+Kmkk);

		dv7dt = inhmkk*kfmkk*rafp*mkkp/(mkkp+Kmkk)-kbmkk*mkkpp/(mkkpp+Kmkk);

		dv8dt = -kferk*mkkpp*erk/(erk+Kmk)
			+kberk*erkp/(erkp+Kmk);

		dv9dt = kferk*mkkpp*erkp/(erkp+Kmk)
			 -kberk*erkpp/(erkpp+Kmk);

		erkact = erkpp;

		dv10dt = -kferk*mkkpp*erks/(erks+Kmk)
			+kberk*erkps/(erkps+Kmk);

		dv11dt = kferk*mkkpp*erkps/(erkps+Kmk)
			 -kberk*erkpps/(erkpps+Kmk)
			 -inhpka*pkaact*knucin*erkpps+knucout*erknuc;

		erkacts = erkpps;

// Reversible import of somatic active ERK into nucleus

		dv12dt = inhpka*pkaact*knucin*erkpps-knucout*erknuc;

// Tag phosphorylation and setting

		dv13dt = kphos1*inhck2*ck2act*(1.0-tagp1) - kdeph1*tagp1;

		dv14dt = inhpka*kphos2*pkaact*(1.0-tagp2) - kdeph2*tagp2;

		dv15dt = inherk*kphos3*erkact*(1.0-tagp3) - kdeph3*tagp3;

		tagsyn = tagp1*tagp2*tagp3;

// Phosphorylation of translation factor sites and synthesis of protein denoted gprod.
		dv16dt = kphos4*inhck4*ck4act*(1.0-pcreb)-kdeph4*pcreb;

		dv17dt = kphos5*inherk*erknuc*(1.0-perk)-kdeph5*perk;
 
		dv18dt = inhgp*ktrans1*(pcreb/(pcreb+Kcreb))
*(perk/(perk+Kerk1))
			 - gprod/taug1 + inhgp*vtrbas1;

/* Synaptic weight dynamics, increase when tagsyn and gprod are up and when limiting factor prew is present, and first order decay */

		dv19dt = kltp*tagsyn*gprod*(prew/(prew+Kprew)) - wsyn/tausyn;

/* Dynamics of limiting protein prew. Is used up during synaptic weight increase, so equation is similar to that for synaptic weight */

		dv20dt = -kltp*tagsyn*gprod*(prew/(prew+Kprew)) + vprew - prew/tauprew;

		dv21dt = 0.0;

		values[1]+=dt*dv1dt;
		values[2]+=dt*dv2dt;
		values[3]+=dt*dv3dt;
		values[4]+=dt*dv4dt;
		values[5]+=dt*dv5dt;
		values[6]+=dt*dv6dt;
		values[7]+=dt*dv7dt;
		values[8]+=dt*dv8dt;
		values[9]+=dt*dv9dt;
		values[10]+=dt*dv10dt;
		values[11]+=dt*dv11dt;
		values[12]+=dt*dv12dt;
		values[13]+=dt*dv13dt;
		values[14]+=dt*dv14dt;
		values[15]+=dt*dv15dt;
		values[16]+=dt*dv16dt;
		values[17]+=dt*dv17dt;
		values[18]+=dt*dv18dt;
		values[19]+=dt*dv19dt;
		values[20]+=dt*dv20dt;
		values[21]+=dt*dv21dt;

// Equilibrate basal synaptic weight, prior to stimulation,
// according to basal kinase activities

		if (tref < -60.0)
		  {
		  wsyn = kltp*tausyn*tagsyn*gprod;
		  values[19] = kltp*tausyn*tagsyn*gprod;
		  }

// Increment time
                time=time+dt;

// END INNER LOOP

                j++;
               } while (j <= delta/dt);

// COMPUTE AND PRINT OUTPUT VARIABLES

            if ((time > recstart) && (time < recend) && (k % recintvl == 0))
              {
		timewrite=tref;

//	OUTPUT CONCENTRATION UNITS ARE uM. SCALING FACTORS ARE FOR
//	EASE OF CONCURRENT VISUALIZATION

               out1.println(timewrite + "\t" + 1.0*cacyt);
               out2.println(timewrite + "\t" + 1.0*raf);
               out3.println(timewrite + "\t" + 5.0*mkkpp);
               out4.println(timewrite + "\t" + 5.0*erkact);
               out5.println(timewrite + "\t" + 1.0*pcreb);
               out6.println(timewrite + "\t" + 0.4*gprod); 
               out7.println(timewrite + "\t" + 1.0*pkaact);
               out8.println(timewrite + "\t" + 110.0*tagsyn);
               out9.println(timewrite + "\t" + 1.0*wsyn);
               out10.println(timewrite + "\t" + 1.0*perk);
               out11.println(timewrite + "\t" + 10.0*ck4act);
               out12.println(timewrite + "\t" + 0.1*ck2act);
               out13.println(timewrite + "\t" + 1.0*tagp1);
               out14.println(timewrite + "\t" + 5.0*erknuc); 
               out15.println(timewrite + "\t" + 1.0*ckkact);
               out16.println(timewrite + "\t" + 1.0*camp);
               out17.println(timewrite + "\t" + 1.0*tagp2);
               out18.println(timewrite + "\t" + 10.0*rafp);
               out19.println(timewrite + "\t" + 1.0*canuc);
               out20.println(timewrite + "\t" + 1.0*tagp3);
               out21.println(timewrite + "\t" + 1.0*prew);
               out22.println(timewrite + "\t" + 5.0*erkacts);

              }

/* Time intervals, with respect to stimulus, used to plot figures in Smolen et al. 2006 are as follows. To set these intervals, go to beginning of file and adjust parameters recstart and recend. For example, recstart set to 20 min before stime (stimulus time) and recend set to 60 min after stime is used to plot Rafp Camkii and Camkiv time courses in response to 3-tetanus protocol */

// FIG 3A, RAFP CAMKII CAMKIV - FROM -20 TO 60
// FIG 3B, ERK ERKNUC TAGSYN - FROM -30 TO 120
// FIG 3C, GPROD WS WS2 TAGSYN - PANEL FROM -50 TO 250

// FIG 4, ERKACT, WSYN2, INHIB WS2 FOR ERK, CK2, CK4, -50 TO 250

// FIG 5A, UNSCALED WS2 FOR ALL FOUR PROTOCOLS, -50 TO 250
// FIG 5B, ERK TAGSYN GENEP WS2 PKAINHIBWS2 FOR HUANG, -50 TO 200

// FIG 6A, ERK CK2 TAGSYN FOR CHEM, -30 TO 120
// FIG 6B, GPROD WS2 WS2 WITH MEK INHIB FOR CHEM, -50 TO 250

// FIG 7, SYN TAGGING, TAG STIM COMES 60 MIN AFTER NUC STIM, 
// -50 to 250, WS2 FOR BOTH, GPROD, CAMKIV, TAGSYN

            k++;
           } while (k <= recend/delta);

// END OF LOOPS AND OF SIMULATION. CLOSE OUTPUT FILES.

        out1.close();
        out2.close();
	out3.close();
	out4.close();
        out5.close();
	out6.close();
	out7.close();
        out8.close();
	out9.close();
	out10.close();
	out11.close();
	out12.close();
	out13.close();
	out14.close();
	out15.close();
	out16.close();
	out17.close();
	out18.close();
	out19.close();
	out20.close();
	out21.close();
	out22.close();

       }
}
