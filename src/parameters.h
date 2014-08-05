/************************
gpwexact5.F Parameters:

Requested at run-time:
sigqf             Standard deviation of fixed charge (cm^-2) (normalized by q)
wtunnel           Width of buried trap distribution (cm)
outfile           Name of file to which results are written

Hard coded:
T                 Temperature (¡C)
nd                Doping (cm^-3) (assumes n-type)
tox               Oxide thickness (cm)
qfbar             Fixed charge (C/cm^2)

Simulation range:
nfreq             Number of frequency data points
ftart             Lowest frequency (hard coded w/ 10 samples/decade log spaced)

nufs              Number of surface potential data points
etmin             Minimum trap level (as measured below Ec)
etmax             Maximum trap level (as measured below Ec)

Interface trap model:
dmin              Value of dit in the middle of the bandgap (cm^-2 eV^-1)
dc0               Extrapolated value of Dit at midgap, upper half of bandgap
betac             Exponential slope factor for Dit, upper half of bandgap
dv0               Extrapolated value of Dit at midgap, lower half of bandgap
betav             Exponential slope factor for Dit, lower half of bandgap
da                Donor-like flag: if 1, upper half donor, lower half acceptor
                                   if 0, lower half donor, lower half acceptor

Capture cross section model:
signmax           Value of capture cross section at midgap
signcb0           Extrapolated value of sigman at midgap, upper half of bandgap
gammacb           Exponential slope factor for sigman, upper half of bandgap                                   
signvb0           Extrapolated value of sigman at midgap, lower half of bandgap
gammavb           Exponential slope factor for sigman, upper half of bandgap                                   
************************/


#define GPW_FILENAME    "data/gpw.csv"
#define DATA_FILENAME   "data/data.csv"
#define PARAM_FILENAME  "data/param.csv"        
#define USVG_FILENAME   "data/usvg.csv"

#define TEMPERATURE     21.0        		/* T                    */

#define DOPING          1e16         		/* nd                   */
#define T_OXIDE         20e-7       		/* tox                  */
#define FIXED_CHARGE    1e12        		/* qfbar                */
#define SIGMA_QF        (FIXED_CHARGE*0.2)	/* sigqf                */

#define TUNNEL_WIDTH    0             		/* wtunnel              */

#define N_FREQ          151         		/* nfreq                */
#define START_FREQ      1e-3        		/* fstart               */

#define N_ET            11          		/* nufs                 */
#define ET_MIN          0.2         		/* etmin                */
#define ET_MAX          0.8         		/* etmax                */

#define DIT_MIDGAP      0.0           		/* dmin                 */
#define DIT_C0          8.48e6        		/* dc0                  */
#define DIT_BETA_C      0.207       		/* betac                */
#define DIT_V0          8.48e6      		/* dv0                  */
#define DIT_BETA_V      0.207       		/* betav                */
#define DIT_DA          1           		/* da                   */

#define SIGN_MAX        1e12             		/* signmax              */
#define SIGN_CB0        1e-16        		/* signcb0              */
#define SIGN_GAMMA_CB   0           		/* gammacb              */
#define SIGN_VB0        1e-16        		/* signcb0              */
#define SIGN_GAMMA_VB   0           		/* gammacb              */