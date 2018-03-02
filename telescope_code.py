# telescopecode
import numpy as np
import math
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import sys
from matplotlib.offsetbox import AnchoredText
import sys
import matplotlib.animation as animation
from scipy import signal


'''
DEFINITIONS OF THE FUNCTIONS USED
'''

def writetxt (list, name):
    
    #open the file to write
    datafile = open(name, 'w+')
    
    #stack lists vertically
    datalist = np.column_stack(list)
    
    #save data in specified file, close it
    np.savetxt(datafile, datalist, fmt = '%s')
    datafile.close()

#function to calculate chi squared - will be used for chi sq vs period plot
def chisq (m_obs, m_mod, err):
    chi = 0
    for i in range (len(m_obs)):
        chi+=((m_obs[i]-m_mod[i])/err[i])**2
    return chi

#period binning. Creates a list of period values with the given interval dT and maximum period value T_max
def period_list(dT, T_max):
    T_list=[]
    n=T_max/dT
    n = int(n)
    for i in range (1,n):
        T_list.append(dT*i)
    return T_list

#fractional part of x, to be used in saw function
def frac (x):
    return x - np.floor(x)


'''
'magnitude' is the objective function

define the function for calculating estimated magnitude values, depends on 3 varied parameters A-amplitude, C-phase shift, M_av-average magnitude

take x of the function to be (1/T)*t, where t is the time of observation and T the period
'''


def magnitude_sin (x, A, C, M_av):
    return A*np.sin((2*np.pi*x)+C)+M_av #sine must be in radians
    


'''
#saw tooth objective function
def magnitude_saw (x, A, C, M_av):
    if (C<0 or C>1 or A<0 ):
        ans = -1000*x
    else:
        ans = A*frac(x+C)+M_av
    return ans

'''
#modified saw tooth function, takes C(width)=1 for perfect sawtooth, unless modified otherwise (signal.sawtooth(x,C))
def magnitude_saw(x,A, M_av):
    return A*signal.sawtooth(x)+M_av


#folding function: takes in the list of observation times and the period. Any points past one period will be moved to an appropriate position between 0 and T
def fold (x_list, T):
    for i in range(len(x_list)):
        if x_list[i] > T:
            x_list[i]=np.mod(x_list[i],T)
    return x_list


#creates a list of magnitudes for given x with given parameter values
def magnitude_list_sin(x_list, A, C, M_av):
    m_mod = []
    for i in range(len(x_list)):
        m = magnitude_sin(x_list[i], A, M_av)
        m_mod.append(m)
    return m_mod

def magnitude_list_st(x_list, A, M_av):
    m_mod = []
    for i in range(len(x_list)):
        m = magnitude_saw(x_list[i], A, M_av)
        m_mod.append(m)
    return m_mod

'''
function that plots real data vs the model with optimal parameters for a given T. 
Usage:
t_list - list of observation times,
ydata - observed magnitudes,
sigma - errors on observed magnitudes,
n - the index in T list that corresponds to wanted period value
folding - 1 for plotting folded function (one period), else - unfolded
'''
def plot(t_list, ydata, sigma, n , folding, r, A_g_err):
    fold_list = list(t_list)    #creates a copy of t_list for folding
    popt = popt_list[n]
    T = T_list[n]
    pcov = pcov_list[n]
    #getting appropriate parameter errors from the covariance matrix - diagonal ements provide variance on the estimates
    if function =="st":
        A_err, M_av_err =np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])
    else:
        A_err, C_err, M_av_err =np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1]), (np.sqrt(pcov[2][2])+A_g_err)
    #if "folding parameter" equals 1, the function is folded, if not, fold_list=t_list
    if folding == 1:
        if function == "st":
            fold(fold_list,T*2*np.pi)
        else:
            fold(fold_list,T)

    if folding == 1:
        name = "folded"
    else:
        name = "unfolded"
   
    plt.errorbar(fold_list, ydata, yerr = sigma, fmt='o')
    #plt.text(0.1, (ydata.min()-0.3), '$A=$'+str("%.2f " % popt[0])+r'$\sigma_{A}$='+str("%.2f" % A_err) + '\n'+'$width=$'+str("%.2f " % popt[1])+r'$\sigma_{width}$='+str("%.2f" % C_err)+'\n'+'$M(av)=$'+str("%.2f " % popt[2])+r'$\sigma_{M(av)}$='+str("%.2f" % M_av_err))
    plt.text(0.1, (ydata.min()-0.3), '$A=$'+str("%.2f " % popt[0])+r'$\sigma_{A}$='+str("%.2f" % A_err) + '\n'+'$M(av)=$'+str("%.2f " % popt[1])+r'$\sigma_{M(av)}$='+str("%.2f" % M_av_err))


    
    plt.xlabel('Time (days)')
    plt.ylabel('Magnitude')
    #plt.title("Light curve for T=%.2f " % (T*2*np.pi) + names[r] )
    if function == "st":
        plt.title("Light curve for T=%.2f$\pm$%r " % (T*2*np.pi,errors[r]) + names[r] )
    else:
        plt.title("Light curve for T=%.2f$\pm$%r " % (T,errors[r]) + names[r] )
    plt.ylim(ydata.min()-0.4,ydata.max()+0.5)


    max_x = np.ceil(max(fold_list)/T)*T  #defining max x
    plt.xlim(0.,max_x)
    time_list = period_list(0.001,max_x)  #reusing period_list function to bin the time for model function x values
    x_arr = np.zeros(len(time_list))          #arrays to hold x and y values for model function
    m_arr = np.zeros(len(time_list))
    for i in range(len(time_list)):
        x_arr[i]=time_list[i]*(1/T)
        if function == "st":
            m_arr[i]=magnitude_saw(x_arr[i], *popt)
        else:
            m_arr[i]=magnitude_sin(x_arr[i], *popt)


    plt.xlabel("Time(days)")
    plt.ylabel("Magnitude")

    #plot light curve (magnitude vs time)
    plt.plot(time_list,m_arr)
    plt.savefig('%s_period_%s.png' %(names[r], name))
    plt.show()

    plt.clf()
    


'''
SET UP
'''

magnitudes = np.loadtxt('/Users/agnesemenaite/Documents/cepheids/data2/magnitudes.txt')
#star = np.int(sys.argv[1])
errors = np.loadtxt('/Users/agnesemenaite/Documents/cepheids/data2/errors.txt')




#observation time (days)
#t_list = [1., 6., 7., 11., 13., 16., 22., 25., 32., 35.]
t_list =[1,18,31,34,39,43,45,56,57]
#an empty numpy array to store the x values for 'magnitude' objective function, t_list = the list of observation times
xdata = np.zeros(len(t_list))
#measured luminosity magnitudes as numpy array
ydata = np.array(magnitudes)
#initial guess for parameter values for sine and saw tooth
x0sin = np.array([1.5, 0.1, 17.])
x0st = np.array([1.5, 17.])


#errors: the length of this array should be the same as the number of measurements (y data)
#sigma = np.array ([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
sigma = np.array(errors)
#names of stars
names = ["dl_cas","fm_cas", "rs_cas", "rw_cam", "rw_cas", "rx_cam", "ry_cas", "su_cas", "sw_cas", "tu_cas", "v636_cas", "v_lac", "vx_per", "andromeda"]
distances = [1688.,10.,1406.,1913.,3041., 882., 3828., 265., 2008., 808., 605., 1631., 2427.,10.]
E_BV = [0.53, 0.31, 0.88, 0.65, 0.42, 0.57, 0.65, 0.29, 0.49, 0.12, 0.7, 0.36, 0.51, 0.06]
errors = [0.09,0.06,0.07,0.4,0.03,0.06,0.1,0.004,0.03,0.003, 0.04,0.001,0.2]
#the objective function: sine for all the cepheids
function = sys.argv[1]
R_V= 3.86
R_V_err = 1.85
cg = 1.101349
cg_err = 0.124
chi0_list = []
T0_list = []
Mav_list = []
Merr_list = []


def periodfind(t_list, xdata,ydata, sigma,r, obf, T_s, E_BV):
    
    #correcting for dust extinction
    A_g = 4.251*E_BV
    A_g_err=0
    ydata = ydata-A_g
    
    #distance correction - converting to absolute magnitudes
    ydata = ydata-5*np.log10(distances[r]/10.)
    print ydata
  
    
    for i in range(len(sigma)):
        if sigma[i] < 0.05:
            sigma[i]=0.05
    
    '''
    T_list = period_list(0.005, 30.) #binning the period, must select a maximum value of period to be checked
    chi_list = []                   #empty list to hold chi square values
    popt_list = []                  #empty list to hold optimal parameters arrays
    chi_0 = 10000.                   #initialising chi_0 for comparison to find several smallest values. The initial value must be set such that it is bigger than any potential chi values from parameter fitting.
    pcov_list = []                  #empty list to hold covariance matrices
    T_0 = 0.                        #initialising T_0 to store a T value corresponding to chi_0. The initial value can be set to any value.
    
    for i in range(len(sigma)):
        if sigma[i] < 0.05:
            sigma[i]=0.05
    '''
            

    '''
        FINDING CHI SQ
        '''
    chi_0 = 10000.
    chi_1 = 10000. #initialising chi_0 for comparison to find several smallest values. The initial value must be set such that it is bigger than any potential chi values from parameter fitting.
    T_0 = 0.                        #initialising T_0 to store a T value corresponding to chi_0. The initial value can be set to any value.


    #this loops over binned periods
    for i in range(len(T_list)):
        T=T_list[i]
        for y in range(len(t_list)):
            xdata[y]=t_list[y]*(1/T)
        #optimization.curve_fit provides least squares fitting. Usage: (objective function, x values, y values, initial guess, data errors)
        #popt - array containing the optimal parameter values [amplitude, phase shift, average magnitude], pcov - 3x3s covariance matrix
        if  obf == "st":
       
            popt, pcov = optimize.curve_fit(magnitude_saw, xdata, ydata, x0st, sigma)

        else:
            popt, pcov = optimize.curve_fit(magnitude_sin, xdata, ydata, x0sin, sigma)
        

        #appending parameters and their errors in appropriate lists
        popt_list.append(popt)
        pcov_list.append(pcov)
        #print popt[0]
        #model magnitude for given T with optimal parameters (list of magnitudes for the given observing times)
        if obf == "st":
            m_mod = magnitude_list_st(xdata, *popt)
        else:
            m_mod = magnitude_list_sin(xdata, *popt)
        #chi squared value for given T appended to a list
        chi = chisq(ydata, m_mod, sigma)
        chi_list.append(chi)
      

        #check and store the values for three lowest chi estimates chi_0<chi_1<chi_2 and their corresponding T values
        if chi < chi_0:
            chi_1 = chi_0
            T_1 = T_0
            chi_0 = chi
            T_0 = T
        elif chi < chi_1:
            chi_2 = chi_1
            T_2 = T_1
            chi_1 = chi
            T_1 = T


    


    if obf =="st":
        plt.plot(np.array(T_list)*2*np.pi, chi_list)
    else:
        plt.plot(np.array(T_list), chi_list)
    plt.xlabel("Period(days)")
    plt.ylabel(r'$\chi^{2}$')
    plt.title("Least squares fitting "+names[r])
    #plt.savefig('%s_chiplot.png' %names[r])
    plt.show()
    plt.clf()
    
    
    
    #plot(t_list, ydata, sigma, np.int(T_s/(0.005*2*np.pi)), 0, r, A_g_err)
    #plot(t_list, ydata, sigma, np.int(T_s/(0.005*2*np.pi)), 1, r, A_g_err)

    
    plot(t_list, ydata, sigma, np.int(T_0/0.005), 0, r, A_g_err)
    plot(t_list, ydata, sigma, np.int(T_0/0.005), 1, r, A_g_err)
    
    chi0_list.append(chi_0)
    T0_list.append(T_0)
    Mav_list.append(popt[1])
    Merr_list.append(np.sqrt(pcov[1][1])+A_g_err)
    
    
    

    print names [r], chi_0, T_0, chi_1, T_1, chi_2, T_2





    
andromeda_m = np.array([17.8203963493, 17.93, 17.88, 18.188640793, 18.0989296827, 18.1443185715, 18.43, 18.512579682, 17.4829407938, 17.6461741271,17.608029682, 18.3042819049, 18.06, 18.4182696827, 17.42, 17.8867296827, 17.9484963493])
andromeda_t =[1., 2., 4., 6., 11., 12., 21., 26.,29., 30., 32., 42., 51., 56., 62., 97., 98]



for i in range(0,13):
    T_list = period_list(0.005, 30.) #binning the period, must select a maximum value of period to be checked
    chi_list = []                   #empty list to hold chi square values
    popt_list = []                  #empty list to hold optimal parameters arrays
    pcov_list = []                  #empty list to hold covariance matrices
    
    

    periodfind(t_list, xdata,ydata[i],sigma[i], i, function, 8.385, E_BV[i])


textinfo = [names[:-1], T0_list, Mav_list, Merr_list, chi0_list]
#writetxt(textinfo, "cepheids_prop_reduced_final.txt")
