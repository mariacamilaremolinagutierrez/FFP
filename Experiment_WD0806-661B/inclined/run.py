import numpy as np
import os, time, sys, signal
import ffp

def stop_code():
    print '\nSTOP'
    import sys
    sys.exit()

# Register an handler for the timeout
# From http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call
def handler(signum, frame):
    raise Exception("-10.0")

def permute(i_bp):

    #Fixed Values
    m0 = 0.58 #MSun
    m_ffp = 7.5 #MJupiter
    e_bp = 0.0

    n_r0_in_rinf = 40.0

    #Numbers of each parameter
    n_ms_bp = 3
    n_as_bp = 1
    n_bs_ffp = 1
    n_lans_bp = 10
    n_phis_bp = 1000
    
    ms_bp = [0.1, 2.822222, 2.822222]
    as_bp = [1.0, 28.222222, 33.666667]
    bs_ffp = [10.0, 25.0, 35.0]

    total_permutations = n_ms_bp*n_as_bp*n_bs_ffp*n_lans_bp*n_phis_bp

    #Variable parameters
    lans_bp = np.linspace(0.0,180.0,n_lans_bp) #degrees

    #Filenames to write info
    i_bp_filename = ('i%.4e')%(i_bp)
    filename_parameters = './particles/'+i_bp_filename+'/parameters_'+i_bp_filename+'.txt'
    filename_stables = './particles/'+i_bp_filename+'/stables_'+i_bp_filename+'.txt'

    #Counting the lines
    try:
        par = open('./particles/'+i_bp_filename+'/parameters_'+i_bp_filename+'.txt','r')
        num_runs = sum(1 for line in par)
        par.close()
    except:
        num_runs = 0
        os.makedirs('./particles/'+i_bp_filename)

    i=1

    for k in range(n_ms_bp):    
        
        m_bp = ms_bp[k]
        a_bp = as_bp[k]
        b_ffp = bs_ffp[k]
    
        for lan_bp in lans_bp:
    
            np.random.seed(12)
    
            for j in range(n_phis_bp):
    
                phi_bp = np.random.random()*360.0
    
                #Check that you don't run the last parameters combination again, it starts on the last one (in case the program has to be restarted)
                if (i>num_runs):
                    
                    #status
                    to_write = str(i)+'/'+str(int(total_permutations))+'\t'+str(m_bp)+'\t'+str(a_bp)+'\t'+str(b_ffp)+'\t'+str(phi_bp)+'\t'+str(i_bp)+'\t'+str(lan_bp)
                    res = os.system('echo "'+to_write+'" >> ./particles/'+i_bp_filename+'/status.txt')
    
                    #Start Handler
                    #Register the signal function handler
                    signal.signal(signal.SIGALRM, handler)
                    #Define a timeout for your function
                    signal.alarm(10)
                    
                    #Time starts
                    start_time = time.time()
    
                    try:
                        #Run it
                        t_end, max_energy_change, is_stable, e_star_ffp, e_star_bp, sma_star_ffp, sma_star_bp, inc_star_ffp, inc_star_bp, lan_star_ffp, lan_star_bp, ap_star_ffp, ap_star_bp = ffp.run_capture(m0_p=m0, m_ffp_p=m_ffp, e_bp_p=e_bp, m_bp_p=m_bp, a_bp_p=a_bp, b_ffp_p=b_ffp, phi_bp_p=phi_bp, inc_bp_p=i_bp, lan_bp_p=lan_bp, n_r0_in_rinf=n_r0_in_rinf)
                    
                        #Time stops
                        duration = time.time()-start_time
                        
                        #Write in files
                        line = ('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\t%f')%(m_bp,a_bp,b_ffp,phi_bp,i_bp,lan_bp,e_star_ffp,e_star_bp,sma_star_ffp,sma_star_bp,inc_star_ffp,inc_star_bp,lan_star_ffp,lan_star_bp,ap_star_ffp,ap_star_bp,duration,max_energy_change,t_end)
                        command = 'echo "'+line+'" >> '+filename_parameters
                        os.system(command)
                        
                        if(is_stable):
                            command = 'echo "'+line+'" >> '+filename_stables
                            os.system(command)
                    
                    #Raised when the program gets stucked -> inestable                            
                    except Exception, exc: 
                        #Write in files
                        line = ('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.4f\t%.4e\t%f')%(m_bp,a_bp,b_ffp,phi_bp,i_bp,lan_bp,10.0,10.0,-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,-10.0,0.0,0.0)
                        command = 'echo "'+line+'" >> '+filename_parameters
                        os.system(command)
                    
                    #Cancel the timer if the function returned before timeout
                    signal.alarm(0)
    
                #Advance counter
                i += 1

if __name__ in ('__main__', '__plot__'):

    #Obtain mass
    i_bp = float(sys.argv[1])
    
    
    #Permute
    permute(i_bp)
