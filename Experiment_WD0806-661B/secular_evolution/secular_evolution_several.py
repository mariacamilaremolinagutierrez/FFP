import os

if __name__ in ('__main__', '__plot__'):

    #Masses Directories
    masses_directories = os.listdir('./particles/')
    
    os.system('mkdir secular_logs/')    
    
    for m_dir in masses_directories:
        command = '(amuse.sh secular_evolution.py '+m_dir+' >> secular_logs/log'+m_dir+'.txt 2>&1 &)'
        os.system(command)
