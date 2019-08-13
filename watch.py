import os, time, subprocess

path_to_watch = '/var/kristoph_flask/outfiles'
path_to_watch2 = '/var/kristoph_flask/outfiles/testScripts'

before = dict ([(f, None) for f in os.listdir (path_to_watch)])
before2 = dict ([(z, None) for z in os.listdir (path_to_watch2)])

while 1:
    
    after = dict ([(f, None) for f in os.listdir (path_to_watch)])
    added = [f for f in after if not f in before]
    
    for i in added:
        if 'testScript_' in i:
            args1 = 'sudo mv /var/kristoph_flask/outfiles/'+i+' /var/kristoph_flask/outfiles/testScripts/'
            subprocess.run(args1, shell=True)
        elif '.fa' in i:
            args2 = 'sudo mv /var/kristoph_flask/outfiles/'+i+' /var/kristoph_flask/outfiles/fa_files/'
            subprocess.run(args2, shell=True)
        else:
            args3 = 'sudo rm /var/kristoph_flask/outfiles/'+i
            subprocess.run(args3, shell=True)
    
    after2 = dict ([(z, None) for z in os.listdir (path_to_watch2)])
    added2 = [z for z in after2 if not z in before2]
    
    if added2:
        for j in added2:
            print(j)
            time.sleep(5)
            args = '/var/www/vhosts/knaggert-vm.mdibl.net/flask_project/paHMM/bin/hmm < /var/kristoph_flask/outfiles/testScripts/'+j
            subprocess.run(args, shell=True)
    
    before = after
    before2 = after2
    added.clear()
    added2.clear()
